function [k, g_weights] = getBroadenKer(x_vals, f_vals, F_handle)
    % INVERT_G Solve f(x) = ∫ F(k,x) g(k) dk for g(k)
    % Inputs:
    %   f_vals   - vector of f(x) values at points x_vals
    %   x_vals   - corresponding x values
    %   F_handle - function handle: F(k, x_vector)
    %
    % Outputs:
    %   k         - symmetric log grid
    %   g_weights - estimated values of g(k)

    %% Parameters
    N = 5e2;
    k_min = 1e-4;
    k_max = 1e-2;

    % Symmetric log grid centered at 0
    k_pos = logspace(log10(k_min), log10(k_max), N);
    k = [-fliplr(k_pos), k_pos];  % size: 1 x (2N)

    % Trapezoidal weights
    w = zeros(1, length(k));
    w(2:end-1) = 0.5 * (k(3:end) - k(1:end-2));
    w(1) = k(2) - k(1);
    w(end) = k(end) - k(end-1);
    w = abs(w);

    %% Build matrix A(i,j) = F(k_j, x_i) * w_j
    A = zeros(length(x_vals), length(k));
    for j = 1:length(k)
        A(:, j) = F_handle(k(j), x_vals(:)) * w(j);
    end

    %% Initial guess for g: normalized Gaussian
    g0 = exp(-k.^2 / (2 * 0.001^2));
    g0 = g0 / (g0 * w');

    %% Objective function
    % Regularization: penalize deviation from target standard deviation
    sigma0 = 1e-3;     % desired width
    lambda = 1e3;     % strength of the penalty
    
    % Define the full objective
    % Objective = data misfit + std-deviation regularization
    objective = @(g) norm(A * g(:) - f_vals(:))^2 + ...
                     lambda * (sum(g(:) .* (k(:).^2) .* w(:)) - sigma0^2);





    %% Normalization constraint
    norm_constraint = @(g) deal([], g * w' - 1);  % no inequality, one equality

    %% Lower bound for g(k): g(k) >= 0
    lb = zeros(size(g0));  % lower bound (g(k) >= 0)

    %% Custom output function to monitor convergence
    convergence_metric = @(x, optimValues, state) monitor_convergence(x, optimValues, state);

    %% Solve with box constraints and convergence monitoring
    opts = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                        'OutputFcn', convergence_metric);
    [g_weights, ~, exitflag] = fmincon(objective, g0, [], [], [], [], lb, [], norm_constraint, opts);

    if exitflag <= 0
        warning('Optimization did not converge.');
    end

    %% Plot
    figure;
    plot(k, g_weights, 'LineWidth', 2);
    xlabel('k'); ylabel('g(k)');
    title('Recovered g(k) with Non-Negativity Constraint');
    grid on;
end

% Output function to monitor convergence
function stop = monitor_convergence(~, optimValues, state)
    stop = false;  % continue optimization by default
    
    % Check if it's the first iteration
    if strcmp(state, 'init')
        disp('Starting optimization...');
    end
    
    % Display the objective function value at each iteration
    if strcmp(state, 'iter')
        fprintf('Iteration: %d, Objective Value: %.4e, Rel. Change: %.4e\n', ...
            optimValues.iteration, optimValues.fval, optimValues.stepsize);
    end
end
