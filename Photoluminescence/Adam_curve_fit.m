function [params, loss_history] = Adam_curve_fit(model_func, x_data, y_data, init_params, options)
    % ADAM_CURVE_FIT Fits a function to data using the Adam optimizer, with bounds and noise injection.
    %
    % Inputs:
    %   - model_func: function handle, y_hat = model_func(params, x_data)
    %   - x_data: input data
    %   - y_data: target/output data
    %   - init_params: initial guess (column vector)
    %   - options: struct with optional fields:
    %       .lr - learning rate (default: 0.01)
    %       .beta1 - decay rate for 1st moment (default: 0.9)
    %       .beta2 - decay rate for 2nd moment (default: 0.999)
    %       .epsilon - small constant (default: 1e-8)
    %       .max_iters - number of iterations (default: 1000)
    %       .loss_type - 'mse' or a custom loss function handle
    %       .lb - lower bounds for parameters (default: -Inf)
    %       .ub - upper bounds for parameters (default: +Inf)
    %
    % Outputs:
    %   - params: optimized parameter vector
    %   - loss_history: loss at each iteration
    %
    % [params, loss_history] = adam_curve_fit(model_func, x_data, y_data, init_params, options)
    % includes periodic random noise to avoid local minima.

    % Default options
    if ~isfield(options, 'lr'), options.lr = 0.01; end
    if ~isfield(options, 'beta1'), options.beta1 = 0.9; end
    if ~isfield(options, 'beta2'), options.beta2 = 0.999; end
    if ~isfield(options, 'epsilon'), options.epsilon = 1e-8; end
    if ~isfield(options, 'max_iters'), options.max_iters = 1000; end
    if ~isfield(options, 'loss_type'), options.loss_type = 'mse'; end
    if ~isfield(options, 'lb'), options.lb = -Inf(size(init_params)); end
    if ~isfield(options, 'ub'), options.ub = Inf(size(init_params)); end

    % Set up loss function
    if isa(options.loss_type, 'function_handle')
        loss_func = options.loss_type;
    elseif strcmpi(options.loss_type, 'mse')
        loss_func = @(yhat, y) mean((yhat - y).^2);
    else
        error('Unsupported loss type');
    end

    % Initialize
    params = init_params(:);
    m = zeros(size(params));
    v = zeros(size(params));
    loss_history = zeros(options.max_iters, 1);
    delta = 1e-5; % finite difference step size

    for t = 1:options.max_iters
        % Compute prediction and loss
        y_hat = model_func(params, x_data);
        loss_history(t) = loss_func(y_hat, y_data);

        % Estimate gradient numerically
        grad = zeros(size(params));
        for i = 1:length(params)
            dp = zeros(size(params));
            dp(i) = delta;
            y1 = model_func(params + dp, x_data);
            y2 = model_func(params - dp, x_data);
            grad(i) = (loss_func(y1, y_data) - loss_func(y2, y_data)) / (2 * delta);
        end

        % Adam updates
        m = options.beta1 * m + (1 - options.beta1) * grad;
        v = options.beta2 * v + (1 - options.beta2) * (grad .^ 2);
        m_hat = m / (1 - options.beta1^t);
        v_hat = v / (1 - options.beta2^t);

        % Parameter update
        params = params - options.lr * m_hat ./ (sqrt(v_hat) + options.epsilon);

        % Inject random noise periodically to help escape local minima
        if mod(t, 100) == 0
            noise = randn(size(params)) * 0.01;
            params = params + noise;
        end

        % Apply parameter bounds (clipping)
        params = max(min(params, options.ub(:)), options.lb(:));
    end
end
