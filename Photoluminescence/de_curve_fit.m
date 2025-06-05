function [best_params, loss_history] = de_curve_fit(model_func, x_data, y_data, init_params, options)
    % DE_CURVE_FIT Fits a function using Differential Evolution (global optimizer)
    % Inputs and outputs match adam_curve_fit.

    % Parameters and bounds
    D = length(init_params);
    if ~isfield(options, 'lb'), options.lb = -10*ones(D,1); end
    if ~isfield(options, 'ub'), options.ub = 10*ones(D,1); end
    if ~isfield(options, 'pop_size'), options.pop_size = 20*D; end
    if ~isfield(options, 'max_iters'), options.max_iters = 1000; end
    if ~isfield(options, 'F'), options.F = 0.8; end % differential weight
    if ~isfield(options, 'CR'), options.CR = 0.9; end % crossover rate
    if ~isfield(options, 'loss_type'), options.loss_type = 'mse'; end

    % Loss function
    if isa(options.loss_type, 'function_handle')
        loss_func = options.loss_type;
    elseif strcmpi(options.loss_type, 'mse')
        loss_func = @(yhat, y) mean((yhat - y).^2);
    else
        error('Unsupported loss type');
    end

    % Initialize population
    pop = rand(options.pop_size, D) .* (options.ub(:)' - options.lb(:)') + options.lb(:)';
    fitness = zeros(options.pop_size, 1);
    for i = 1:options.pop_size
        y_hat = model_func(pop(i,:)', x_data);
        fitness(i) = loss_func(y_hat, y_data);
    end

    loss_history = zeros(options.max_iters, 1);

    % Main loop
    for iter = 1:options.max_iters
        for i = 1:options.pop_size
            % Mutation: select 3 random indices
            idxs = randperm(options.pop_size, 3);
            while any(idxs == i)
                idxs = randperm(options.pop_size, 3);
            end
            a = pop(idxs(1), :);
            b = pop(idxs(2), :);
            c = pop(idxs(3), :);

            % Mutation and crossover
            mutant = a + options.F * (b - c);
            cross = rand(1, D) < options.CR;
            if ~any(cross), cross(randi(D)) = true; end
            trial = pop(i, :);
            trial(cross) = mutant(cross);

            % Clip to bounds
            trial = max(min(trial, options.ub(:)'), options.lb(:)');

            % Selection
            y_trial = model_func(trial', x_data);
            trial_fitness = loss_func(y_trial, y_data);
            if trial_fitness < fitness(i)
                pop(i, :) = trial;
                fitness(i) = trial_fitness;
            end
        end

        % Record best loss
        loss_history(iter) = min(fitness);
    end

    % Return best
    [~, best_idx] = min(fitness);
    best_params = pop(best_idx, :)';
end
