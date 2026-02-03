function results = micro_macro_moment_demo()
    rng(1234);

    settings = struct();
    settings.T = 120;
    settings.N = 400;
    settings.pi = [0.55 0.45];
    settings.simul_drop = 20;
    settings.mod_name = 'micro_macro_two_type';
    settings.base_dir = fileparts(mfilename('fullpath'));
    settings.param_guess = struct('beta', [0.96 0.93], 'sigma', [1.2 1.4], ...
        'h', [0.5 0.35], 'omega', [0.6 0.45], 'barF', [0.12 0.2], ...
        'chi', [0.7 0.3], 'rho_y', [0.9 0.85], 'sigma_eps', [0.3 0.45]);

    true_params = struct('beta', [0.97 0.94], 'sigma', [1.2 1.4], ...
        'h', [0.55 0.35], 'omega', [0.65 0.45], 'barF', [0.12 0.22], ...
        'chi', [0.75 0.25], 'rho_y', [0.9 0.85], 'sigma_eps', [0.3 0.45]);

    data = simulate_dynare_micro(true_params, settings);
    targets = compute_moments(data, settings);

    theta0 = pack_params(settings.param_guess);
    objective = @(theta) moment_objective(theta, settings, targets);
    options = optimset('Display', 'iter', 'TolX', 1e-4, 'TolFun', 1e-4, ...
        'MaxIter', 300, 'MaxFunEvals', 1500);
    [theta_hat, fval] = fminsearch(objective, theta0, options);

    params_hat = unpack_params(theta_hat);
    model_hat = simulate_dynare_micro(params_hat, settings);
    moments_hat = compute_moments(model_hat, settings);

    results = struct();
    results.true_params = true_params;
    results.estimated_params = params_hat;
    results.targets = targets;
    results.moments_hat = moments_hat;
    results.objective = fval;

    disp(format_param_table(true_params, params_hat));
    disp(table(targets(:), moments_hat(:), 'VariableNames', {'Target', 'Model'}));
end

function data = simulate_dynare_micro(params, settings)
    dynare_dir = fullfile(settings.base_dir, 'dynare');
    dynare_file = fullfile(dynare_dir, settings.mod_name);
    write_dynare_params(params, settings, dynare_dir);
    addpath(settings.base_dir, dynare_dir);
    dynare(dynare_file, 'noclearall', 'nolog');

    oo_ = evalin('base', 'oo_');
    M_ = evalin('base', 'M_');
    vmap = containers.Map(M_.endo_names, 1:numel(M_.endo_names));

    C1 = extract_series(oo_, vmap, 'C1', settings);
    C2 = extract_series(oo_, vmap, 'C2', settings);
    F1 = extract_series(oo_, vmap, 'F1', settings);
    F2 = extract_series(oo_, vmap, 'F2', settings);
    Y1 = extract_series(oo_, vmap, 'Y1', settings);
    Y2 = extract_series(oo_, vmap, 'Y2', settings);
    pF = extract_series(oo_, vmap, 'pF', settings);

    T = numel(C1);
    N = settings.N;
    data = struct();
    data.C = zeros(T, N, 2);
    data.F = zeros(T, N, 2);
    data.Y = zeros(T, N, 2);
    data.participation = false(T, N, 2);
    data.pF = pF;
    data.Cg = [C1 C2];
    data.Fg = [F1 F2];

    for g = 1:2
        eps_income = params.sigma_eps(g) .* randn(T, N);
        Yg = (g == 1) * Y1 + (g == 2) * Y2;
        Yi = Yg .* exp(eps_income);
        Ci = (g == 1) * C1 + (g == 2) * C2;
        Fi = (g == 1) * F1 + (g == 2) * F2;
        Ci = Ci .* (Yi ./ Yg);
        Fi = Fi .* (Yi ./ Yg);

        pr = logistic(-0.6 + params.chi(g) .* log(Yi));
        participation = rand(T, N) < pr;

        data.Y(:, :, g) = Yi;
        data.C(:, :, g) = Ci;
        data.F(:, :, g) = Fi;
        data.participation(:, :, g) = participation;
    end
end

function series = extract_series(oo_, vmap, name, settings)
    idx = vmap(name);
    raw = oo_.endo_simul(idx, :)';
    raw = raw(1 + settings.simul_drop:end);
    series = raw(1:settings.T);
end

function write_dynare_params(params, settings, dynare_dir)
    param_file = fullfile(dynare_dir, 'micro_macro_params.m');
    fid = fopen(param_file, 'w');
    fprintf(fid, 'pi1 = %.6f;\n', settings.pi(1));
    fprintf(fid, 'pi2 = %.6f;\n', settings.pi(2));
    fprintf(fid, 'beta1 = %.6f;\n', params.beta(1));
    fprintf(fid, 'beta2 = %.6f;\n', params.beta(2));
    fprintf(fid, 'sigma1 = %.6f;\n', params.sigma(1));
    fprintf(fid, 'sigma2 = %.6f;\n', params.sigma(2));
    fprintf(fid, 'h1 = %.6f;\n', params.h(1));
    fprintf(fid, 'h2 = %.6f;\n', params.h(2));
    fprintf(fid, 'omega1 = %.6f;\n', params.omega(1));
    fprintf(fid, 'omega2 = %.6f;\n', params.omega(2));
    fprintf(fid, 'barF1 = %.6f;\n', params.barF(1));
    fprintf(fid, 'barF2 = %.6f;\n', params.barF(2));
    fprintf(fid, 'rho_y1 = %.6f;\n', params.rho_y(1));
    fprintf(fid, 'rho_y2 = %.6f;\n', params.rho_y(2));
    fprintf(fid, 'Rbar = 1.01;\n');
    fprintf(fid, 'pFbar = 1.0;\n');
    fprintf(fid, 'rho_r = 0.7;\n');
    fprintf(fid, 'rho_pf = 0.8;\n');
    fprintf(fid, 'sigma_r = 0.01;\n');
    fprintf(fid, 'sigma_pf = 0.02;\n');
    fprintf(fid, 'sigma_y1 = 0.02;\n');
    fprintf(fid, 'sigma_y2 = 0.02;\n');
    fprintf(fid, 'Y1bar = 1.0;\n');
    fprintf(fid, 'Y2bar = 0.8;\n');
    fclose(fid);
end

function moments = compute_moments(data, settings)
    moments = [];
    for g = 1:2
        C = data.C(:, :, g);
        Y = data.Y(:, :, g);
        F = data.F(:, :, g);
        pF = data.pF;
        food_share = (pF .* F) ./ (C + pF .* F);
        log_cir = log(C ./ Y);
        part_rate = mean(data.participation(:, :, g), 'all');
        moments_g = [mean(food_share, 'all'); mean(log_cir, 'all'); var(log_cir, 0, 'all'); part_rate];
        moments = [moments; moments_g];
    end

    C_agg = settings.pi(1) .* data.Cg(:, 1) + settings.pi(2) .* data.Cg(:, 2);
    F_agg = settings.pi(1) .* data.Fg(:, 1) + settings.pi(2) .* data.Fg(:, 2);
    pF = data.pF;
    food_share_agg = (pF .* F_agg) ./ (C_agg + pF .* F_agg);
    c_growth = diff(log(C_agg));
    moments = [moments; mean(c_growth); mean(food_share_agg)];
end

function value = moment_objective(theta, settings, targets)
    params = unpack_params(theta);
    data = simulate_dynare_micro(params, settings);
    model_moments = compute_moments(data, settings);
    weights = diag([repmat([10 5 2 10], 1, 2) 5 5]);
    diff = targets - model_moments;
    value = diff' * weights * diff;
end

function theta = pack_params(params)
    theta = [logit(params.beta(1)); logit(params.beta(2)); log(params.sigma(1)); log(params.sigma(2)); ...
        logit(params.h(1)); logit(params.h(2)); logit(params.omega(1)); logit(params.omega(2)); ...
        log(params.barF(1)); log(params.barF(2)); logit(params.chi(1)); logit(params.chi(2)); ...
        logit(params.rho_y(1)); logit(params.rho_y(2)); log(params.sigma_eps(1)); log(params.sigma_eps(2))];
end

function params = unpack_params(theta)
    params = struct();
    params.beta = [logistic(theta(1)); logistic(theta(2))];
    params.sigma = [exp(theta(3)); exp(theta(4))];
    params.h = [logistic(theta(5)); logistic(theta(6))];
    params.omega = [logistic(theta(7)); logistic(theta(8))];
    params.barF = [exp(theta(9)); exp(theta(10))];
    params.chi = [logistic(theta(11)); logistic(theta(12))];
    params.rho_y = [logistic(theta(13)); logistic(theta(14))];
    params.sigma_eps = [exp(theta(15)); exp(theta(16))];
end

function y = logistic(x)
    y = 1 ./ (1 + exp(-x));
end

function y = logit(x)
    x = min(max(x, 1e-6), 1 - 1e-6);
    y = log(x ./ (1 - x));
end

function tbl = format_param_table(true_params, params_hat)
    names = {'beta_1', 'beta_2', 'sigma_1', 'sigma_2', 'h_1', 'h_2', 'omega_1', 'omega_2', ...
        'barF_1', 'barF_2', 'chi_1', 'chi_2', 'rho_y1', 'rho_y2', 'sigma_eps1', 'sigma_eps2'};
    true_values = [true_params.beta(1); true_params.beta(2); true_params.sigma(1); true_params.sigma(2); ...
        true_params.h(1); true_params.h(2); true_params.omega(1); true_params.omega(2); ...
        true_params.barF(1); true_params.barF(2); true_params.chi(1); true_params.chi(2); ...
        true_params.rho_y(1); true_params.rho_y(2); true_params.sigma_eps(1); true_params.sigma_eps(2)];
    est_values = [params_hat.beta(1); params_hat.beta(2); params_hat.sigma(1); params_hat.sigma(2); ...
        params_hat.h(1); params_hat.h(2); params_hat.omega(1); params_hat.omega(2); ...
        params_hat.barF(1); params_hat.barF(2); params_hat.chi(1); params_hat.chi(2); ...
        params_hat.rho_y(1); params_hat.rho_y(2); params_hat.sigma_eps(1); params_hat.sigma_eps(2)];
    tbl = table(true_values, est_values, 'RowNames', names, 'VariableNames', {'True', 'Estimated'});
end
