function cfg = ncb_default_cfg(cfg)
% Fill defaults for NCB / CAFNCB / VCG-CAFNCB, runnable.
% Eq.11/Eq.12 are replaced by plug-in operators:
%   - macro_step() (Newton/BFGS/SGD + Armijo)
%   - node_update_fn() (gate-passed node update)
% Micro level uses cfg.micro_fn (wrapped) which counts coherence successes.

    must = {'T','J','R','d','k','Nbrs','W_edge'};
    for i=1:numel(must)
        if ~isfield(cfg, must{i})
            error("cfg.%s is required", must{i});
        end
    end

    % windows
    if ~isfield(cfg,'delta1'), cfg.delta1 = 2; end
    if ~isfield(cfg,'delta2'), cfg.delta2 = 2; end

    % CAF
    if ~isfield(cfg,'beta'), cfg.beta = 4.0; end
    if ~isfield(cfg,'eta'),  cfg.eta  = 0.25; end

    % 0th-order metric: projection Π_S
    if ~isfield(cfg,'H'), cfg.H = ones(1,cfg.d); end
    if ~isfield(cfg,'proj_lambda'), cfg.proj_lambda = 1.0; end
    if ~isfield(cfg,'project_fn')
        cfg.project_fn = @(w) project_linear_stabilizer(w, cfg.H, cfg.proj_lambda);
    end

    % model hooks
    if ~isfield(cfg,'f_fn')
        if isfield(cfg,'A')
            cfg.f_fn = @(t,w,r,cfg_) (cfg_.A * w(:));
        else
            cfg.f_fn = @(t,w,r,cfg_) w(:);
        end
    end

    if ~isfield(cfg,'target_fn')
        if isfield(cfg,'Yobs')
            cfg.target_fn = @(t,r,cfg_) cfg_.Yobs(:,t,r);
        else
            cfg.target_fn = @(t,r,cfg_) zeros(size(cfg.f_fn(t,zeros(cfg.d,1),r,cfg_)));
        end
    end

    if ~isfield(cfg,'residual_fn')
        cfg.residual_fn = @(t,r,yhat,cfg_) mse_residual_default(t,r,yhat,cfg_);
    end

    % macro optimizer (Eq.11 replacement)
    if ~isfield(cfg,'macro_method'), cfg.macro_method = "newton"; end % "newton"|"bfgs"|"sgd"
    if ~isfield(cfg,'macro_lr'), cfg.macro_lr = 0.05; end
    if ~isfield(cfg,'macro_max_ls'), cfg.macro_max_ls = 20; end

    if ~isfield(cfg,'obj_fn')
        cfg.obj_fn = @(t1,r,w,cfg_) default_obj_mse(t1,r,w,cfg_);
    end
    if ~isfield(cfg,'grad_fn')
        cfg.grad_fn = @(t1,r,w,cfg_) finite_diff_grad(@(x) cfg_.obj_fn(t1,r,x,cfg_), w);
    end
    if ~isfield(cfg,'hess_fn')
        cfg.hess_fn = @(t1,r,w,cfg_) finite_diff_hess(@(x) cfg_.obj_fn(t1,r,x,cfg_), w);
    end

    % micro level (counts coherence)
    if ~isfield(cfg,'delta0'), cfg.delta0 = cfg.delta1; end
    if ~isfield(cfg,'delta_schedule_fn')
        cfg.delta_schedule_fn = @(delta0,iter,maxIter) max(1, round(delta0*(1 - 0.8*iter/maxIter)));
    end

    if ~isfield(cfg,'micro_counter'), cfg.micro_counter = MicroCounter(); end

    if ~isfield(cfg,'micro_core_fn')
        cfg.micro_core_fn = @(w,delta,cfg_) cross_v(w,delta);  % core op (no counting)
    end

    % Success criterion: localSum decreases over the whole vector (coherence-style).
    if ~isfield(cfg,'micro_success_fn')
        cfg.micro_success_fn = @(w_before,w_after) (localSum(w_after,1:numel(w_after)) < localSum(w_before,1:numel(w_before)) - 1e-12);
    end

    cfg.micro_fn = @(w,delta,cfg_) micro_counting_wrapper(w,delta,cfg_);

    % node update (Eq.12 replacement)
    if ~isfield(cfg,'alpha_node'), cfg.alpha_node = 0.25; end
    if ~isfield(cfg,'node_update_fn')
        cfg.node_update_fn = @(w_old,w_assigned,t2,r,cfg_) default_node_update(w_old,w_assigned,t2,r,cfg_);
    end

    % VCG (curvature governance on increments)
    if ~isfield(cfg,'vcg_m'), cfg.vcg_m = 2; end
    if ~isfield(cfg,'vcg_lambda'), cfg.vcg_lambda = 0.6; end
    if ~isfield(cfg,'kappa_limit'), cfg.kappa_limit = 5.0; end
    if ~isfield(cfg,'rollback_gamma'), cfg.rollback_gamma = 0.7; end
    if ~isfield(cfg,'eps0'), cfg.eps0 = 1e-12; end

    % Density governance (optional)
    if ~isfield(cfg,'density_enable'), cfg.density_enable = true; end
    if ~isfield(cfg,'density_bandwidth'), cfg.density_bandwidth = []; end
    if ~isfield(cfg,'density_floor'), cfg.density_floor = 0.0; end
    if ~isfield(cfg,'density_ceil'),  cfg.density_ceil  = 1e12; end

    % ================== Edge Length Normalization (for weights) ==================
    if ~isfield(cfg,'min_edge_length'), cfg.min_edge_length = 0.5; end   % Minimum allowed edge length
    if ~isfield(cfg,'max_edge_length'), cfg.max_edge_length = 1.5; end   % Maximum allowed edge length

    % ================== Density Value Normalization (for rho) ==================
    if ~isfield(cfg,'min_density'), cfg.min_density = 0.5; end           % Minimum allowed density value
    if ~isfield(cfg,'max_density'), cfg.max_density = 1.5; end           % Maximum allowed density value

    if ~isfield(cfg,'density_lambda'), cfg.density_lambda = 0.8; end
    if ~isfield(cfg,'density_m'), cfg.density_m = 2; end
    if ~isfield(cfg,'density_kappa_limit'), cfg.density_kappa_limit = 6.0; end

    if ~isfield(cfg,'density_iter'), cfg.density_iter = 6; end
    if ~isfield(cfg,'density_shrink'), cfg.density_shrink = 0.92; end
    if ~isfield(cfg,'density_expand'), cfg.density_expand = 1.08; end
    if ~isfield(cfg,'density_scale_min'), cfg.density_scale_min = 0.2; end
    if ~isfield(cfg,'density_scale_max'), cfg.density_scale_max = 5.0; end

    % CAF corr function
    if ~isfield(cfg,'corr_range_fn')
        cfg.corr_range_fn = @(t,r,rr,cfg_,Wcenter_guess) corr_range_by_f(t,r,rr,cfg_,Wcenter_guess);
    end
end
