function w_next = macro_step(t1,r,w,cfg,iter,maxIter)
    % ================================================================
    % 1. Eliminating Macro-Micro Tug-of-War (Most Important)
    % ================================================================
    % Reason: micro_fn was interfering with the macro optimization direction.
    % Solution: Remove premature micro correction here.
    % Micro QEC should only act as a "post-audit" in Step 3 (default_node_update).

    delta = cfg.delta_schedule_fn(cfg.delta0, iter, maxIter);
    w0 = w(:);                    % Use current center weight directly
    % w0 = cfg.micro_fn(w(:), delta, cfg);   % ← Commented out

    method = string(cfg.macro_method);

    switch method
        case "newton"
            g = cfg.grad_fn(t1,r,w0,cfg);
            H = cfg.hess_fn(t1,r,w0,cfg);
            H = (H+H')/2 + 1e-8*eye(numel(w0));
            p = - H \ g;
            alpha = armijo_backtracking(@(x) cfg.obj_fn(t1,r,x,cfg), w0, g, p, cfg.macro_max_ls);
            w_next = w0 + alpha*p;

        case "bfgs"
            persistent Binv_map;
            key = sprintf("t%d_r%d", t1, r);
            if isempty(Binv_map), Binv_map = containers.Map(); end
            n = numel(w0);
            if ~isKey(Binv_map,key) || iter == 1
                Binv_map(key) = eye(n);
            end
            Binv = Binv_map(key);

            g0 = cfg.grad_fn(t1,r,w0,cfg);
            p  = -Binv*g0;
            alpha = armijo_backtracking(@(x) cfg.obj_fn(t1,r,x,cfg), w0, g0, p, cfg.macro_max_ls);
            w1 = w0 + alpha*p;
            g1 = cfg.grad_fn(t1,r,w1,cfg);

            s = w1 - w0;
            y = g1 - g0;
            if (s'*y) > 1e-12
                rho = 1/(y'*s);
                I = eye(n);
                Binv = (I - rho*s*y')*Binv*(I - rho*y*s') + rho*(s*s');
                Binv_map(key) = Binv;
            end
            w_next = w1;

        otherwise % "sgd"
            g = cfg.grad_fn(t1,r,w0,cfg);
            w_next = w0 - cfg.macro_lr * g;
    end
end
