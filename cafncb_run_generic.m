function out = cafncb_run_generic(cfg)
    cfg = ncb_default_cfg(cfg);

    T = cfg.T; J = cfg.J; R = cfg.R; d = cfg.d; k = cfg.k;

    out.W_center = zeros(d,T,R);
    out.W_node   = zeros(d,T,R);

    out.logs.gate_total = 0;
    out.logs.gate_hit   = 0;
    out.logs.mse_trace  = zeros(J,1);

    for t=1:T
        for r=1:R
            out.W_node(:,t,r) = cfg.project_fn(mean(cfg.W_edge{t,r},2));
        end
    end

    for j=1:J
        for t=1:T
            for r=1:R
                t1_list = window_t1(t,T,cfg.delta1);
                Wedge = cfg.W_edge{t,r};
                iter = 0;
                for t1 = t1_list
                    iter = iter + 1;
                    w = mean(Wedge,2);
                    w_next = macro_step(t1,r,w,cfg,iter,numel(t1_list));
                    Wedge = Wedge + repmat((w_next - w), 1, k);
                end
                cfg.W_edge{t,r} = Wedge;
            end

            Wcenter_guess = zeros(d,R);
            for r=1:R
                Wcenter_guess(:,r) = mean(cfg.W_edge{t,r},2);
            end

            for r=1:R
                nbrs = cfg.Nbrs{r};
                if isempty(nbrs)
                    w_fused = Wcenter_guess(:,r);
                else
                    rho = zeros(1,numel(nbrs));
                    for kk=1:numel(nbrs)
                        rr = nbrs(kk);
                        rho(kk) = cfg.corr_range_fn(t, r, rr, cfg, Wcenter_guess);
                    end
                    alpha = softmax_row(cfg.beta * rho);

                    w_mix = zeros(d,1);
                    for kk=1:numel(nbrs)
                        rr = nbrs(kk);
                        w_mix = w_mix + alpha(kk) * Wcenter_guess(:,rr);
                    end
                    w_fused = (1 - cfg.eta) * Wcenter_guess(:,r) + cfg.eta * w_mix;
                end
                out.W_center(:,t,r) = cfg.project_fn(w_fused);
            end

            for r=1:R
                t2_list = window_t2(t,T,cfg.delta2);
                for t2 = t2_list
                    y_old = cfg.f_fn(t2, out.W_node(:,t2,r), r, cfg);
                    w_assigned = out.W_center(:,t,r);
                    y_new = cfg.f_fn(t2, w_assigned, r, cfg);

                    r_old = cfg.residual_fn(t2,r,y_old,cfg);
                    r_new = cfg.residual_fn(t2,r,y_new,cfg);

                    out.logs.gate_total = out.logs.gate_total + 1;
                    if r_new < r_old
                        out.logs.gate_hit = out.logs.gate_hit + 1;
                        w_upd = cfg.node_update_fn(out.W_node(:,t2,r), w_assigned, t2, r, cfg);
                        out.W_node(:,t2,r) = cfg.project_fn(w_upd);
                    end
                end
            end
        end
        out.logs.mse_trace(j) = mean_mse_all(cfg, out.W_node);
    end

    out.metrics.final_mse = mean_mse_all(cfg, out.W_node);
    out.metrics.gate_rate = out.logs.gate_hit / max(1,out.logs.gate_total);
end
