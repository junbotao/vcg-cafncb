function m = mean_mse_all(cfg, Wnode)
    T=cfg.T; R=cfg.R;
    acc = 0; n = 0;
    for t=1:T
        for r=1:R
            y = cfg.f_fn(t, Wnode(:,t,r), r, cfg);
            acc = acc + cfg.residual_fn(t,r,y,cfg);
            n = n + 1;
        end
    end
    m = acc / max(1,n);
end
