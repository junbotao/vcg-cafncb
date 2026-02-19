function val = default_obj_mse(t1,r,w,cfg)
    yhat = cfg.f_fn(t1,w,r,cfg);
    y    = cfg.target_fn(t1,r,cfg);
    e = yhat(:) - y(:);
    val = mean(e.^2);
    if isnan(val), val = Inf; end
end
