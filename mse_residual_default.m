function val = mse_residual_default(t,r,yhat,cfg)
    y = cfg.target_fn(t,r,cfg);
    e = yhat(:) - y(:);
    val = mean(e.^2);
    if isnan(val), val = Inf; end
end
