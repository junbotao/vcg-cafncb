function [w_out, rho_now, rho_target, kappa_rho, alpha_rho, hit] = density_govern_state(w_in, rho_prev, rho_prev2, cfg)
% ------------------------------------------------------------
% Cold-start guard (state-level):
% If density history is uninitialized, skip 2nd-diff curvature audit.
% ------------------------------------------------------------
if abs(rho_prev) < 1e-12 && abs(rho_prev2) < 1e-12
    w_out = w_in;
    rho_now = density_scalar(w_in(:), cfg.density_bandwidth);
    rho_target = rho_now;
    kappa_rho = 0;
    alpha_rho = 0;
    hit = false;
    return;
end
    w = w_in(:);
    rho_now = density_scalar(w, cfg.density_bandwidth);
    rho_now = min(max(rho_now, cfg.density_floor), cfg.density_ceil);

    dr1 = rho_now - rho_prev;
    dr2 = rho_prev - rho_prev2;
    kappa_rho = abs(dr1 - dr2) / (abs(dr2) + cfg.eps0+ 1e-3);

    A = sqrt(1 + (kappa_rho.^cfg.density_m));
    alpha_rho = 1 / (1 + cfg.density_lambda * A);
    rho_target = rho_prev + alpha_rho * (rho_now - rho_prev);
    
% ------------------------------------------------------------
% HIT gating: only apply state density scaling when curvature exceeds limit
% ------------------------------------------------------------
hit = (kappa_rho > cfg.density_kappa_limit);

if ~hit
    w_out = w_in;
    return;
end

% ---------- Only executed if hit == true ----------
mu = mean(w);
v = w - mu;

scale = 1.0;
w_tmp = w;

for it=1:cfg.density_iter
    rho_tmp = density_scalar(w_tmp, cfg.density_bandwidth);
    if rho_tmp < rho_target
        scale = max(cfg.density_scale_min, scale*cfg.density_shrink);
    else
        scale = min(cfg.density_scale_max, scale*cfg.density_expand);
    end
    w_tmp = mu + scale*v;
end

w_out = w_tmp;

% Enforce limits
w_out = max(w_out, cfg.min_edge_length);
w_out = min(w_out, cfg.max_edge_length);

rho_now = max(rho_now, cfg.min_density);
rho_now = min(rho_now, cfg.max_density);

end
