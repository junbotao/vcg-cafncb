function [dw_out, rho_dw_now, rho_dw_target, kappa_dw, alpha_dw, hit] = density_govern_increment(dw_in, rho_prev, rho_prev2, cfg)
% ------------------------------------------------------------
% Cold-start guard (increment-level):
% No valid history for 2nd-difference audit.
% If both previous densities are ~0, skip increment governance.
% ------------------------------------------------------------
if abs(rho_prev) < 1e-12 && abs(rho_prev2) < 1e-12
    dw_out = dw_in;
    rho_dw_now = density_scalar(dw_in(:), cfg.density_bandwidth);
    rho_dw_target = rho_dw_now;
    kappa_dw = 0;
    alpha_dw = 0;
    hit = false;
    return;
end
    dw = dw_in(:);
    rho_dw_now = density_scalar(dw, cfg.density_bandwidth);
    rho_dw_now = min(max(rho_dw_now, cfg.density_floor), cfg.density_ceil);

    dr1 = rho_dw_now - rho_prev;
    dr2 = rho_prev - rho_prev2;
    kappa_dw = abs(dr1 - dr2) / (abs(dr2) + cfg.eps0+ 1e-3);

    A = sqrt(1 + (kappa_dw.^cfg.density_m));
    alpha_dw = 1 / (1 + cfg.density_lambda * A);
    rho_dw_target = rho_prev + alpha_dw * (rho_dw_now - rho_prev);
    
% ------------------------------------------------------------
% HIT gating: only apply density scaling when curvature exceeds limit
% ------------------------------------------------------------
hit = (kappa_dw > cfg.density_kappa_limit);

if ~hit
    dw_out = dw_in;
    return;
end

% ---------- Only executed if hit == true ----------
mu = mean(dw);
v = dw - mu;

scale = 1.0;
dw_tmp = dw;

for it=1:cfg.density_iter
    rho_tmp = density_scalar(dw_tmp, cfg.density_bandwidth);
    if rho_tmp < rho_dw_target
        scale = max(cfg.density_scale_min, scale*cfg.density_shrink);
    else
        scale = min(cfg.density_scale_max, scale*cfg.density_expand);
    end
    dw_tmp = mu + scale*v;
end

dw_out = dw_tmp;

% Enforce edge limits
dw_out = max(dw_out, cfg.min_edge_length);
dw_out = min(dw_out, cfg.max_edge_length);

%     % Enforce density value limits (for rho)
%     rho_dw_now = max(rho_dw_now, cfg.min_density);
%     rho_dw_now = min(rho_dw_now, cfg.max_density);
%     if isfield(cfg,'debug_density') && cfg.debug_density
%     fprintf('[INC] kappa_dw=%.4g  hit=%d  alpha_dw=%.4g  rho_now=%.4g  rho_prev=%.4g  rho_prev2=%.4g\n', ...
%         kappa_dw, hit, alpha_dw, rho_dw_now, rho_prev, rho_prev2);
%     end

end
