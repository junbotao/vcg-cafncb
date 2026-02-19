function out = vcg_cafncb_run_generic(cfg)
    % VCG-CAFNCB: Gen-4 Variational Curvature Governance Engine
    % Unifying Macro-scale Metric Reconstruction and Micro-scale QEC.
    % Based on Onsager Dissipative Principle and Rayleighian Minimization.

    cfg = ncb_default_cfg(cfg);

    T = cfg.T; J = cfg.J; R = cfg.R; d = cfg.d; k = cfg.k;

    % Output structures for both raw and governed states
    out.W_center = zeros(d,T,R);
    out.W_node   = zeros(d,T,R);
    out.W_center_governed = zeros(d,T,R);
    out.W_node_governed   = zeros(d,T,R);

    out.logs.gate_total = 0;
    out.logs.gate_hit   = 0;

    % Governance hit logs (The "Willow" Mechanism index)
    out.logs.willow_hit_w   = false(J,T,R);
    out.logs.willow_hit_rho = false(J,T,R);

% ===== log initialization =====
out.logs.kappa = zeros(J,T,R);

if ~isfield(out.logs,'kappa_smooth')
    out.logs.kappa_smooth = zeros(J,T,R);
end
    
    % Step 0: Initialize Manifold via Projection
    for t=1:T
        for r=1:R
            out.W_node(:,t,r) = cfg.project_fn(mean(cfg.W_edge{t,r},2));
        end
    end
    out.W_node_governed = out.W_node;

    % Historical state tracking for Curvature (kappa) extraction
    Wc_prev  = out.W_center;   Wc_prev2 = out.W_center;
    Wn_prev  = out.W_node_governed; Wn_prev2 = out.W_node_governed;

    % Density tracking for Scale-Flow Governance
    rhoC_prev  = zeros(T,R); rhoC_prev2 = zeros(T,R);
    rhoN_prev  = zeros(T,R); rhoN_prev2 = zeros(T,R);
    rhoDWC_prev = zeros(T,R); rhoDWC_prev2 = zeros(T,R);
    rhoDWN_prev = zeros(T,R); rhoDWN_prev2 = zeros(T,R);

    for t=1:T
        for r=1:R
            rhoC_prev(t,r) = density_scalar(Wc_prev(:,t,r), cfg.density_bandwidth);
            rhoN_prev(t,r) = density_scalar(Wn_prev(:,t,r), cfg.density_bandwidth);
        end
    end
    rhoC_prev2 = rhoC_prev; rhoN_prev2 = rhoN_prev;

    
% --- Curvature smoothing buffer (for VCG) ---
if ~isfield(out.logs,'kappa_smooth')
    out.logs.kappa_smooth = zeros(J, T, R);
end

% --- One-step lag smoothing coefficient for kappa ---
if isfield(cfg,'eta_smooth')
    eta_smooth = cfg.eta_smooth;
else
    eta_smooth = 0.7; % default
end
    
    % Main Evolution Loop
    for j=1:J
        for t=1:T
            % --- Part A: Macro-step Propagation (Newton/BFGS) ---
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

            % --- Part B: Cross-Area Fusion (CAF) via Coherence Correlation ---
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
                        w_mix = w_mix + alpha(kk)*Wcenter_guess(:,rr);
                    end
                    w_fused = (1-cfg.eta)*Wcenter_guess(:,r) + cfg.eta*w_mix;
                end
                out.W_center(:,t,r) = cfg.project_fn(w_fused);
            end

            % --- Part C: Residual Gating (SEC Audit) ---
            for r=1:R
                t2_list = window_t2(t,T,cfg.delta2);
                for t2 = t2_list
                    y_old = cfg.f_fn(t2, out.W_node_governed(:,t2,r), r, cfg);
                    w_assigned = out.W_center(:,t,r);
                    y_new = cfg.f_fn(t2, w_assigned, r, cfg);

                    r_old = cfg.residual_fn(t2,r,y_old,cfg);
                    r_new = cfg.residual_fn(t2,r,y_new,cfg);

                    out.logs.gate_total = out.logs.gate_total + 1;
                    if r_new < r_old
                        out.logs.gate_hit = out.logs.gate_hit + 1;
                        w_upd = cfg.node_update_fn(out.W_node_governed(:,t2,r), w_assigned, t2, r, cfg);
                        out.W_node_governed(:,t2,r) = cfg.project_fn(w_upd);
                    end
                end
            end

            % --- Part D: Gen-4 Variational Curvature Governance (VCG) ---
            for r=1:R
                % 1. Extract Curvature Gradient (kappa) for State and Node
                wC = out.W_center(:,t,r);
                dw1 = wC - Wc_prev(:,t,r);
                dw2 = Wc_prev(:,t,r) - Wc_prev2(:,t,r);
                kappa = norm(dw1 - dw2, 2) / (norm(dw2, 2) + cfg.eps0+ 1e-3);
                
% --- kappa smoothing uses previous outer-iteration history (j-1) ---
if j > 1
    kappa_prev = out.logs.kappa_smooth(j-1, t, r);
else
    kappa_prev = 0;
end

kappa_smooth = (1 - eta_smooth) * kappa_prev + eta_smooth * kappa;
out.logs.kappa_smooth(j, t, r) = kappa_smooth;
                
% ------------------------------------------------------------
% Small-step guard:
% If the reference step magnitude is too small, the curvature
% ratio becomes noise-dominated and should not trigger governance.
% ------------------------------------------------------------
if norm(dw2,2) < 1e-2 || norm(dw1-dw2,2) < 1e-6
    kappa = 0;
end



                % ================================================================
                % 3. VCG Governance Signal Relaxation (Smoothing)
                % ================================================================
                % Reason: Instant kappa is too sensitive and causes false Willow triggers.
                % Solution: Add 1st-order lag smoothing to filter numerical noise.
                % This allows us to safely lower kappa_limit back to 5~8.

% ================================================================
% 3. VCG Governance Signal Relaxation (Smoothing)
% ================================================================
% Reason: instantaneous kappa can be noise-dominated near convergence.
% Solution: 1st-order lag smoothing using previous outer iteration (j-1).

if ~isfield(out.logs,'kappa_smooth')
    out.logs.kappa_smooth = zeros(J,T,R);
end

if isfield(cfg,'eta_smooth')
    eta_smooth = cfg.eta_smooth;
else
    eta_smooth = 0.7;  % default fallback
end

% Use (j-1) history, not (j)
if j > 1
    kappa_prev = out.logs.kappa_smooth(j-1,t,r);
else
    kappa_prev = 0;
end

kappa_smooth = (1-eta_smooth)*kappa_prev + eta_smooth*kappa;
out.logs.kappa_smooth(j,t,r) = kappa_smooth;


                % Use smoothed kappa for governance
                A_k = sqrt(1 + kappa_smooth.^cfg.vcg_m);
                alpha_step = 1 / (1 + cfg.vcg_lambda * A_k);
                wC_g = Wc_prev(:,t,r) + alpha_step * (wC - Wc_prev(:,t,r));

                wN = out.W_node_governed(:,t,r);
                dnw1 = wN - Wn_prev(:,t,r);
                dnw2 = Wn_prev(:,t,r) - Wn_prev2(:,t,r);
                kappaN = norm(dnw1 - dnw2, 2) / (norm(dnw2, 2) + cfg.eps0);

                A_kN = sqrt(1 + kappaN.^cfg.vcg_m);
                alpha_stepN = 1 / (1 + cfg.vcg_lambda * A_kN);
                wN_g = Wn_prev(:,t,r) + alpha_stepN * (wN - Wn_prev(:,t,r));

                % 2. Execute Willow Phase-Transition (Emergency Rollback)
if kappa_smooth > cfg.kappa_limit || kappaN > cfg.kappa_limit
out.logs.willow_hit_w(j,t,r) = true;

g = cfg.rollback_gamma;  % Rollback strength in [0,1]. Larger g => stronger rollback.

% Roll back current state toward the last stable reference (previous outer-iteration state)
wC_g = g * Wc_prev(:,t,r) + (1-g) * wC_g;
wN_g = g * Wn_prev(:,t,r) + (1-g) * wN_g;
end

                % 3. Density-driven Scale-Flow Governance
                if cfg.density_enable
                    [wC_d, rhoC_now, ~, kC, ~, hitC] = density_govern_state(wC_g, rhoC_prev(t,r), rhoC_prev2(t,r), cfg);
                    [wN_d, rhoN_now, ~, kN, ~, hitN] = density_govern_state(wN_g, rhoN_prev(t,r), rhoN_prev2(t,r), cfg);

                    dwC = wC_d - Wc_prev(:,t,r);
                    dwN = wN_d - Wn_prev(:,t,r);

                    [dwC_d, rhoDWC_now, ~, kDWC, ~, hitDWC] = density_govern_increment(dwC, rhoDWC_prev(t,r), rhoDWC_prev2(t,r), cfg);
                    [dwN_d, rhoDWN_now, ~, kDWN, ~, hitDWN] = density_govern_increment(dwN, rhoDWN_prev(t,r), rhoDWN_prev2(t,r), cfg);

                    wC_d2 = Wc_prev(:,t,r) + dwC_d;
                    wN_d2 = Wn_prev(:,t,r) + dwN_d;

                    % Scale-flow Willow Mechanism
                    if hitC || hitN || hitDWC || hitDWN || (kC>cfg.density_kappa_limit) || (kN>cfg.density_kappa_limit)
                        out.logs.willow_hit_rho(j,t,r) = true;
                        g = cfg.rollback_gamma;
                        wC_d2 = (1-g)*wC_d2 + g*Wc_prev(:,t,r);
                        wN_d2 = (1-g)*wN_d2 + g*Wn_prev(:,t,r);
                    end
                    
                    % Final State Assignment after Scale-flow
                    wC_final = wC_d2;
                    wN_final = wN_d2;
                    
                    % Update Density History
                    rhoC_prev2(t,r) = rhoC_prev(t,r); rhoC_prev(t,r) = rhoC_now;
                    rhoN_prev2(t,r) = rhoN_prev(t,r); rhoN_prev(t,r) = rhoN_now;
                    rhoDWC_prev2(t,r) = rhoDWC_prev(t,r); rhoDWC_prev(t,r) = rhoDWC_now;
                    rhoDWN_prev2(t,r) = rhoDWN_prev(t,r); rhoDWN_prev(t,r) = rhoDWN_now;
                else
                    wC_final = wC_g;
                    wN_final = wN_g;
                end

                % --- [CRITICAL AUDIT FIX]: Unified Micro-Macro Convergence ---
                % Re-injecting Micro-scale QEC after Governance to refine residual ripples.
                % This ensures micro_counter stats (success_rate) are active and MSE is minimized.
                delta_vcg = cfg.delta_schedule_fn(cfg.delta0, j, J);
                
                out.W_center_governed(:,t,r) = cfg.micro_fn(wC_final, delta_vcg, cfg);
                out.W_node_governed(:,t,r)   = cfg.micro_fn(wN_final, delta_vcg, cfg);
                
                % Enforce Metric Manifold Constraints (Stabilizer Projection)
                out.W_center_governed(:,t,r) = cfg.project_fn(out.W_center_governed(:,t,r));
                out.W_node_governed(:,t,r)   = cfg.project_fn(out.W_node_governed(:,t,r));
            end
        end

        % Update Time-step Evolution History
        Wc_prev2 = Wc_prev;  Wc_prev = out.W_center_governed;
        Wn_prev2 = Wn_prev;  Wn_prev = out.W_node_governed;
    

  
%     disp([j, kappa, kappa_smooth, norm(dw2,2)]);
  
    end

    % Final Performance Audit
    out.metrics.final_mse = mean_mse_all(cfg, out.W_node_governed);
    out.metrics.gate_rate = out.logs.gate_hit / max(1,out.logs.gate_total);
    out.metrics.willow_w_count   = nnz(out.logs.willow_hit_w);
    out.metrics.willow_rho_count = nnz(out.logs.willow_hit_rho);
end
