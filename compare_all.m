function results = compare_all()
    % ================================================================
    % compare_all.m - One-click comparison of NCB, CAFNCB, VCG-CAFNCB
    % ================================================================
    % All tunable parameters are centralized here for easy adjustment.
    % Do NOT delete the line "cfg = ncb_default_cfg(cfg);"

    rng(0);   % For reproducibility

    % ================== Basic Configuration ==================
    cfg.T  = 30;
    cfg.J  = 6;
    cfg.R  = 4;
    cfg.d  = 12;
    cfg.k  = 5;
    cfg.dy = 3;

    cfg.A = randn(cfg.dy, cfg.d);

    % Generate toy data
    Wtrue = zeros(cfg.d, cfg.T, cfg.R);
    for t=1:cfg.T
        for r=1:cfg.R
            base = sin(t/3) + 0.3*cos(r);
            Wtrue(:,t,r) = base * linspace(0.5,1.5,cfg.d)' + 0.1*randn(cfg.d,1);
        end
    end
    cfg.Yobs = zeros(cfg.dy, cfg.T, cfg.R);
    for t=1:cfg.T
        for r=1:cfg.R
            cfg.Yobs(:,t,r) = cfg.A * Wtrue(:,t,r) + 0.05*randn(cfg.dy,1);
        end
    end

    % Neighbors (ring topology)
    cfg.Nbrs = cell(cfg.R,1);
    for r=1:cfg.R
        cfg.Nbrs{r} = unique([mod(r-2,cfg.R)+1, mod(r,cfg.R)+1]);
    end

    % Edge weights initialization
    cfg.W_edge = cell(cfg.T, cfg.R);
    for t=1:cfg.T
        for r=1:cfg.R
            cfg.W_edge{t,r} = 0.5*randn(cfg.d, cfg.k);
        end
    end

    cfg.H = ones(1,cfg.d);
    cfg.proj_lambda = 1.0;

    cfg.delta1 = 2;
    cfg.delta2 = 2;
    cfg.beta   = 4.0;
    cfg.eta    = 0.25;
    
    % Call default configuration first
    cfg = ncb_default_cfg(cfg);

    % ================================================================
    % Tunable Parameters - Centralized Here (Recommended Balanced Version)
    % ================================================================

    cfg.macro_method        = "bfgs";    % Macro optimizer: "bfgs" is more stable for this setup

    % VCG Curvature Governance
    cfg.vcg_lambda          = 0.20;      % Braking strength (lower = gentler braking, allows more optimization)
    cfg.kappa_limit         = 5;      % Curvature threshold (higher = more tolerant to fluctuations)

    % Density Governance
    cfg.density_enable      = true;
    cfg.density_lambda      = 0.3;      % Density sensitivity (lower = less reactive to density changes)
    cfg.density_kappa_limit = 30;      % Density curvature threshold (higher = less frequent interventions)
    cfg.density_shrink      = 0.98;     % Shrink factor when density is too high (milder contraction)
    cfg.density_expand      = 1.02;     % Expand factor when density is too low (milder expansion)
    cfg.density_iter        = 1; 

    % Bandwidth & Normalization
    cfg.density_bandwidth   = 0.8;       % KDE bandwidth (0.8 is a balanced smoothing value)

    % Edge Length Normalization (for weights)
    cfg.min_edge_length     = 0.8;       % Minimum allowed edge length (prevents collapse to zero)
    cfg.max_edge_length     = 1.2;       % Maximum allowed edge length (prevents explosion)

    % Density Value Normalization (for rho)
    cfg.min_density         = 0.9;       % Minimum allowed density value
    cfg.max_density         = 1.1;       % Maximum allowed density value
    
cfg.eta_smooth = 0.6;
% One-step curvature smoothing coefficient.
% Defines how strongly the current curvature signal influences
% the smoothed governance signal.
% Lower values damp short-term fluctuations and reduce
% unintended Willow triggers.
    
    %cfg.density_enable = false;
    cfg.density_audit_stride = 3;% Set to 3 or another value depending on your preference
    cfg.debug_density = true;

    % ================== Run the three frameworks ==================
    % NCB
    cfg1 = cfg;
    cfg1.micro_counter = MicroCounter();
    cfg1.micro_counter.reset();
    out_ncb = ncb_run_generic(cfg1);

    % CAFNCB
    cfg2 = cfg;
    cfg2.micro_counter = MicroCounter();
    cfg2.micro_counter.reset();
    out_caf = cafncb_run_generic(cfg2);

    % VCG-CAFNCB
    cfg3 = cfg;
    cfg3.micro_counter = MicroCounter();
    cfg3.micro_counter.reset();
    out_vcg = vcg_cafncb_run_generic(cfg3);

    % ================== Display Results ==================
    results = struct();
    results.NCB_MSE  = out_ncb.metrics.final_mse;
    results.CAF_MSE  = out_caf.metrics.final_mse;
    results.VCG_MSE  = out_vcg.metrics.final_mse;

    results.NCB_qecr = out_ncb.metrics.gate_rate;
    results.CAF_qecr = out_caf.metrics.gate_rate;
    results.VCG_qecr = out_vcg.metrics.gate_rate;

  %  results.NCB_micro_success_rate = cfg1.micro_counter.succ / max(1, cfg1.micro_counter.att);
  %  results.CAF_micro_success_rate = cfg2.micro_counter.succ / max(1, cfg2.micro_counter.att);
  %  results.VCG_micro_success_rate = cfg3.micro_counter.succ / max(1, cfg3.micro_counter.att);

    results.VCG_willow_w   = out_vcg.metrics.willow_w_count;
    results.VCG_willow_rho = out_vcg.metrics.willow_rho_count;

    disp("=== COMPARISON (NCB vs CAFNCB vs VCG-CAFNCB) ===");
 %   fprintf("NCB:    MSE=%.4f, Gate=%.4f, MicroSucc=%.4f\n", results.NCB_MSE, results.NCB_gate, results.NCB_micro_success_rate);
 %   fprintf("CAF:    MSE=%.4f, Gate=%.4f, MicroSucc=%.4f\n", results.CAF_MSE, results.CAF_gate, results.CAF_micro_success_rate);
 %   fprintf("VCG:    MSE=%.4f, Gate=%.4f, MicroSucc=%.4f, Willow=%d+%d\n", ...
 %       results.VCG_MSE, results.VCG_gate, results.VCG_micro_success_rate, ...
 %       results.VCG_willow_w, results.VCG_willow_rho);
 
    fprintf("NCB:    MSE=%.4f, Qecr=%.4f\n", results.NCB_MSE, results.NCB_qecr);
    fprintf("CAF:    MSE=%.4f, Qecr=%.4f\n", results.CAF_MSE, results.CAF_qecr);
    fprintf("VCG:    MSE=%.4f, Qecr=%.4f, Willow=%d+%d\n", ...
           results.VCG_MSE, results.VCG_qecr,  ...
           results.VCG_willow_w, results.VCG_willow_rho);
end