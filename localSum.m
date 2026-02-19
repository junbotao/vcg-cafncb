function s = localSum(w, idx)
    % ================================================================
    % 2. Correcting Micro Success Criterion
    % ================================================================
    % Original: sum(w.^2) tried to drive energy to zero.
    % Problem: This conflicted with the main MSE objective.
    % New: Use variance (structure consistency) instead.
    % This aligns with "maintain structural consistency without injecting 
    % non-physical energy".

    w_local = w(idx);
    s = sum((w_local - mean(w_local)).^2);   % Variance proxy
end