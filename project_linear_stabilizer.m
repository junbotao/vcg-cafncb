function w_proj = project_linear_stabilizer(w, H, lambda)
    if nargin < 3, lambda = 1.0; end
    w = w(:);
    n = numel(w);
    A = speye(n) + lambda*(H'*H);
    w_proj = A \ w;
    if any(isnan(w_proj)), w_proj = w; end
end
