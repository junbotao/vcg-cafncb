function H = finite_diff_hess(f, x)
    x = x(:);
    n = numel(x);
    H = zeros(n,n);
    h = 1e-4;
    g0 = finite_diff_grad(f, x);
    for i=1:n
        xi = x; xi(i)=xi(i)+h;
        gi = finite_diff_grad(f, xi);
        H(:,i) = (gi - g0)/h;
    end
    H = (H+H')/2;
    H(isnan(H)) = 0;
end
