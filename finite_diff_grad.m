function g = finite_diff_grad(f, x)
    x = x(:);
    n = numel(x);
    g = zeros(n,1);
    fx = f(x);
    h = 1e-6;
    for i=1:n
        xp = x; xp(i)=xp(i)+h;
        g(i) = (f(xp)-fx)/h;
    end
    g(isnan(g)) = 0;
end
