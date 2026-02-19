function alpha = armijo_backtracking(f, x, g, p, maxIter)
    alpha = 1.0;
    c = 1e-4;
    rho = 0.5;
    f0 = f(x);
    gp = g(:)'*p(:);
    for k=1:maxIter
        x_new = x + alpha*p;
        f_new = f(x_new);
        if f_new <= f0 + c*alpha*gp
            return;
        end
        alpha = rho*alpha;
        if alpha < 1e-12
            return;
        end
    end
end
