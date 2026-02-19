function a = softmax_row(x)
    x = x(:)';
    m = max(x);
    ex = exp(x - m);
    s  = sum(ex);
    if s <= 0 || isnan(s)
        a = ones(size(x))/numel(x);
    else
        a = ex/s;
    end
end
