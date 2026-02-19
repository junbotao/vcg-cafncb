function idx = window_t2(t,T,delta2)
    lo = max(1, round(t - delta2/1.5));
    hi = min(T, round(t + delta2));
    idx = lo:hi;
end
