function idx = window_t1(t,T,delta1)
    lo = max(1, t-delta1);
    hi = min(T, t+delta1);
    idx = lo:hi;
end
