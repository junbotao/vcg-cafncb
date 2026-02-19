function w2 = micro_counting_wrapper(w,delta,cfg)
% Wrapper around cfg.micro_core_fn that increments coherence counters.

    w1 = w(:);
    w2 = cfg.micro_core_fn(w1,delta,cfg);

    cfg.micro_counter.att = cfg.micro_counter.att + 1;

    ok = false;
    try
        ok = cfg.micro_success_fn(w1,w2);
    catch
        ok = (norm(w2-w1) > 1e-12);
    end
    if ok
        cfg.micro_counter.succ = cfg.micro_counter.succ + 1;
    end
end
