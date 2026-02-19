function w_new = default_node_update(w_old,w_assigned,t2,r,cfg)
    w_new = (1 - cfg.alpha_node)*w_old(:) + cfg.alpha_node*w_assigned(:);
    delta = cfg.delta_schedule_fn(cfg.delta0, 1, 1);
    w_new = cfg.micro_fn(w_new, delta, cfg);
end
