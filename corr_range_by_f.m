function rho = corr_range_by_f(t,r,rr,cfg,Wcenter_guess)
    T = cfg.T;
    win = window_t2(t,T,cfg.delta1);
    fr = [];
    fi = [];
    for tt = win
        yr  = cfg.f_fn(tt, Wcenter_guess(:,r),  r,  cfg);
        yi  = cfg.f_fn(tt, Wcenter_guess(:,rr), rr, cfg);
        fr  = [fr; yr(:)]; %#ok<AGROW>
        fi  = [fi; yi(:)]; %#ok<AGROW>
    end

    if std(fr) < 1e-12 || std(fi) < 1e-12
        rho = 0.0;
    else
        C = corrcoef(fr,fi);
        rho = C(1,2);
        if isnan(rho), rho = 0.0; end
    end
end
