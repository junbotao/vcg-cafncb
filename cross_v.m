function wH = cross_v(wH,delta)
    wH = wH(:);
    n  = size(wH,1);
    for i = 1:n
        idx   = max(1, i-delta) : min(n, i+delta);
        Ri    = localSum(wH, idx);

        bestR = Ri; bestJ = i;
        for j = idx
            if j==i, continue; end
            wTemp = wH; wTemp(i) = wH(j);
            Rj = localSum(wTemp, idx);
            if Rj < bestR
                bestR = Rj; bestJ = j;
            end
        end

        if bestJ ~= i
            wH(i) = wH(i) + (wH(bestJ)-wH(i))*rand;
        end
    end
end
