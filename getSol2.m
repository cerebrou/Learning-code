function sol = getSol2(v, labels, nodes)

%minRatio = 0.2;

minRatio = eps; %0.00000000000001;

[maxConfidence, bestAssig] = max(v);

minConfidence = maxConfidence * minRatio;

n = length(v);

sol = zeros(n,1);

while 1

    [m, ind] = max(v);

    if m <= minConfidence
        break;
    end

    sol(ind) = 1;

    f = find(nodes == nodes(ind));

    v(f) = 0;

    f = find(labels == labels(ind));

    v(f) = 0;

end

sol = sol/norm(sol);

return