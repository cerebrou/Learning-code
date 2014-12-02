function fE = getFeatE(points, i, j)
    attE = points(j,:) - points(i,:);
    fE = makeFeatBin(attE);
    fE = fE(:);
end