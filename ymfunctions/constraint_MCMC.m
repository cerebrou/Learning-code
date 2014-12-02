function yhat = constraint_MCMC(param, model, M, y)

    % warning: only for same number of outliers on both graphs
    nMatch = length(M);
    nP = sqrt(nMatch);
    E12 = ones(nP,nP);
    [L12(:,1), L12(:,2)] = find(E12);
    [group1, group2] = make_group12(L12);

    lambda = model.w;
    yhat = MCMC_aug(M, group1, group2, E12,L12, lambda, y);
    yhat = yhat(:);
end