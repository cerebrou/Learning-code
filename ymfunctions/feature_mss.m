function psi = feature_mss(param,M,y)
    psi = zeros(3,1);
    psi(1) = y(:)'*M*y(:);
    psi(2) = -sum(y);
    psi(3) = -sum(y)^2;
    psi = sparse(psi);
end