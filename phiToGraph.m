function g = psiToGraph(psi,nP)
% INPUT
% psi: (dFeat by 1) vector, where dFeat = dV + dE
%
% OUTPUT
% g.fV: (dV by nV) matrix
% g.fE: (dE by nE) matrix
% g.edges: (2 by nE) matrix
% g.nNodes

% g.fV = psi(1:nP);
% g.fE = psi(nP+1:end);

% disp('in psiToGraph');

dFeat = 31;
g.fV = [];
g.fE = reshape(psi,[dFeat, nP^2]);
g.edges = [repmat((1:nP)',nP,1), kron((1:nP)',ones(nP,1))];
g.nNode = nP;

end