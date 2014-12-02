function psi = graphToPsi(g)
% INPUT
% g.fV: (dV by nV) matrix
% g.fE: (dE by nE) matrix
% g.edges: (2 by nE) matrix
% g.nNodes
%
% OUTPUT
% psi: (dFeat by 1) vector, where dFeat = dV + dE

psi = [g.fV(:); g.fE(:)];

end

