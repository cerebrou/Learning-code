
%author: Marius Leordeanu
%date:   Feb 05, 2008

%computes the approximate derivative of the principal eigenvector of M with
%respect to a specific parameter (indirectly known through dM, which is the
%derivative of M with respect to that parameter). The approximation is
%based on the power method of computing the eigenvector. 

%Input:
%M         - the original matrix
%dM        - the derivative of M w.r.t the desired parameter
%nSteps    - number of steps used by the Power Method to compute the
%            principal eigenvector (recommended value: between 30 and 100 
             % for the spectral matching algorithm)


function [v, dv] = getEigenDerivative(M, dM, nSteps)

% disp(' computing derivative of principal eigenvector ');

% tic;

n = size(M,1);

v = ones(n,1);

v = M*v;
c = norm(v);
v = v/c;

dMn =  (dM*ones(n,1))/c;

for i = 2:nSteps
   
    dMn = dM*v + M*dMn; 
    
    v = M*v;
    c = norm(v);
    v = v/c;

    dMn = dMn/c;
   
end

dv = dMn - v*(v'*dMn);

% toc;

end