function samples = makeSamples(nSample,nP,maxSigma,maxOutlier)

sigma = floor(maxSigma*rand);
nOut = floor(maxOutlier*rand);

gparam.sigma = sigma;
gparam.nOut = nOut;
gparam.ratioFill = 1;
gparam.nInlier = nP;
gparam.scale_2D = 0.1;

for iSample = 1 : nSample
    samples(iSample) = makeSample(gparam);
end

end

function sample = makeSample(gparam)

deformation = gparam.sigma;
nOut = gparam.nOut;
ratioFill = gparam.ratioFill;
nInlier = gparam.nInlier;
scale_2D = gparam.scale_2D;

nP1 = nInlier + nOut;
nP2 = nInlier + nOut;
seq = randperm(nP2);

G1 = tril(rand(nP1),-1); % lower triangular graph
G1 = G1+G1';
P = tril(rand(nP1),-1);
P = P+P';
P = P > ratioFill;
G1(P) = NaN;

N = deformation*tril(randn(nP2),-1);
N = N+N';

G2 = tril(rand(nP2),-1);
G2 = G2+G2';
P = tril(rand(nP2),-1);
P = P+P';
P = P > ratioFill;
G2(P) = NaN;
G2(seq(1:nInlier),seq(1:nInlier)) = G1(1:nInlier,1:nInlier);
G2 = G2+N;

E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);

M = (repmat(G1, nP2, nP2)-kron(G2,ones(nP1))).^2;
M = exp(-M./scale_2D);
% M = max(2-M,0);
M(isnan(M)) = 0;
M(1:(length(M)+1):end)=0;
M = M+M';

GT.seq = seq(1:nInlier);
GT.matrix = zeros(nP1, nP2);
for i = 1:nInlier, GT.matrix(i,seq(i)) = 1; end
GT.bool = GT.matrix(:);

%% Return results
sample.nP1 = nP1;
sample.nP2 = nP2;
sample.G1 = G1;
sample.G2 = G2;

sample.E12 = E12;
sample.L12 = L12;
sample.affinityMatrix = M;
% sample.productMatrix = M2;
sample.group1 = group1;
sample.group2 = group2;
sample.GTbool = GT.bool;


end