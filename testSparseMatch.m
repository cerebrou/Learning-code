footpath = 'e:/codes/fgm-master';
addpath(genpath([footpath '/src']));
addpath(genpath([footpath '/lib']));


clear variables;
prSet(1);

%% src parameter
tag = 'house';
nIn = [30 30]; % randomly remove 2 nodes
parKnl = st('alg', 'cmum'); % type of affinity: only edge distance

for iFs = 1 : 100
pFs = [1 iFs]; % frame index

%% algorithm parameter
[pars, algs] = gmPar(2);

%% src
wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);
asgT = wsSrc.asgT;

%% feature
parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
wsFeat = cmumAsgFeat(wsSrc, parG, parF, 'svL', 1);
[gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

%% affinity
[nV1,nE1] = size(gphs{1}.G);
[nV2,nE2] = size(gphs{2}.G);
[KP, KQ] = conKnlGphPQU(gphs, parKnl);
% KP = zeros(nV1,nV2);
% fE1 = makeFeatBin([gphs{1}.dsts; gphs{1}.angs]);
% fE2 = makeFeatBin([gphs{2}.dsts; gphs{2}.angs]);
% KQ = double(fE1' * fE2);
% KQ = KQ(1:nE1, 1:nE2);
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));

E12 = ones(nV1,nV2);
[L12(:,1) L12(:,2)] = find(E12);
[group1 group2] = make_group12(L12);
yraw = RRWM(double(K), group1, group2);
ysol = greedyMapping(yraw, group1, group2);
ysol = reshape(ysol, [nV1,nV2]);
[y, ~] = find(ysol');

[i,j] = find(asgT.X');
ytrue = zeros(nV1,1);
for ii = 1 : length(j)
    ytrue(j(ii)) = i(ii);
end
accdist(iFs) = 1 - sum(ytrue ~=  y) / nV1;

end

figure, plot(accdist,'linewidth',2,'color','r');
ylim([0 1.05]);
xlabel('frame');
ylabel('accuracy');
hold on, plot(acc,'linewidth',2);
legend('diff','dot');