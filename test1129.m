initpath;

nSample = 50;
nP = 10;
maxSigma = 0.2;
maxOutlier = 5;

%% Make samples
samples = makeSamples(nSample, nP, maxSigma, maxOutlier);

% Divide graphs to train/test

%% Train graphs
parm.patterns = {samples.affinityMatrix};
parm.labels = {samples.GTbool};
parm.lossFn = @lossHamming;
parm.constraintFn = @constraint_MCMC;
parm.oracleFn = @constraint_MCMC;
parm.featureFn = @feature_mss;
parm.nP = nP;
parm.dimension = 3; % for feature_mss

args = ' -c 1 -o 2 -v 3 -e 0.00002 -# 30 -w 3';  % HARG:0.001, dw:0.0003, sw:0.0001
model = svm_struct_learn(args, parm);

%czesc
%abc

%anyoenghaseyo

%tests2
% conflict

<<<<<<< HEAD
=======
%Hi Yumin
>>>>>>> e90015eea38a601f8f02188757ca98b48fb88e68

%adding something


