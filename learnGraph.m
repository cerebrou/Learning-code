function model = learnGraph(graphs, matches, nP, dFeat, refGraph, learningType)

% if(strcmp(learningType,'sw-marius')
%     model = ;
%     return model;
% end

parm.patterns = graphs;
parm.labels = matches;
parm.lossFn = @lossCB;
parm.nP = nP;
parm.dFeat = dFeat;
parm.gm = refGraph;

options = [];
options.lambda = 1;
options.gap_threshold = 0.1; % duality gap stopping criterion
options.num_passes = 100; % max number of passes through data
options.do_line_search = 1;
options.debug = 0; % for displaying more info (makes code about 3x slower)

switch(learningType)
    case 'non'
        model = [];
        return;
    case 'HARG'
        parm.constraintFn = @constraintCB;
        parm.oracleFn = @constraintCB;
        parm.featureFn = @featureCB;
        parm.dimension = dFeat * nP^2;
    case 'dw'
        parm.constraintFn = @constraintCB_dw;
        parm.oracleFn = @constraintCB_dw;
        parm.featureFn = @featureCB_dw;
        parm.dimension = nP^2;
    case 'sw'
        parm.constraintFn = @constraintCB_sw;
        parm.oracleFn = @constraintCB_sw;
        parm.featureFn = @featureCB_sw;
        parm.dimension = 2;
    case 'sw-marius'
        Ms = makeAffSet(graphs,refGraph,dFeat);
        model.w = marius_sw(Ms,parm);
        return;
    case 'dw-marius'
        Ms = makeAffSet(graphs,refGraph,dFeat);
%         model.w = marius_sw(Ms,parm);
        model.w = marius_dw(Ms,parm);
        return;
end

% args = ' -c 1 -o 2 -v 3 -e 0.00002 -# 30 -w 3';  % HARG:0.001, dw:0.0003, sw:0.0001
% model = svm_struct_learn(args, parm);
% [model, progress] = solverBCFW(parm, options);
[model, progress] = solverFW(parm, options);
% [model, progress] = solverSSG(parm, options);

end

function psi = featureCB(param, g, y)
%     disp('in featureCB');
    gperm = fitGraphToModel(g,y);
    psi = graphToPhi(gperm);
%     psi = sparse(double(psi));
end

function psi = featureCB_dw(param, g, y)
%     disp('in featureCB');
    % load gmodel*
    gm = param.gm;
    gperm = fitGraphToModel(g,y);
    wV = [];
%     wV = sum(gperm.fV .* gmodel.fV,1);
    wE = sum(gperm.fE .* gm.fE,1);
    psi = [wV(:); wE(:)];
%     psi = sparse(double(psi));
end

function psi = featureCB_sw(param, g, y)
%     disp('in featureCB');
    % load gmodel*
    gm = param.gm;
    gperm = fitGraphToModel(g,y);
    wV = [];
%     wV = sum(gperm.fV .* gmodel.fV,2);
    wE = sum(gperm.fE .* gm.fE,2);
    wE1 = sum(wE(1:13));
    wE2 = sum(wE(14:end));
    psi = [wV(:); wE1; wE2];
%     psi = sparse(double(psi));
end

function delta = lossCB(param, y, ybar)
    delta = double(sum(y~=ybar)) / length(y);
end

% function yhat = constraintCB(param, model, g, y)
function yhat = constraintCB(param, model, g, varargin)
    if(length(varargin) > 3)
        y = varargin{1};
    else
        y = [];
    end
%     figure, plot(model.w);
%     disp('in constraintCB');
    %
    gm = phiToGraph(model.w, param.nP);
    affinityMatrix = makeAff(gm, g, y);
    nP = gm.nNode;
    nQ = g.nNode;
    %
    E12 = ones(nP,nQ);
    [L12(:,1) L12(:,2)] = find(E12);
    [group1 group2] = make_group12(L12);
    %
    yraw = RRWM(double(affinityMatrix), group1, group2);
%     yraw = SM(double(affinityMatrix));
%     [binaryAss,~] = hungarian(reshape(yraw,[nP,nQ]));
    ysol = greedyMapping(yraw, group1, group2);
    ysol = reshape(ysol, [nP,nQ]);
    [yhat, ~] = find(ysol');
end

% function yhat = constraintCB_dw(param, model, g, y)
function yhat = constraintCB_dw(param, model, g, varargin)
    if(length(varargin) > 3)
        y = varargin{1};
    else
        y = [];
    end
%     figure, plot(model.w);
%     disp('in constraintCB');
    %
    % load gmodel*
    gm = param.gm;
    gm.fE = gm.fE .* repmat(model.w(:)',[param.dFeat,1]);
    affinityMatrix = makeAff(gm, g, y);
    nP = gm.nNode;
    nQ = g.nNode;
    %
    E12 = ones(nP,nQ);
    [L12(:,1) L12(:,2)] = find(E12);
    [group1 group2] = make_group12(L12);
    %
    yraw = RRWM(double(affinityMatrix), group1, group2);
%     yraw = SM(double(affinityMatrix));
%     [binaryAss,~] = hungarian(reshape(yraw,[nP,nQ]));
    ysol = greedyMapping(yraw, group1, group2);
    ysol = reshape(ysol, [nP,nQ]);
    [yhat, ~] = find(ysol');
end

% function yhat = constraintCB_sw(param, model, g, y)
function yhat = constraintCB_sw(param, model, g, varargin)
    if(length(varargin) > 3)
        y = varargin{1};
    else
        y = [];
    end
%     figure, plot(model.w);
%     disp('in constraintCB');
    %
    % load gmodel*
    gm = param.gm;
    gm.fE(1:13,:) = gm.fE(1:13,:) .* repmat(model.w(1),[13,param.nP^2]);
    gm.fE(14:end,:) = gm.fE(14:end,:) .* repmat(model.w(2),[18,param.nP^2]);
    affinityMatrix = makeAff(gm, g, y);
    nP = gm.nNode;
    nQ = g.nNode;
    %
    E12 = ones(nP,nQ);
    [L12(:,1) L12(:,2)] = find(E12);
    [group1 group2] = make_group12(L12);
    %
    yraw = RRWM(double(affinityMatrix), group1, group2);
%     yraw = SM(double(affinityMatrix));
%     [binaryAss,~] = hungarian(reshape(yraw,[nP,nQ]));
    ysol = greedyMapping(yraw, group1, group2);
    ysol = reshape(ysol, [nP,nQ]);
    [yhat, ~] = find(ysol');
end