function model = yumin_ssvm_cutting(input)
% Cutting plane learning for graph matching
%
% INPUT
% X: (nTrain by 1) cell
% X.xy: (dFeat, nFeat) matrix
% Y: (nTrain by 1) cell
%
% OUTPUT
% model

x = input.patterns;
y = input.labels;
parm.dimension = input.dimension;
parm.dFeat = input.dFeat;
parm.nP = input.nP;
parm.gm = input.gm;
nTrain = length(x);

itermax = 100;
thres = 0.5;

workingSet = [];
for iter = 1 : itermax
    workingSetChanged = false;
    z = solvecvx(parm,workingSet,x,y);
    for iTrain = 1 : nTrain
        yhat = constraintCB_dw(parm,z(1:parm.dimension),x{iTrain},y{iTrain});
%         yhat = constraintCB(parm,z(1:parm.dimension),x{iTrain},y{iTrain});
        if(lossCB(parm,yhat,y{iTrain})>thres)
            if(isempty(workingSet) || ~ismember(iTrain,[workingSet.iTrain]) ...
                    || min(lossCBs(parm,workingSet([workingSet.iTrain]==iTrain),yhat)) ~= 0)
                tmp.iTrain = iTrain;
                tmp.yhat = yhat;
                workingSet = [workingSet,tmp];
                workingSetChanged = true;
            end
        end
    end
    if(~workingSetChanged)
        break;
    end
end

model.w = z;

end


function z = solvecvx(parm,workingSet, x, y)
    % one slack
    nActive = length(workingSet);
    nTrain = length(x);
    dFeat = parm.dFeat;
    nP = parm.nP;
    dimension = parm.dimension;
    H = zeros(dimension+nTrain);
    H(1:dimension,1:dimension) = eye(dimension);
    C = 30;
%     C = 1/4;
    
    cvx_begin
        variables z(dimension+nTrain,1)
%         minimize( z'*H*z + C*[zeros(1,dimension), ones(1,nTrain)]*z )
        minimize( norm(z(1:dimension),1) + C/nTrain*ones(1,nTrain)*z((dimension+1):end) )
%         minimize( sum(norms(reshape(z(1:dimension),dFeat,nP^2),2)) + C/nTrain*ones(1,nTrain)*z((dimension+1):end) )
%         minimize( sum(norms(reshape(z(1:dimension),dFeat,nP^2),2)) + C/nTrain*ones(1,nTrain)*z((dimension+1):end) )
        subject to
            z((dimension+1):end) >= 0;
            for iActive = 1 : nActive
                iTrain = workingSet(iActive).iTrain;
                yhat = workingSet(iActive).yhat;
                z(1:dimension)'*(featureCB_dw(parm, x{iTrain},y{iTrain})-featureCB_dw(parm, x{iTrain},yhat))  >= lossCB(parm,y{iTrain},yhat) - z(dimension+iTrain);
%                 z(1:dimension)'*(featureCB(parm, x{iTrain},y{iTrain})-featureCB(parm, x{iTrain},yhat))  >= lossCB(parm,y{iTrain},yhat) - z(dimension+iTrain);
            end
    cvx_end
    
    % n-slack
    % cvx/TFOCS/quadprog
end

function delta = lossCB(param, y, ybar)
    delta = double(sum(y~=ybar)) / length(y);
end

function deltas = lossCBs(param, y, ybar)
    n = length(y);
    deltas = zeros(n,1);
    for i = 1 : n
        deltas(i) = double(sum(y(i).yhat~=ybar)) / length(y);
    end
end

function psi = featureCB(param, g, y)
%     disp('in featureCB');
    gperm = fitGraphToModel(g,y);
    psi = graphToPhi(gperm);
    psi = sparse(double(psi));
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
    psi = sparse(double(psi));
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
    psi = sparse(double(psi));
end


function yhat = constraintCB(param, w, g, y)
% function yhat = constraintCB(param, w, g, varargin)
%     if(length(varargin) > 3)
%         y = varargin{1};
%     else
%         y = [];
%     end
%     figure, plot(model.w);
%     disp('in constraintCB');
    %
    gm = phiToGraph(w, param.nP);
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

function yhat = constraintCB_dw(param, w, g, y)
% function yhat = constraintCB_dw(param, w, g, varargin)
%     if(length(varargin) > 3)
%         y = varargin{1};
%     else
%         y = [];
%     end
%     figure, plot(model.w);
%     disp('in constraintCB');
    %
    % load gmodel*
    gm = param.gm;
    gm.fE = gm.fE .* repmat(w(:)',[param.dFeat,1]);
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

function yhat = constraintCB_sw(param, w, g, y)
% function yhat = constraintCB_sw(param, w, g, varargin)
%     if(length(varargin) > 3)
%         y = varargin{1};
%     else
%         y = [];
%     end
%     figure, plot(model.w);
%     disp('in constraintCB');
    %
    % load gmodel*
    gm = param.gm;
    gm.fE(1:13,:) = gm.fE(1:13,:) .* repmat(w(1),[13,param.nP^2]);
    gm.fE(14:end,:) = gm.fE(14:end,:) .* repmat(w(2),[18,param.nP^2]);
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