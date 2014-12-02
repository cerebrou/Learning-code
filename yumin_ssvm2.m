function model = yumin_ssvm2(input)
% Subgradient descent learning for graph matching
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
parm.nP = input.nP;
parm.gm = input.gm;
parm.dFeat = input.dFeat;

% dFeat = size(x{1}.xy,1);
% nFeat = size(x{1}.xy,2);
nTrain = length(x);
iterMax = 300;
C = 1; %(HARG-point)C=1, learningRate=0.1, iterMax=300
%(DW-point-l1)C=3, learningRate=0.1, iterMax=300
%(DW-CMU-l2)C=1, learningRate=0.1, iterMax=300
%(DW-CMU-l1)C=3, learningRate=0.01, iterMax=300
%(HARG-CMU-l1)C=3 or 1, learningRate=0.01, iterMax=300
%(HARG-CMU-wl21)C=3, learningRate=0.01, iterMax=300
lambda = 0.3;   %(1,0.3,0.003)  %(0.1,0.03)
learningRate = 0.003;

w = 0.01*(rand(input.dimension,1)-0.5);
yhat = cell(nTrain,1);
converge = false;
grad = zeros(nTrain, length(w));
objective = zeros(iterMax,1);
for iter = 1 : iterMax    
    if(mod(iter,10)==1),    disp([num2str(iter) '-th iteration']);  end
    if(converge),   break;  end
    
%     w = w - learningRate/iter * (2*w+C/nTrain*sum(grad,1)');
%     w = w - learningRate/iter * (1/2*sign(w)+C/nTrain*sum(grad,1)');
    %
    ww = reshape(w,input.dFeat,input.nP^2);
    gradReg = zeros(size(ww));
%     sparseInd = sqrt(sum(ww.^2,1)) > 0.0000001;
%     gradReg(:,sparseInd) = bsxfun(@rdivide, ww(:,sparseInd), sqrt(sum(ww(:,sparseInd).^2,1)));
    gradReg = bsxfun(@rdivide, ww, sqrt(sum(ww.^2,1)));
% %     w = w - learningRate/iter * (lambda*gradReg(:)+C/nTrain*sum(grad,1)');
    w = w - learningRate/sqrt(iter) * (2*w + lambda*gradReg(:) + C/nTrain*sum(grad,1)');
    %
    totalDelta = 0;
    switch(input.learningType)
        case 'HARG'
            parfor iTrain = 1 : nTrain
                yhat{iTrain} = constraintCB(parm,w,x{iTrain},y{iTrain});
                grad(iTrain,:) = featureCB(parm,x{iTrain},yhat{iTrain}) - featureCB(parm,x{iTrain},y{iTrain});
                totalDelta = totalDelta + lossCB([],y{iTrain},yhat{iTrain});
            end
        case 'dw'
            parfor iTrain = 1 : nTrain
                yhat{iTrain} = constraintCB_dw(parm,w,x{iTrain},y{iTrain});
                grad(iTrain,:) = featureCB_dw(parm,x{iTrain},yhat{iTrain}) - featureCB_dw(parm,x{iTrain},y{iTrain});
                totalDelta = totalDelta + lossCB([],y{iTrain},yhat{iTrain});
            end
        case 'sw'
            parfor iTrain = 1 : nTrain
                yhat{iTrain} = constraintCB_sw(parm,w,x{iTrain},y{iTrain});
                grad(iTrain,:) = featureCB_sw(parm,x{iTrain},yhat{iTrain}) - featureCB_sw(parm,x{iTrain},y{iTrain});
                totalDelta = totalDelta + lossCB([],y{iTrain},yhat{iTrain});
            end
    end
%     objective(iter) = w(:)'*w(:) + C/nTrain*totalDelta;
%     objective(iter) = sqrt(w(:)'*w(:)) + C/nTrain*totalDelta;
%     objective(iter) = lambda*sum(sqrt(sum(w.^2,1))) + C/nTrain*totalDelta;
    objective(iter) = w(:)'*w(:) + sum(lambda*sqrt(sum(w.^2,1))) + C/nTrain*totalDelta;
    % convergence check
%     if(iter>=2 && dist(objective(iter-1), objective(iter))<input.e)
%         converge = true;
%     end
figure(1), plot(objective(max(1,iter-50):iter));
% figure(2), plot(w);
end

figure; plot(objective);

% model.w = reshape(w,size(x{1}));
model.w = w;

end


function delta = lossCB(param, y, ybar)
    delta = double(sum(y~=ybar)) / length(y);
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