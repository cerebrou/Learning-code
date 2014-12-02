function model = yumin_ssvm(input)
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

% dFeat = size(x{1}.xy,1);
% nFeat = size(x{1}.xy,2);
nTrain = length(x);
iterMax = 300;
C = 1; %C=1, learningRate=0.1, iterMax=300
lambda = 0.1;   %(1,0.3)  %(0.1,0.03)
learningRate = 0.1;

w = 0.01*(rand(input.dimension,1)-0.5);
yhat = cell(nTrain,1);
converge = false;
grad = zeros(nTrain, length(w));
objective = zeros(iterMax,1);
for iter = 1 : iterMax    
    if(mod(iter,10)==1),    disp([num2str(iter) '-th iteration']);  end
    if(converge),   break;  end
    
%     w = w - learningRate/iter * (2*w+C/nTrain*sum(grad,1)');
    w = w - learningRate/iter * (1/2*sign(w)+C/nTrain*sum(grad,1)');
    %
%     ww = reshape(w,input.dFeat,input.nP^2);
%     gradReg = zeros(size(ww));
%     sparseInd = sqrt(sum(ww.^2,1)) > 0.0000001;
%     gradReg(:,sparseInd) = bsxfun(@rdivide, ww(:,sparseInd), sqrt(sum(ww(:,sparseInd).^2,1)));
% %     w = w - learningRate/iter * (lambda*gradReg(:)+C/nTrain*sum(grad,1)');
%     w = w - learningRate/sqrt(iter) * (2*w + lambda*gradReg(:) + C/nTrain*sum(grad,1)');
    %
    totalDelta = 0;
    parfor iTrain = 1 : nTrain
        yhat{iTrain} = constraintCB(parm,w,x{iTrain},y{iTrain});
        grad(iTrain,:) = featureCB(parm,x{iTrain},yhat{iTrain}) - featureCB(parm,x{iTrain},y{iTrain});
        totalDelta = totalDelta + lossCB([],y{iTrain},yhat{iTrain});
    end
%     objective(iter) = w(:)'*w(:) + C/nTrain*totalDelta;
    objective(iter) = sqrt(w(:)'*w(:)) + C/nTrain*totalDelta;
%     objective(iter) = sum(lambda*sqrt(sum(w.^2,1))) + C/nTrain*totalDelta;
%     objective(iter) = w(:)'*w(:) + sum(lambda*sqrt(sum(w.^2,1))) + C/nTrain*totalDelta;
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

function yhat = constraintCB(param, w, g, y)
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

% function w = featureCB(param, x, y)
%     w = sparse(y*x(:)/2);
% end
% 
% function yhat = constraintCB(param, w, x, y)
%     if dot(y*x(:), w) > 1
%         yhat = y;
%     else
%         yhat = -y;
%     end
% end
