function model = yumin_ssvm_cutting_test(input)
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
parm.dFeat = input.dFeat;
nTrain = length(x);

itermax = 100;
thres = 0.5;

workingSet = [];
for iter = 1 : itermax
    workingSetChanged = false;
    z = solvecvx(parm,workingSet,x,y);
    for iTrain = 1 : nTrain
        yhat = constraintCB(parm,z(1:parm.dFeat),x{iTrain},y{iTrain});
        if(lossCB(parm,yhat,y{iTrain})>thres)
            if(isempty(workingSet) || ~ismember(iTrain,[workingSet.iTrain]) ...
                    || lossCB(parm,workingSet([workingSet.iTrain]==iTrain).yhat,yhat) ~= 0)
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
    H = zeros(dFeat+nTrain);
    H(1:dFeat,1:dFeat) = eye(dFeat);
    C = 10;
%     C = 1/4;
    
    cvx_begin
        variables z(dFeat+nTrain,1)
%         minimize( z'*H*z + C*[zeros(1,dFeat), ones(1,nTrain)]*z )
        minimize( norm(z(1:dFeat),1) + C/nTrain*ones(1,nTrain)*z((dFeat+1):end) )
        subject to
            z((dFeat+1):end) >= 0;
            for iActive = 1 : nActive
                iTrain = workingSet(iActive).iTrain;
                yhat = workingSet(iActive).yhat;
                z(1:dFeat)'*(featureCB(parm, x{iTrain},y{iTrain})-featureCB(parm, x{iTrain},yhat))  >= lossCB(y{iTrain},yhat) - z(dFeat+iTrain);
            end
    cvx_end
    
    % n-slack
    % cvx/TFOCS/quadprog
end

function delta = lossCB(parm, y, ybar)
    delta = double(y ~= ybar) ;
end

function w = featureCB(parm, x, y)
    w = sparse(y*x) ;
    w = w(:);
end

function yhat = constraintCB(parm, w, x, y)
%     dFeat = parm.dFeat;
% 	if dot(y*x, w) > 1, yhat = y ; else yhat = - y ; end
	if dot(y*x, w) > dot(-y*x, w), yhat = y ; else yhat = - y ; end
end