function w = trainWeight_marius(M,y,nP1,dFeat)
% INPUT
% M: (nTrain by 1) cell
% M{iTrain}: (nMatches by nMatches by dFeat) tensor
% y: (nTrain by 1) cell
% y{iTrain}: (nMatches by 1) binary vector
% 
% OUTPUT

nMatches = size(M{1},1);
nTrain = length(M);
w = ones(2,1);
itermax = 300;
eta = 0.1;
N = 10;
nGroup = 2;
Y = cell(nTrain,1);
for iTrain = 1 : nTrain
    Y{iTrain} = zeros(nMatches,1);
    Y{iTrain}((1:nP1)' + (y{iTrain}(:)-1)*nP1) = 1;
    M{iTrain} = M{iTrain} /dFeat;
end
groups{1} = 1:13;
groups{2} = 14:31;
wgroup = zeros(dFeat,1);
wgroup(groups{1}) = w(1);
wgroup(groups{2}) = w(2);

for iter = 1 : itermax
    grad = zeros(size(w));
    
    for iGroup = 1 : nGroup
        for iTrain = 1 : nTrain
            Mtot = zeros(nMatches);
            for iFeat = 1 : dFeat
                Mtot = Mtot + M{iTrain}(:,:,iFeat)*wgroup(iFeat);
            end
            Mn1 = Mtot^N * ones(nMatches,1);
            Mn1norm = norm(Mn1,2);
            Mn1pPrev = sum(M{iTrain}(:,:,groups{iGroup}),3)*ones(nMatches,1);
            for n = 2 : N
                Mn1p = M{iTrain}(:,:,iFeat)*Mtot^(n-1)*ones(nMatches,1) ...
                        + Mtot*Mn1pPrev;
                Mn1pPrev = Mn1p;
            end
            vp = (Mn1p*Mn1norm - Mn1/Mn1norm*Mn1'*Mn1p) / Mn1norm^2;
            grad(iGroup) = grad(iGroup) + Y{iTrain}'*vp;
        end
    end
    w = w + eta*grad;
    figure(1),plot(w);
    iter
end

end

% 
% function w = trainWeight_marius(M,y)
% % INPUT
% % M: (nTrain by 1) cell
% % M{iTrain}: (nMatches by nMatches by dFeat) tensor
% % y: (nTrain by 1) cell
% % y{iTrain}: (nMatches by 1) binary vector
% %
% % OUTPUT
% 
% nMatches = size(M,1);
% dFeat = size(M,3);
% w = zeros(dFeat,1);
% itermax = 300;
% eta = 0.1;
% 
% for iter = 1 : itermax
%     grad = zero(size(w));
%     
%     for iFeat = 1 : dFeat
%         for iTrain = 1 : nTrain
%             Mtot = zeros(nMatches,nMatches);
%             for i = 1 : dFeat
%                 Mtot = Mtot + w(i)*M{iTrain}(:,:,i);
%             end
%             Mtot = exp(Mtot);
%             Mn1 = Mtot^N * ones(nMatches,1);
%             Mn1norm = norm(Mn1,2);
%             Mn1pPrev = sum(Mtot.*M{iTrain}(:,:,iFeat),2);
%             for n = 2 : N
%                 Mn1p = M{iTrain}(:,:,iFeat)*M^(n-1)*ones(nMatches,1) ...
%                         + M{iTrain}*Mn1pPrev;
%                 Mn1pPrev = Mn1p;
%             end
%             vp = (Mn1p*Mn1norm - Mn1/Mn1norm*Mn1'*Mn1p) / Mn1norm^2;
%             grad = grad + y{iTrain}(:)'*vp;
%         end
%         w(iFeat) = w(iFeat) + eta*grad;
%     end
% end
% 
% end