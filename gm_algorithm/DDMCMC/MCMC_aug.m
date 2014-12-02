%% Solve matching by MCMC sampling
% Data Driven Markov Chain Monte Carlo Sampling for matching
% by Jungmin Lee
function [X,currentW] = MCMC_aug( M, group1, group2, E12,L12, lambda, XGT )


[list(:,1) list(:,2)] = find(E12);
n1 = size(E12,1);
n2 = size(E12,2);

% Assume exact graph matching
% Assume fully connected graph
% Data Driven : Use spectral information

%% parameters
startTemperature = 10; %50;
finishTemperature = 5; %5;
nSamplePerTemp = 5;
TempDecreasingFactor = 0.9;

%% Get the first eigenVector from affinity matrix
options.disp = 0;
% M = sparse(M);
[eigenVectorMat eigenValue] = eigs(M, 1, 'lm', options); % eigensolver for sparse matrix
eigenVector = abs(eigenVectorMat(:,1)); % take first eigenvector

%% initial state (suppose numA = numB)
numMatch = length(unique(list(:,1)));

% e = zeros(n1,1);
% e(randperm(n1,10)) = 1;
% e(1:10) = 1;
% e = sum(reshape(XGT,n1,n2),2);
%     [Xtmp, currentEnergy] = graphmatching(M,group1,group2,E12);
%     tmp = find(Xtmp);
%     rel = sum(M(logical(Xtmp),logical(Xtmp)));
%     [sol,sortid] = sort(rel,'descend');
%     Xin = zeros(size(Xtmp));
%     Xin(tmp(sortid(1:40))) = 1;
%     Xin = reshape(Xin,n1,n2);
%     e = sum(Xin,2);
e = ones(n1,1);
currentE = e;
[Xcurr, currentEnergy] = graphmatching1129(M.*repmat(currentE*currentE',n2,n2),group1,group2,E12);
currentEnergy = lambda(:)' * [currentEnergy; -sum(currentE); -sum(currentE)^2];

temperature = startTemperature;
count = 0;

bestE = currentE;
bestEnergy = currentEnergy;
fprintf('bestEnergy: %f\n',bestEnergy);
nUpdateE = 0;
nBestUpdateE = 0;

if(all(e))
    X = Xcurr;
    return;
end
iter = 0;
% cntAccept = 0;
% cntBestAccept = 0;
%% Sampling Loop
while temperature > finishTemperature 
    iter = iter + 1;
    proposeE = currentE;
    k = 1;    
    sel = randperm(n1,k);
    proposeE(sel) = 1-proposeE(sel);
%     inlist = find(currentE);
%     outlist = find(~currentE);
%     iDel = inlist(randperm(length(inlist),k));
%     iAdd = outlist(randperm(n1-length(inlist),k));
%     proposeE(iDel) = 0;
%     proposeE(iAdd) = 1;

    probGo = 1;
    probBack = 1;

%     [Xprop,proposeEnergy] = graphmatching(M.*repmat(proposeE*proposeE',n2,n2),group1,group2,E12);
%     [Xprop,proposeEnergy] = graphmatching_ipfp(M.*repmat(proposeE*proposeE',n2,n2),group1,group2,E12,L12,Xcurr);
    [Xprop,proposeEnergy] = graphmatching1129_partial(M,group1,group2,L12,proposeE,Xcurr);
%     proposeEnergy = proposeEnergy - lambda*sum(proposeE)^2;
    proposeEnergy = lambda(:)' * [proposeEnergy; -sum(proposeE); -sum(proposeE)^2];
    
    deltaEnergy = proposeEnergy - currentEnergy;
    acceptRatio = exp(deltaEnergy/temperature)*probBack/probGo;
    count = count + 1;
%     if(~isempty(Xprop))
%         nCorrectIter(iter) = Xprop'*XGT;    
%         meanEIter(iter) = (proposeEnergy+lambda*sum(proposeE)^2)/sum(Xprop)^2;
%         eee(iter) = sum(proposeE);
%         egt(iter) = sum(reshape(XGT,n1,n2),2)'*proposeE;
%         figure(1), plot(nCorrectIter);  drawnow;
%         hold on; plot(eee,'r');
%         plot(egt,'g');
%         drawnow;    hold off;
%         figure(2), plot(meanEIter);   drawnow;
%     end

    % Apply MH rule
    if rand < min(acceptRatio, 1)
        currentEnergy = proposeEnergy;
        currentE = proposeE;
        nUpdateE = nUpdateE + 1;
        Xcurr = Xprop;
%         cntAccept = cntAccept + 1;
%         meanE(cntAccept) = (currentEnergy+lambda*sum(currentE)^2) / sum(Xcurr)^2;
%         nCorrect(cntAccept) = XGT'*Xcurr;
%         figure(3), plot(meanE); drawnow;
%         figure(4), plot(nCorrect); drawnow;
        if currentEnergy > bestEnergy
            tmpscore = Xprop'*(M.*repmat(currentE*currentE',n2,n2))*Xprop;
            fprintf('temperature: %f, #nodes = %d, score: %f, scoreReg1: %f\n',temperature,sum(currentE), tmpscore,tmpscore-lambda*sum(currentE)^2);
            bestE = currentE;
            bestEnergy = currentEnergy;
            nBestUpdateE = nBestUpdateE + 1;
%             cntBestAccept = cntBestAccept + 1;
%             meanBestE(cntBestAccept) = meanE(cntAccept);
%             figure(5), plot(meanBestE); drawnow;
        end
    end
           
    if count > nSamplePerTemp*sqrt(numMatch)
        temperature = temperature*TempDecreasingFactor;
        count = 0;
    end
    
%     pause
end
% disp(sprintf('acceptance ratio: %f, best update ratio: %f',nUpdateE/iter, nBestUpdateE/iter));
% [X,sco] = graphmatching(M.*repmat(bestE*bestE',n2,n2),group1,group2,E12);
% sco
if(bestEnergy==0)
    [X, scoIn] = graphmatching1129(M.*repmat(bestE*bestE',n2,n2),group1,group2,E12);
    return;
end
[X,scoIn] = graphmatching1129_partial(M,group1,group2,L12,bestE,[]);
% scoIn = X'*M*X;
% fprintf(sprintf('[DDMCMCE] scoreTot: %f, scoreIn: %f, scoreInOut: %f, scoreOut: %f\n', scoreTot, scoreIn,scoreInOut,scoreOut));
fprintf(sprintf('[DDMCMCE] scoreIn: %f, #iter: %d\n', scoIn,iter));
