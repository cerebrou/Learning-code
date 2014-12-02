function ww = marius_sw(Ms, param)

nEigenIterations = 10;
nF = param.nP;
n1 = nF;
n2 = nF;

thresh.dist = 10000;
thresh.ratio = 10000;
thresh.alpha = 10000;
thresh.theta = 100;

weights = 0.01*ones(1,2);

% nPair = length(vec);
nPair = length(Ms);
bestScore = 0;
maxIterations = 200;
coef = 0.0001; % 0.00001
totalTime = 0;
outliersRate=0; nExp = 1;

%%
w = weights;
featrange = [1,13; 14,31];

softError = zeros(1, maxIterations);
hardScore = zeros(1, maxIterations);
correlationScore = zeros(1, maxIterations);
selfCorrelation = zeros(1, maxIterations);
lambda1 = zeros(1, maxIterations);
lambda2 = zeros(1, maxIterations);
eigengap = zeros(1,maxIterations);
errorBound = zeros(1, maxIterations);
pRatio = zeros(1, maxIterations);
vRatio = zeros(1, maxIterations);
vRatio_ideal = zeros(1, maxIterations);

var2ind = find(~tril(ones(n1)));

for pair = 1:nPair
    nodes{pair} = repmat(1:n1,1,n2)';
    labels{pair} = kron(1:n2,ones(1,n1))';
    Mvec{pair,1} = sum(Ms{pair}(:,:,featrange(1,1):featrange(1,2)),3);
    Mvec{pair,2} = sum(Ms{pair}(:,:,featrange(2,1):featrange(2,2)),3);
end

for iter = 1:maxIterations

    dw = zeros(1, numel(w));

    for pair = 1:nPair

        gt = zeros(length(labels{pair}), 1);
        gt(find(nodes{pair} == param.labels{pair}(labels{pair}))) = 1;
        gt(find(nodes{pair} > nF | labels{pair} > nF)) = 0;
%         gt = gt/norm(gt);         
        nCorrect = length(find(gt));
        
%         [M,spc_affinities] = getMatchingMatrix(vec{pair}, spc_i_indices{pair}, spc_j_indices{pair}, labels{pair}, nodes{pair}, 0.5, w);
        M = Mvec{pair,1}*w(1) + Mvec{pair,2}*w(2);
        
        for var = 1:numel(w)
%             dM = spconvert([spc_i_indices{pair} spc_j_indices{pair} spc_affinities.*(-vec{pair}(:,var))]);
%             dM = M.*sum(Mvec(:,:,featrange(var,1):featrange(var,2)),3);
            dM = Mvec{pair,var};
            
            [v, dv] = getEigenDerivative(M, dM, nEigenIterations);
            aux_vec = getSol2(v, labels{pair}, nodes{pair});
            dw(var) = dw(var) + (v - aux_vec)'*dv;
        end
        
        selfCorrelation(iter) = selfCorrelation(iter) + v'*aux_vec;
        
        aux_vec(find(nodes{pair} > nF | labels{pair} > nF))=0;
        aux_vec = aux_vec/(norm(aux_vec)+eps);
             
        v(find(nodes{pair} > nF | labels{pair} > nF))=0;
        v = v/(norm(v)+eps);
        %monitor eigenvalues, error matrices ....
 
        [eve, eva] = eigs(M,2);
        eigenvec{iter} = v;
        gtruth{iter} = gt;
        correlationScore(iter) = correlationScore(iter) + v'*gt;
        softError(iter) = softError(iter) + norm(v - gt);
         
        f = find(gt);
        f0 = find(gt == 0);
        fMatch  = find(aux_vec);
        fMatch0 = find(aux_vec==0);
        f2 = find(nodes{pair}(fMatch) == labels{pair}(fMatch));
        score_curr = length(f2)/nCorrect;
        hardScore(iter) = hardScore(iter) + score_curr;
        
        %monitor eigenvalues, error matrices ....
        eigengap(iter) = eigengap(iter) + (eva(1,1) - eva(2,2))/mean2(M);
        p1 = mean2(M(f,f));
        p0 = (mean2(M(f0,f0))*length(f0)^2 + 2*mean2(M(f0,f))*length(f0)*length(f))/(length(f0)^2 + 2*length(f0)*length(f));
        pRatio(iter) = pRatio(iter) + p0/(p1 + eps);
        v1 = mean(v(f));
        v0 = mean(v(f0));
        vRatio(iter) = vRatio(iter) + v0/(v1 + eps);
        outlierRate = length(gt)/nCorrect;
        kr = 1 - (outlierRate - 1)*p0/(p1 + eps);
        vRatio_ideal(iter) = vRatio_ideal(iter) + ...
                    (-kr + sqrt(kr^2 + 4*(outlierRate-1)*(p0/(p1 + eps))^2)) ...
                    /(2*(outlierRate-1)*(p0/(p1 + eps)));      
        lambda1(iter) = lambda1(iter) + eva(1,1);
        lambda2(iter) = lambda2(iter) + eva(2,2);
        
    end
   
    eigengap(iter) = eigengap(iter)/nPair;
    pRatio(iter) = pRatio(iter)/nPair;
    vRatio(iter) = vRatio(iter)/nPair;
    vRatio_ideal(iter) = vRatio_ideal(iter)/nPair;
    lambda1(iter) = lambda1(iter)/nPair;
    lambda2(iter) = lambda2(iter)/nPair;
    hardScore(iter) = hardScore(iter)/nPair;
    correlationScore(iter) = correlationScore(iter)/nPair;
    selfCorrelation(iter)  = selfCorrelation(iter)/nPair;
    
%     dw = dw + 2*w;
    w = w - coef/(iter+1)*dw;
%     w = w / (norm(w)+eps);
    warch{iter} = w;
    bv(iter) = aux_vec' * v;
end

figure, plot([warch{:}]);
figure, plot(eigengap); title('eigengap');
figure, plot(pRatio); title('pRatio');
figure, plot(vRatio); title('vRatio');
figure, plot(bv);

ww = w;
end