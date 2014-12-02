function learnSpectralMatching_Cars(vec, spc_i_indices, spc_j_indices, labels, nodes, nF)

nTrainingPairs = 6;

inputDir = './Cars';

outputDir = './Res_Cars';

mkdir(outputDir);

nEigenIterations = 10;

inputD = dir(inputDir);

train_seq = 2+(1:nTrainingPairs);

%load ../Databases/Caltech101/Data_Training/PairClassifier_final.mat

weights = 0.01*ones(1,5);

thresh.dist = 10000;
thresh.ratio = 10000;
thresh.alpha = 10000;
thresh.theta = 100;

% vec = [];
% spc_i_indices = [];
% spc_j_indices = [];
% labels = [];
% nodes = [];
% gTruth2 = [];

nF = [];

nPair = 0;

bestScore = 0;

maxIterations = 50;

coef = 0.1;

totalTime = 0;

outliersRate=0; nExp = 1;

% for pair = 1:length(train_seq)
% 
%     load(sprintf('%s/%s', inputDir, inputD(train_seq(pair)).name));
%     
%     outliersRate1 = 0;
%     
%     minFeat1 = min(size(features1,1), round(length(gTruth)*(1+outliersRate1)));
%         
%     r1 = [1:nF1, nF1+randperm(min(round(nF1*outliersRate1), size(features1,1)-nF1))];
%        
%        
%     features1 = features1(r1(1:minFeat1),:);
%     
%     minFeat2 = min(size(features2,1), round(length(gTruth)*(1+outliersRate)));
%         
%     r2 = [1:nF2, nF2+randperm(min(round(nF2*outliersRate), size(features2,1)-nF2))];
%        
%     features2 = features2(r2(1:minFeat2),:);
%   
%     edges1 = zeros(minFeat1);
%     
%     tri = delaunay(features1(:,1), features1(:,2));
%                
%     for q = 1:size(tri,1)
%        
%         edges1(tri(q,:), tri(q,:)) = 1;
%         
%     end
%     
%     
%     edges2 = zeros(minFeat2);
%     
%     tri = delaunay(features2(:,1), features2(:,2));
%                
%     for q = 1:size(tri,1)
%        
%         edges2(tri(q,:), tri(q,:)) = 1;
%         
%     end
%      
% %    [u ind] = unique(gTruth);
%     
% %     f = find(u);
% %     
% %     u = u(f);
% %     
% %     ind = ind(f);
% %     
% %     features1 = features1(ind,:);
% %     
% %     features2 = features2(u,:);
%     
% %     r = randperm(size(features2,1));
% %     
% %     gTruth(r) = 1:length(r);
% %     
% %     features2 = features2(r,:);
%         
%     [vec_aux, spc_i_indices_aux, spc_j_indices_aux, labels_aux, nodes_aux] = getPairWiseInfo_Faces_2(features1, features2, edges1, edges2, thresh);
% 
%     nPair = nPair+1;
% 
%     vec{nPair} = vec_aux;
%     spc_i_indices{nPair} = spc_i_indices_aux;
%     spc_j_indices{nPair} = spc_j_indices_aux;
%     labels{nPair} = labels_aux;
%     nodes{nPair}  = nodes_aux;
% 
%     nF(nPair) = nF1;
%     
%     
%     %gTruth2{pair} = gTruth;
% end

w = weights;

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

param = [];

for iter = 1:maxIterations


    dw = zeros(1, length(w));

    for pair = 1:nPair

        iter
        
        pair

        gt = zeros(length(labels{pair}), 1);
        gt(find(nodes{pair} == labels{pair})) = 1;;
        
        gt(find(nodes{pair} > nF(pair) | labels{pair} > nF(pair))) = 0;
        
        
        gt = gt/norm(gt);
         
        nCorrect = length(find(gt));
        
        [M,spc_affinities] = getMatchingMatrix(vec{pair}, spc_i_indices{pair}, spc_j_indices{pair}, labels{pair}, nodes{pair}, 0.5, w);
 
        %keyboard;
       
        %tic;
        
        for var = 1:length(w)
            
            dM = spconvert([spc_i_indices{pair} spc_j_indices{pair} spc_affinities.*(-vec{pair}(:,var))]);
            dM = dM + dM';
            
            [v, dv] = getEigenDerivative(M, dM, nEigenIterations);
            
            aux_vec = getSol2(v, labels{pair}, nodes{pair});
            
            dw(var) = dw(var) + (v - aux_vec)'*dv;

             
        end
        
        selfCorrelation(iter) = selfCorrelation(iter) + v'*aux_vec;
        
        aux_vec(find(nodes{pair} > nF(pair) | labels{pair} > nF(pair)))=0;
        aux_vec = aux_vec/(norm(aux_vec)+eps);
             
        v(find(nodes{pair} > nF(pair) | labels{pair} > nF(pair)))=0;
        v = v/(norm(v)+eps);
        %monitor eigenvalues, error matrices ....
 
        [eve, eva] = eigs(M,2);
        
        eigenvec{iter} = v;
        
        gtruth{iter} = gt;
        
        correlationScore(iter) = correlationScore(iter) + v'*gt;
       
        softError(iter) = softError(iter) + norm(v - gt);
         
        %[sol, score_aux, V_aux] = eigen_solution(M, labels{pair}, nodes{pair}, 0);
        
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
        
        vRatio_ideal(iter) = vRatio_ideal(iter) + (-kr + sqrt(kr^2 + 4*(outlierRate-1)*(p0/(p1 + eps))^2))/(2*(outlierRate-1)*(p0/(p1 + eps)));
       
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
        
    w = w - coef*dw;

    param{iter} = w;
    
    save(sprintf('%s/learning_unsup_%f_%d.mat', outputDir, outliersRate, nExp), ...
        'softError', 'hardScore', 'param', 'pRatio', 'vRatio', 'vRatio_ideal', 'w', 'totalTime', 'correlationScore', 'gtruth', 'eigenvec', ...
         'selfCorrelation', 'eigengap', 'lambda1', 'lambda2', 'train_seq');
end

