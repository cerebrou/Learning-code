function [ X score ] = RRWM_v2( problem, config )
% Reweighted Random Walk Matching 
% by Minsu Cho, Jungmin Lee, and Kyoung Mu Lee
%
% input:   affinityMatrix  - affinity matrix (the larger, the more similar)
%          conflictMatrix1 - confliction matrix of feature set 1
%          conflictMatrix2 - confliction matrix of feature set 2
%
% output:  X               - nMatch by 1 indicator vector of the solution  
%          score           - nMatch by 1 scalar vector

% settings
isDisplay = 0;             % flag for visualization
c = 0.2;                   % prob. for walk or reweighted jump?
amp_max = 30;              % maximum value for amplification procedure
iterMax = 300;             % maximum iterations for random walks 
thresConvergence1=1e-30;    % convergence threshold for random walks
tolC=1e-3;                 % convergence threshold for the Sinkhorn method
convergengeType = 1;

% load data
%M = exp(-cdata.affinityMatrix/100);
M = problem.affinityMatrix;
nMatch = length(M);

% get groups for bistochastic normalization
[idx1 ID1] = make_groups(problem.conflictMatrix1);
[idx2 ID2] = make_groups(problem.conflictMatrix2);

% eliminate conflicting elements to prevent conflicting walks
conf = problem.conflictMatrix1 | problem.conflictMatrix2;
M = M.*(~conf);

% note that this matrix is column-wise stochastic
d = sum(M, 1); % degree : col sum
maxD = max(d);
M = M ./ maxD; % nomalize by the max degree

% initialize answer
prev_score = ones(nMatch,1)/nMatch; % buffer for the previous score
prev_score2 = prev_score;         % buffer for the two iteration ahead
prev_assign = ones(nMatch,1)/nMatch; % buffer for Sinkhorn result of prev_score

% for convergence check of power iteration
thresConvergence2 = length(prev_score)*norm(M,1)*eps;
la = prev_score'*M*prev_score;

bCont = 1;  iter_i = 0;

%show_script_for_RWPRM;
if isDisplay,    figure(1); end

%% start main iteration
while bCont && iter_i < iterMax
    
    iter_i = iter_i + 1;
    
    %% random walking with reweighted jumps
    cur_score = M * ( c*prev_score + (1-c)*prev_assign );
    
    sumCurScore = sum(cur_score); % normalization of sum 1
    if sumCurScore>0, cur_score = cur_score./sumCurScore; end
    
    %% update reweighted jumps
    cur_assign = cur_score;
    % attenuate small values and amplify large values
    amp_value = amp_max/ max(cur_assign);  % compute amplification factor
    cur_assign = exp( amp_value*cur_assign );  % amplification 
    
    % Sinkhorn method of iterative bistocastic normalizations
    cur_assign = mexBistocNormalize_match(cur_assign, int32(idx1), int32(ID1), int32(idx2), int32(ID2), tolC);
      
    sumCurAssign = sum(cur_assign); % normalization of sum 1
    if sumCurAssign>0, cur_assign = cur_assign./sumCurAssign; end
    
    %% Check the convergence of random walks
    if convergengeType == 1
        diff1 = sum((cur_score-prev_score).^2);
        diff2 = sum((cur_score-prev_score2).^2); % to prevent oscillations
        diff_min = min(diff1, diff2);
        if diff_min < thresConvergence1
            bCont = 0;
        end
    else
        normed_cur_score = cur_score/norm(cur_score);
        if norm(M*normed_cur_score - la*normed_cur_score,1) < thresConvergence2
            bCont = 0;
        end
    end
    
    prev_score2 = prev_score;
    prev_score = cur_score;
    prev_assign = cur_assign;
end % end of main iteration

% fprintf(' %d',iter_i);
score = cur_score;

%% discretize the solution    
%X = greedyMapping(cur_assign, conf);
X = greedyMapping(cur_score, conf);
% X = feval(AlgConfig.postProcess, problem, X2);

% %     X2=computeXorthonormal(X2);
% %X = computeDiscreteCorrespondancesGreedy(X2,E12);
% X = discretisationMatching_hungarian(X2,E12);
% X = reshape(X(:), n1, n2);


end

function display_result(X,idxFig)
if nargin>1
    fig = figure(idxFig);
    set(fig, 'MenuBar', 'none');
    set(fig, 'Toolbar', 'none');
end

%imagesc(X);
bar3(X,'detachted');
drawnow;
%caxis([0,1]);
%pause(0.01);
%disp2('bounds(X(:))');
%colormap(gray);

end