function [M,  spc_affinities]  = getMatchingMatrix(vec, spc_i_indices, spc_j_indices, labels, nodes, minContribution, weights)

%tic;

%MIN_G = 2;

%nEntries = length(spc_i_indices);

%spc_affinities = zeros(nEntries,1);

%w = weights; %[weights(1), weights(2), weights(3), weights(4)];
 
%weights = weights/2;
 
%weights = [6.1334 -2.3655 -2.2047 -3.2587 -0.2215 -0.1957 -0.0565 -0.8437];

%spc_affinities = 1./(1 + exp(-vec*weights'));


spc_affinities = exp(-vec*weights');


% for i = 1:nEntries
%     
%     if spc_i_indices(i) == spc_j_indices(i)
%         continue;
%     end
%       
%     gScore = w*vec(i,:)';
%     
% %     if gScore > MIN_G
% %         continue;
% %     end
%     
%     spc_affinities(i) =  exp(-gScore);
%     
%     %keyboard;
%     
% end


%fprintf('Calling spconvert to create the sparse matrix M to its transpose\n');


%keyboard;

M = spconvert([spc_i_indices spc_j_indices spc_affinities]);

%size(M)

%fprintf('Adding M to its transpose\n');
M = M + M';

%fprintf(' matrix fullness: %f \n', (200*nz_k)/nAssig^2);


%% finding the assignments 

% [sol, score, V] = eigen_solution(M, labels, nodes, minContribution);
% f = find(sol);
% 
% labels = labels(f);
% nodes = nodes(f);
% 
% %toc;

return