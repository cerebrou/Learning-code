function s = kfold(L,k,i)
% INPUT
% L: (N by 1) index list
% k: k-partition
% i: i-th group
%
% OUTPUT
% i-th set of k-fold when datasize is N

    assert(i>=1 && i<=k);
    N = length(L);
    step = floor(N/k*i);
    starti = 1 + (i-1)*step;
    endi = starti + step -1;
    if(k==i),   endi = N;   end
    s = L(starti:endi);

end