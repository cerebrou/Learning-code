function G = edgeidx2G(edgeidx1, n)

m = size(edgeidx1,1);
G = zeros(n,m);
i = edgeidx1 + n*(repmat((1:m)',1,2)-1);
G(i) = 1;

end