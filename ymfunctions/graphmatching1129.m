function [ysol,sco] = graphmatching1129(M,group1,group2,E12)

nP = size(group1,2);
nQ = size(group2,2);
yraw = RRWM(double(M), group1, group2);
ysol = greedyMapping(yraw, group1, group2);
ysol = reshape(ysol, [nP,nQ]);
% [y, ~] = find(ysol');
sco = ysol(:)'*M*ysol(:);

end