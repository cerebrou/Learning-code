ww = sum(ww.^2);
maxthres = max(ww(:));
a = 0 : maxthres/100 : maxthres;
n = length(a);
num = zeros(n,1);
for i = 1 : n
    num(i) = sum(ww>a(i));
end
figure, plot(num, 'linewidth', 2, 'color','r')
hold on; 
legend('99.4','92.3');
title('(l2 + l21)-norm regularization');