clear; clc;

% ------------------------------------------------------------------
%                                                      Generate data
% ------------------------------------------------------------------

th = pi/3 ;
c = cos(th) ;
s = sin(th) ;

patterns = {} ;
labels = {} ;
for i=1:100
    patterns{i} = diag([2 .5]) * randn(2, 1) ;
    labels{i}   = 2*(randn > 0) - 1 ;
    patterns{i}(2) = patterns{i}(2) + labels{i} ;
    patterns{i} = [c -s ; s c] * patterns{i}  ;
end

% ------------------------------------------------------------------
%                                                    Run SVM struct
% ------------------------------------------------------------------
  
model = yumin_ssvm(patterns, labels);

% ------------------------------------------------------------------
%                                                              Plots
% ------------------------------------------------------------------

w = model.w;
figure(1) ; clf ; hold on ;
x = [patterns{:}] ;
y = [labels{:}] ;
plot(x(1, y>0), x(2,y>0), 'g.') ;
plot(x(1, y<0), x(2,y<0), 'r.') ;
set(line([0 w(1)], [0 w(2)]), 'color', 'y', 'linewidth', 4) ;
xlim([-3 3]) ;
ylim([-3 3]) ;
set(line(10*[w(2) -w(2)], 10*[-w(1) w(1)]), ...
  'color', 'y', 'linewidth', 2, 'linestyle', '-') ;
axis equal ;
set(gca, 'color', 'b') ;
w