function points = makeSourcePoints(nP,varmax)

points.xy = rand(nP,2);
points.vars = varmax*rand(nP,1);

end