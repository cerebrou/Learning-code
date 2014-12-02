function showPointsBoth(source, points)
    figure();
    subplot(1,2,1); scatter(source.xy(:,1), source.xy(:,2));
    subplot(1,2,2); scatter(points.xy(:,1), points.xy(:,2));
    
end