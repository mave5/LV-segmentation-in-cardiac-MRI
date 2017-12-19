% compute metrics

function [dm1,dm2,PD1,PD2,HD1,HD2]=eval_metrics(I,m_cnt,auto_seg,manualPoints,para)

    % resize mask in the original image size   
    auto_seg_r=remap_mask(auto_seg,m_cnt,I);
    
    % convex hull of resized mask
    c_auto_seg_r=bwconvhull(auto_seg_r);
    
    % convert mask to contour
    autoPoints1=contourc(double(auto_seg_r), [0 0]); autoPoints1=autoPoints1(:,2:end)';
    autoPoints2=contourc(double(c_auto_seg_r), [0 0]);autoPoints2=autoPoints2(:,2:end)';
    
if ~isempty(autoPoints1) && ~isempty(autoPoints2)
% Dice Metric
dm1 = calc_dm(autoPoints1,manualPoints,para);
dm2 = calc_dm(autoPoints2,manualPoints,para);

% Perpendicular Distance
PD1 = calc_dist(autoPoints1,manualPoints,para);
PD2 = calc_dist(autoPoints2,manualPoints,para);

% Hausdorff distance
HD1 = hausdorff( manualPoints, autoPoints1); 
HD2 = hausdorff( manualPoints, autoPoints2); 
else
dm1=0;
dm2=0;
PD1=100;
PD2=100;
HD1=100;
HD2=100;
end

end