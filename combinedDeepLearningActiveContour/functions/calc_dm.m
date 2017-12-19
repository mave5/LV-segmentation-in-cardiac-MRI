% calculate dice metric from contours
function dm = calc_dm(autoPoints,manualPoints,para)
    %calc dice metric
    auto_mask = poly2mask (autoPoints(:,1),autoPoints(:,2),double(para.width),double(para.height));
    manual_mask = poly2mask (manualPoints(:,1),manualPoints(:,2),double(para.width),double(para.height));

    auto_size = sum(auto_mask(:)>0);
    manual_size = sum(manual_mask(:)>0);
    intersect_size = sum((auto_mask(:) + manual_mask(:))==2);
    dm = 2 * intersect_size / (auto_size + manual_size);
end