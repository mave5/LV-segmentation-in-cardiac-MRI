function t_yLV_clean=clean_segs(t_yLVh)
%%
for k=1:size(t_yLVh,3)
    
    BW=t_yLVh(:,:,k);
    BW2= imfill(BW,'holes');
    cc = bwconncomp(BW2,4);
    stats = regionprops(cc, 'Area');
    [area_max,idx] = max([stats.Area]);
    temp= ismember(labelmatrix(cc), idx);
    t_yLV_clean(:,:,k) = imfill(temp,'holes');
    t_yLV_clean(:,:,k)=bwconvhull(t_yLV_clean);
end


