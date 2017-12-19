function t_yLV_clean=clean_segs(t_yLVh)
%%
for k=1:size(t_yLVh,3)
    BW=t_yLVh(:,:,k);
    cc = bwconncomp(BW,4);
    stats = regionprops(cc, 'Area');
    [area_max,idx] = max([stats.Area]);
    ylvc_ch = ismember(labelmatrix(cc), idx);
    t_yLV_clean(:,:,k)=bwconvhull(ylvc_ch);
end


