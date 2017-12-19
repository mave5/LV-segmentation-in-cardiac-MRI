%-- re-map small mask to big mask

function out_mask = remap_mask(in_mask,m_cnt,I)

    [y_max, x_max]=size(I);
    
    M=size(in_mask,1);
       
    % center
    m_cnt_x=m_cnt(1);    
    m_cnt_y=m_cnt(2);

    
    % top left corner
    x1=m_cnt_x-M/2;
    y1=m_cnt_y-M/2;

    % bottom right corner
    x4=m_cnt_x+M/2-1;
    y4=m_cnt_y+M/2-1;

    out_mask=zeros(x_max,y_max);
    out_mask(y1:y4,x1:x4)=in_mask;

    
end  