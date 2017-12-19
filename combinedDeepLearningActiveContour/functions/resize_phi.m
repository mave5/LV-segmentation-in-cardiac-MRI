%-- converts a mask to a Signed Distance Function (SDF)

function [phi_r,mask_r] = resize_phi(phi,m_cnt,I)

    [x_max, y_max]=size(I);
    
    mask=phi<=0;
    M=size(mask,1);
       
    % center
    m_cnt_x=m_cnt(1);    
    m_cnt_y=m_cnt(2);

    
    % top left corner
    x1=m_cnt_x-M/2;
    y1=m_cnt_y-M/2;

    % bottom right corner
    x4=m_cnt_x+M/2-1;
    y4=m_cnt_y+M/2-1;

    mask_r=zeros(x_max,y_max);
    %mask_r(x1:x4,y1:y4)=mask;
    mask_r(y1:y4,x1:x4)=mask;
    
    % convert mask to phi
    phi_r=mask2phi(mask_r);
    
end  