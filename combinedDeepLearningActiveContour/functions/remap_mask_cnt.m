%-- roi mask will be converted to the original image size

function out_mask = remap_mask_cnt(in_mask,I,m_cnt)

    [y_max, x_max]=size(I);
    
    M=size(in_mask,1);
    M2=floor(M/2);       
    
    % center
    m_cnt_x=m_cnt(1);    
    m_cnt_y=m_cnt(2);

    
    % top left corner
    x1=max(m_cnt_x-M2,1);
    y1=max(m_cnt_y-M2,1);

    % bottom right corner
    if x1==1
        x4=M;
    else
        x4=min(m_cnt_x+M2,x_max);
        if x4==x_max
            x1=x_max-M+1;
        end
    end
    if y1==1
        y4=M;
    else
        y4=min(m_cnt_y+M2,y_max);
        if y4==y_max
            y1=y_max-M+1;
        end
    end
    out_mask=zeros(y_max,x_max);
    out_mask(y1:y4,x1:x4)=in_mask;
    
end  
