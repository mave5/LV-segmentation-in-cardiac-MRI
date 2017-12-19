%-- converts a mask to a Signed Distance Function (SDF)

function phi = mask2phi(mask,m_cnt,I)

if nargin==3
    [x_max, y_max]=size(I);
    M=size(mask,1);
    
    mask2=zeros(x_max,y_max);
    
% coordinates of the top and bottom corners
    m_cnt_y=m_cnt(1);
    m_cnt_x=m_cnt(2);    
    
% top left corner
    x1=m_cnt_x-M/2;
    y1=m_cnt_y-M/2;

% bottom right corner
    x4=m_cnt_x+M/2-1;
    y4=m_cnt_y+M/2-1;

    mask2(y1:y4,x1:x4)=mask;
    mask=mask2;
end    
    phi=bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;
 
end  