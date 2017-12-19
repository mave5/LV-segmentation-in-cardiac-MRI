function D_clipped = ImageClipping(D,zRange,max_res) 
iz = floor(mean (zRange)) ; % sample image
I = squeeze(D(:,:,iz)) ; 

disp('Next: clipping the image'); 
imshow(I,[0,max_res]); 
str = sprintf('MRI No.= %f',iz); % str = sprintf('Plot with frequency = %f , and wavelength = %f',number,number2);
title(str);
waitfor(msgbox(' please locate top left point'))
LT_point =  ginput(1); % left top
LT_point = floor(LT_point) ; 
waitfor(msgbox('located bottom right point')); 
RB_point = ginput(1);
RB_point = floor(RB_point) ; 
%I = I(LT_point(2): RB_point(2), LT_point(1): RB_point(1));
D_clipped = D(LT_point(2): RB_point(2), LT_point(1): RB_point(1),zRange);
