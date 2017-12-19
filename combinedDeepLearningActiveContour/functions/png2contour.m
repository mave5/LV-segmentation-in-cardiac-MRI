
function [contour,I1]=png2contour(I)

A=rgb2gray(I);
[x1,y1]=size(A);
I1=A(:,501:end-502);

[x,y]=size(I1);
count=0;
for i=1:x
    for j=1:y
        if I1(i,j)==150
            count=count+1;
            new_x(count)=j;
            new_y(count)=i;
        end
    end
end

contour(:,1)=new_x*256/x1;
contour(:,2)=new_y*256/x1;

end