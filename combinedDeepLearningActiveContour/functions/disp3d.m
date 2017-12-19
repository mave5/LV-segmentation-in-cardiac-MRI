
function disp3d(V,filter_size)

if nargin==1
    filter_size=9;
end

if filter_size==0
    D=V;
else
    D = smooth3(V,'box',filter_size);
    %D = smooth3(V,'gaussian',filter_size);
end


%patch(isocaps(V,.5),...
%  'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(D,.5),...
   'FaceColor','red','EdgeColor','none');

isonormals(D,p1)

view(35,30) 
axis vis3d tight
%camlight left; 
camlight right;
lighting phong

%set(gcf,'Renderer','zbuffer'); lighting phong
%set(p1,'SpecularColorReflectance',0,'SpecularExponent',50)

end