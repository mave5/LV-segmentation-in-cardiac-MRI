% compute gradient of internal energy function and edge factor 
% note that gradient of internal energy is equal to contour curvature

function [curvature EdgeTerm] = get_curvature_edge(phi,idx,g,gx,gy,gz)
    [dimy, dimx, dimz] = size(phi);        
    [y x z] = ind2sub([dimy,dimx,dimz],idx);  % get subscripts

    %-- get subscripts of neighbors
    ym1 = y-1; xm1 = x-1; zm1= z-1 ; 
    yp1 = y+1; xp1 = x+1; zp1= z+1 ; 
   
    %-- bounds checking  
    ym1(ym1<1) = 1; xm1(xm1<1) = 1; zm1(zm1<1)=1 ;               
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx; zp1(zp1>dimz) = dimz ;     

    %%-- get indexes for 8 neighbors
    %up
    idup = sub2ind(size(phi),yp1,x,z);
    idupf = sub2ind(size(phi),yp1,x,zp1);   % forward
    idupb = sub2ind(size(phi),yp1,x,zm1);  % behind
    %down
    iddn = sub2ind(size(phi),ym1,x,z);
    iddnf = sub2ind(size(phi),ym1,x,zp1);   % forward
    iddnb = sub2ind(size(phi),ym1,x,zm1);  % behind
    %left    
    idlt = sub2ind(size(phi),y,xm1,z);
    idltf = sub2ind(size(phi),y,xm1,zp1);   % forward
    idltb = sub2ind(size(phi),y,xm1,zm1);  % behind
    % right    
    idrt = sub2ind(size(phi),y,xp1,z);
    idrtf = sub2ind(size(phi),y,xp1,zp1);   % forward
    idrtb = sub2ind(size(phi),y,xp1,zm1);  % behind
    % up left
    idul = sub2ind(size(phi),yp1,xm1,z);
    idulf = sub2ind(size(phi),yp1,xm1,zp1);   % forward
    idulb = sub2ind(size(phi),yp1,xm1,zm1);  % behind
    %up right
    idur = sub2ind(size(phi),yp1,xp1,z);
    idurf = sub2ind(size(phi),yp1,xp1,zp1);   % forward
    idurb = sub2ind(size(phi),yp1,xp1,zm1);  % behind
    %down left
    iddl = sub2ind(size(phi),ym1,xm1,z);
    iddlf = sub2ind(size(phi),ym1,xm1,zp1);   % forward
    iddlb = sub2ind(size(phi),ym1,xm1,zm1);  % behind
    %down right
    iddr = sub2ind(size(phi),ym1,xp1,z);
    iddrf = sub2ind(size(phi),ym1,xp1,zp1);   % forward
    iddrb = sub2ind(size(phi),ym1,xp1,zm1);  % behind
    % front 
    idfr = sub2ind(size(phi),y,x,zp1);  % 
    % back 
    idbk = sub2ind(size(phi),y,x,zm1);  %   
    
    
    %-- get central derivatives of SDF at x,y,z
    dxplus =  -phi(idx) + phi(idrt) ; 
    dyplus = -phi(idx) + phi(idup) ; 
    dzplus = -phi(idx) + phi(idfr) ; 
    
    dxminus = +phi(idx) - phi(idlt) ; 
    dyminus = +phi(idx) - phi(iddn) ; 
    dzminus = +phi(idx) - phi(idbk) ; 
    
    dx = 1/2*( -phi(idlt)+phi(idrt) );
    dy = 1/2*( -phi(iddn)+phi(idup) );
    dz = 1/2*( -phi(idbk)+phi(idfr) );

    dxplusy = 1/2*(phi(idur) - phi(idul));
    dxminusy = 1/2*(phi(iddr) - phi(iddl));
    
    dxplusz = 1/2*(phi(idrtf) - phi(idltf)); 
    dxminusz = 1/2*(phi(idrtb) - phi(idltb)); 
    
    dyplusx = 1/2*(phi(idur) - phi(iddr)); 
    dyminusx = 1/2*(phi(idul) - phi(iddl));
    
    dyplusz = 1/2*(phi(idupf) - phi(iddnf)); 
    dyminusz = 1/2*(phi(idupb) - phi(iddnb));
    
    dzplusx = 1/2*(phi(idrtf) - phi(idrtb)); 
    dzminusx= 1/2*(phi(idltf) - phi(idltb)); 
    
    dzplusy = 1/2*(phi(idupf) - phi(idupb)); 
    dzminusy = 1/2*(phi(iddnf) - phi(iddnb)); 
    
    %--- calculating normals. 
    nplusx = dxplus./sqrt ( eps+(dxplus.^2 )+ ... 
        (( dyplusx+dy ) / 2 ).^2 +  (( dzplusx+dz ) / 2 ).^2) ;
    
    nplusy = dyplus./sqrt ( eps+(dyplus.^2 )+ ... 
        (( dxplusy+dx ) / 2 ).^2 + (( dzplusy+dz ) / 2 ).^2 ) ;
    nplusz = dzplus./sqrt ( eps+(dzplus.^2 )+ ... 
        (( dxplusz+dx ) / 2 ).^2 + (( dyplusz+dy ) / 2 ).^2 ) ;
    
    
    
    nminusx = dxminus./sqrt ( eps+(dxminus.^2 )+ ... 
        (( dyminusx+dy ) / 2 ).^2 +  (( dzminusx+dz ) / 2 ).^2) ;
    
    nminusy = dyplus./sqrt ( eps+(dyminus.^2 )+ ... 
        (( dxminusy+dx ) / 2 ).^2 + (( dzminusy+dz ) / 2 ).^2 ) ;
    nminusz = dzminus./sqrt ( eps+(dzminus.^2 )+ ... 
        (( dxminusz+dx ) / 2 ).^2 + (( dyminusz+dy ) / 2 ).^2 ) ;
    
    %-- calculating curvature
    curvature=((nplusx- nminusx)+(nplusy - nminusy) +(nplusz - nminusz) )/2 ;
    
    
    %-- calculate the edge term    
    phi_x = 1/2*(phi(idrt) - phi(idlt)) ;
    phi_y = 1/2*(phi(idup) - phi(iddn)) ;
    phi_z = 1/2*(phi(idfr) - phi(idbk)) ;
    
    s = sqrt(phi_x.^2 + phi_y.^2 + phi_z.^2); 
    Nx=phi_x./(s+eps); % add a small positive number to avoid division by zero
    Ny=phi_y./(s+eps);  
    Nz=phi_z./(s+eps);
   %curvature = div(Nx,Ny);
   %curvature2 = get_curvature(phi,idx);  % force from curvature penalty

    EdgeTerm = gx.*Nx+gy.*Ny + gz.*Nz+ g.*curvature ; %(idx);

    
