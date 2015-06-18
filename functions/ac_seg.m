%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I : input image
% init_mask: initial mask
% max_its : maximum iterations
% lengthEweight : weight of the length energy 
% shapeEweight  : weight of the shape energy
% display : display 

% This code was edited for LV segmentation in MRI from the original code by Shawn Lankton

function [seg,phi] = ac_seg(I,init_mask,max_its,lengthEweight,shapeEweight,display)

  %-- default behavior is to display intermediate outputs
  if(~exist('display','var'))
    display = true;
  end
  
  %-- ensures image is 2D double matrix
  I = im2graydouble(I);    
  
  %-- Create a signed distance map (SDF) from mask
  phi = mask2phi(init_mask);
  
  phi_prior=phi;
  
  %--main loop
  for its = 1:max_its   % Note: no automatic convergence test

    idx = find(phi <= 1.2 & phi >= -1.2);  %get the curve's narrow band
    
    %-- find interior and exterior mean
    upts = find(phi<=0);                 % interior points
    vpts = find(phi>0);                  % exterior points
    u = sum(I(upts))/(length(upts)+eps); % interior mean
    v = sum(I(vpts))/(length(vpts)+eps); % exterior mean
    
    
    F = (I(idx)-u).^2-(I(idx)-v).^2; % region-based force from image information
    curvature = get_curvature(phi,idx);  % force from curvature penalty
    
    % derivite of energy function for prior shape
    dEdl= -(phi(idx)-phi_prior(idx));
    
    % note that d phi/dt= - d E/d phi
    dphidt = F./max(abs(F)) + lengthEweight*curvature+shapeEweight*dEdl;  % gradient descent to minimize energy
    
    %-- maintain the CFL condition
    dt = .45/(max(dphidt)+eps);
    
    %-- evolve the curve
    phi(idx) = phi(idx) + dt*dphidt;

    indicator=0;
     
    %-- Keep SDF smooth
    phi = sussman(phi, .5);

    %-- intermediate output
    if((display>0)&&(mod(its,20) == 0)) 
      showCurveAndPhi(I,phi,its);  
    end
    if indicator==1
        disp(['stopped at:',num2str(its)])
        break
    end
  end
  
  %-- final output
  if(display)
    showCurveAndPhi(I,phi,its);
  end
  
  %-- make mask from SDF
  seg = phi<=0; %-- Get mask from levelset

  
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
  
  
%-- Displays the image with curve superimposed
function showCurveAndPhi(I, phi, i)
   max_range=min(200,max(I(:)));
  imshow(I,'initialmagnification',200,'displayrange',[0 max_range]); hold on;
  contour(phi, [0 0], 'g','LineWidth',4);
  contour(phi, [0 0], 'k','LineWidth',2);
  hold off; title([num2str(i) ' Iterations']); drawnow;
  
%-- converts a mask to a SDF
function phi = mask2phi(init_a)
  phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
  
%-- compute curvature along SDF
function curvature = get_curvature(phi,idx)
    [dimy, dimx] = size(phi);        
    [y x] = ind2sub([dimy,dimx],idx);  % get subscripts

    %-- get subscripts of neighbors
    ym1 = y-1; xm1 = x-1; yp1 = y+1; xp1 = x+1;

    %-- bounds checking  
    ym1(ym1<1) = 1; xm1(xm1<1) = 1;              
    yp1(yp1>dimy)=dimy; xp1(xp1>dimx) = dimx;    

    %-- get indexes for 8 neighbors
    idup = sub2ind(size(phi),yp1,x);    
    iddn = sub2ind(size(phi),ym1,x);
    idlt = sub2ind(size(phi),y,xm1);
    idrt = sub2ind(size(phi),y,xp1);
    idul = sub2ind(size(phi),yp1,xm1);
    idur = sub2ind(size(phi),yp1,xp1);
    iddl = sub2ind(size(phi),ym1,xm1);
    iddr = sub2ind(size(phi),ym1,xp1);
    
    %-- get central derivatives of SDF at x,y
    phi_x  = -phi(idlt)+phi(idrt);
    phi_y  = -phi(iddn)+phi(idup);
    phi_xx = phi(idlt)-2*phi(idx)+phi(idrt);
    phi_yy = phi(iddn)-2*phi(idx)+phi(idup);
    phi_xy = -0.25*phi(iddl)-0.25*phi(idur)...
             +0.25*phi(iddr)+0.25*phi(idul);
    phi_x2 = phi_x.^2;
    phi_y2 = phi_y.^2;
    
    %-- compute curvature (Kappa)
    curvature = ((phi_x2.*phi_yy + phi_y2.*phi_xx - 2*phi_x.*phi_y.*phi_xy)./...
              (phi_x2 + phi_y2 +eps).^(3/2)).*(phi_x2 + phi_y2).^(1/2);        
  
%-- Converts image to one channel (grayscale) double
function img = im2graydouble(img)    
  [dimy, dimx, c] = size(img);
  if(isfloat(img)) % image is a double
    if(c==3) 
      img = rgb2gray(uint8(img)); 
    end
  else           % image is a int
    if(c==3) 
      img = rgb2gray(img); 
    end
    img = double(img);
  end

%-- level set re-initialization by the sussman method
function D = sussman(D, dt)
  % forward/backward differences
  a = D - shiftR(D); % backward
  b = shiftL(D) - D; % forward
  c = D - shiftD(D); % backward
  d = shiftU(D) - D; % forward
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;
  
%-- whole matrix derivatives
function shift = shiftD(M)
  shift = shiftR(M')';

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];

function shift = shiftR(M)
  shift = [ M(:,1) M(:,1:size(M,2)-1) ];

function shift = shiftU(M)
  shift = shiftL(M')';
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    

  




