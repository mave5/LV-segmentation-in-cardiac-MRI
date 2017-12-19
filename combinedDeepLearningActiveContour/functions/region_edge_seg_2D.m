

function [seg,phi] = region_edge_seg_2D(I,phi,max_its,option,display)
  
% Initialization %-- default value for parameter alpha is .1
if(~exist('display','var'))
    display = true;
end

% load parameters of Deep Learning Netowrks
load DBNparams.mat;
Mroi=100;

%  Edge detector parameters. 
sigma = .8;   
EdgeFactor = 1;

% convolution of image with Gauusian kernel
G=fspecial('gaussian',15,sigma);
D_smooth=conv2(I,G,'same');  

% gradient of image
[Dx,Dy]=gradient(D_smooth);
f=Dx.^2+Dy.^2;

% Edge function : g(I)=1/(1+Delta(G*I))
Edge_fun=1./(1+EdgeFactor*f);

% combiningg weights
Edgeweight = option.Edgeweight ;   
Regionweight = option.Regionweight;
interalWeight=option.CurvatureWeight;
dlnWeight1=option.DLNWeight1;
dlnWeightN=option.DLNWeightN;

%-- ensures image is 2D double matrix
I = im2graydouble(I);    

phi_LV=phi;

% compute gradient of the edge-factor function g(I)=1/(1+G*I)
[Edge_fun_x, Edge_fun_y] = gradient(Edge_fun); 
  
  %% %--main loop
  for its = 1:max_its   % Note: no automatic convergence test

    % find sub-image using mask from previous iteration
    mask_tm1 = phi<=0 ;
    [subI,m_cnt]=mask2subImage(I,mask_tm1,Mroi);
      
    % run Deep Learning network to find the segmentation of LV
    %yLV=DLN(subI,stackedAEOptTheta,inputSize,hiddenSizeL1,hiddenSizeL2,outputSize,netconfig);
        
    % find SDF from the LV mask
    %phi_LV=mask2phi(yLV,m_cnt,I);
    
    % weights 
    dlnWeightn(its)=dlnWeight1+(dlnWeightN-dlnWeight1)/cosh(10*(its/max_its-1));
   
    
   %get the curve's narrow band
    idx = find(phi <= 1.2 & phi >= -1.2);  
  % -- exit if there is no zero level. 
    if (sum(idx) == 0) 
        break  
    else 
  %-- find interior and exterior mean
    int_pts = find(phi<=0);% interior points
    ext_pts = find(phi>0); % exterior points
    int_mean = sum(I(int_pts))/(length(int_pts)+eps); % interior mean
    ext_mean = sum(I(ext_pts))/(length(ext_pts)+eps); % exterior mean
    
    % gradient of region-based term: in external energy function
    dF_region = (I(idx)-int_mean).^2-(I(idx)-ext_mean).^2;         
    dF_region = dF_region./max(abs(dF_region)); %normalize
    
    % gradient of edge term and internal energy function
    [dF_curvature, dF_edge] = get_curvature_edge(phi,idx,Edge_fun(idx),Edge_fun_x(idx),Edge_fun_y(idx),0) ;
        
    
    % gradient of total external energy function 
    dE_ext = Regionweight*dF_region + Edgeweight*dF_edge ; % - PriorShapeWeight* DiffShape % gradient descent to minimize energy
    
    % dervitive of engergy of DL-LV
    dE_lv=2*((-phi_LV(idx)+phi(idx)));

    % gradient of total energy 
    dE_tot=dE_ext+ + interalWeight*dF_curvature+dlnWeightn(its)*dE_lv;
    
    %-- maintain the CFL condition
    dt = .45/(max(dE_tot)+eps);
        
    %-- evolve the curve
    phi(idx) = phi(idx) + dt.*dE_tot;

    
    
    %-- Keep SDF smooth 
    if (mod(its,4) ==1 ) 
         phi = sussman(phi, .5);
    end
    
    %-- make mask from SDF
    seg(:,:,its) = phi<=0 ;   %-- Get mask from levelset
   % imshow( seg(:,:,its));
   % profile viewer 
     its
    end
  end
  

 
%% ---------------------------------------------------------------------
%---------------------------------------------------------------------
%-- AUXILIARY FUNCTIONS ----------------------------------------------
%---------------------------------------------------------------------
%---------------------------------------------------------------------
   
%-- converts a mask to a SDF
%function phi = mask2phi(init_a)
%  phi=bwdist(init_a)-bwdist(1-init_a)+im2double(init_a)-.5;
  
function DiffShape = get_DiffShape(phi,thai,idx) 
    Hphi = zeros(size(phi)); 
    ppts = find(phi>0);   % exterior  where H(phi) = 1 
    Hphi(ppts) = 1 ; 
    
    Hthai = zeros(size(thai)); 
    tpts = find(thai>0) ;   % exterior where H(thai) =1 
    Hthai(tpts) = 1 ; 
    DiffShape = 2*(Hphi(idx) - Hthai(idx)) ;   % a 3D matrix in this case. 
%-- compute curvature along SDF

%-- Converts image to one channel (grayscale) double
function img = im2graydouble(img)    
  [dimy, dimx, dimz, c] = size(img);
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
  e = D - shiftF(D); 
  f = shiftB(D) - D ; 
  
  a_p = a;  a_n = a; % a+ and a-
  b_p = b;  b_n = b;
  c_p = c;  c_n = c;
  d_p = d;  d_n = d;
  e_p = e ; e_n = e;
  f_p = f;  f_n = f; 
  
  
  a_p(a < 0) = 0;
  a_n(a > 0) = 0;
  b_p(b < 0) = 0;
  b_n(b > 0) = 0;
  c_p(c < 0) = 0;
  c_n(c > 0) = 0;
  d_p(d < 0) = 0;
  d_n(d > 0) = 0;
  e_p(e < 0) = 0;
  e_n(e > 0) = 0; 
  f_p(f < 0) = 0; 
  f_n(f > 0) = 0 ; 
  
  dD = zeros(size(D));
  D_neg_ind = find(D < 0);
  D_pos_ind = find(D > 0);
  dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                       + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2) ...
                       + max(e_p(D_pos_ind).^2, f_n(D_pos_ind).^2))- 1;
  dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                       + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2) ...
                       + max(e_n(D_neg_ind).^2, f_p(D_neg_ind).^2)) - 1;
  
  D = D - dt .* sussman_sign(D) .* dD;  % supposedly this procedure must be 
  % continuting until convergence but it seems only one iteration is
  % enough.
  
  
  
  %% 
  function EdgeTerm =  get_EdgeTerm(phi,g,gx,gy,gz,curvature,idx) 
   
 [phi_x,phi_y,phi_z] = gradient(phi);
%    phi_x = phi - shiftL(phi); 
%    phi_y = phi - shiftU(phi); 
%    phi_z = phi - shiftB(phi); 
   
   
   s=sqrt(phi_x(idx).^2 + phi_y(idx).^2 + phi_z(idx).^2);
   Nx=phi_x(idx)./(s+eps); % add a small positive number to avoid division by zero
   Ny=phi_y(idx)./(s+eps);  
   Nz=phi_z(idx)./(s+eps);
   %curvature = div(Nx,Ny);
   curvature2 = get_curvature(phi,idx);  % force from curvature penalty

   EdgeTerm = gx.*Nx+gy.*Ny + gz.*Nz+ g.*curvature ; %(idx);
  
    

%%-- whole matrix derivatives
function shift = shiftD(M)
  %shift = shiftR(M')';
  shift = [ M(1,:,:) ; M(1: size(M,1)-1,:,:) ];

function shift = shiftL(M)
  shift = [ M(:,2:size(M,2),:) M(:,size(M,2),:) ];

function shift = shiftR(M)
  shift = [ M(:,1,:) M(:,1:size(M,2)-1,:) ];

function shift = shiftU(M)
  %shift = shiftL(M')';
  shift = [ M(2:size(M,1),:,:) ; M(size(M,1),:,:) ];
  
function shift = shiftF(M)
  %shift = shiftL(M')';
  shift = cat(3, M(:,:,1), M(:,:,1:size(M,3)-1));  
  
function shift = shiftB(M)
%shift = shiftL(M')';
shift = cat(3, M(:,:,2:size(M,3)), M(:,:,size(M,3)));  
  
function S = sussman_sign(D)
  S = D ./ sqrt(D.^2 + 1);    

  




