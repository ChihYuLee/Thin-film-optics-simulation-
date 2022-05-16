%% THIS PROGRAM CALCULATES COLOR CONTRAST IN THE SWEEP THICKNESS OF PCM AND ITO
%% HOW TO USE:
% 1. Input the parameters in the block 1.
% 2. Be sure that the refractive index of the material you want to use is included in the file import_refr_in. 
% 3. Enter the refractive index of each layer in order (in block #3). First entry corresponds to first layer.
% 4. Input the sweep thickness of PCM and ITO layers(min:step:max) in block 4, other layers in block 6. 
%    The number of entries must be consistent with the number of layers in parameters.
% 5. Run.

%% By Carlos Rios, August 2015
clear all
close all
clc

%% Revised By Yi-Siou Huang, January 2022: Subplot
% Figure: Color contrast

%% 1. PARAMETERS
theta_inc=0;        %% Incident angle in degrees
number_layers=5;    %% Total number of layers (without counting air)
                    %  ATTENTION: number_layers must be consistent with the # of entries in d and n,k.
global min_lambda   %  Just to be sure min_lambda is taken inside import_refr_in
min_lambda=300;     %% Minimum lambda to run simulation
max_lambda=799;     %% Maximum lambda to run simulation / Beware: some materials n and k are given up to 800nm
n0=1;               %% refraction index of incident medium
pol='TE';           %% Choose either TE or TM polarization (relevant under oblique incidence)

%% 2. IMPORTS ALL THE REFRACTIVE INDEX FROM THE FILE import_refr_in
import_refr_in      
k_SiO2=zeros(1,1700);
n_SiO2=zeros(1,1700);
n_SiO2(1,:)=1.44;

%% 3. REFRACTIVE INDEX OF EACH LAYER first corresponds to upper layer
% This step defines which material goes in which layer
n_PCMcry=n_SbSecry;
k_PCMcry=k_SbSecry;
n_PCMam=n_SbSeam;
k_PCMam=k_SbSeam;

n=struct('n_val',{n_PCMcry, n_Ag, n_ITO, n_Ag, ones(500)});
k=struct('k_val',{k_PCMcry, k_Ag, k_ITO, k_Ag, zeros(500)});

n2=struct('n_val',{n_PCMam, n_Ag, n_ITO, n_Ag, ones(500)});
k2=struct('k_val',{k_PCMam, k_Ag, k_ITO, k_Ag, zeros(500)});

%% 4. LAYER THICKNESSES first corresponds to upper layer
% Please input the thickness of PCM and ITO layers here, other layers in block 6.: d=1e-9*[thickPCM,30,thickITO,100,1000] 

% thickPCM=tl;
tl_min=10;
tl_max=20;
tl_step=1;
a= [(tl_max- tl_min)/tl_step]+1 

% thickITO=dl;
dl_min=50;
dl_max=150;
dl_step=10;
b=[(dl_max- dl_min)/dl_step]+1

% d=1e-9*[thickPCM,30,thickITO,100,1000];

%% 5. TRANSFER MATRIX CALCULATIONS
% theta0=theta_inc*pi/180;  % incident angle in radians
% [X,A(:,tl/tl_step,dl/dl_step),R(:,tl/tl_step,dl/dl_step),T(:,tl/tl_step,dl/dl_step),Phi(:,tl/tl_step,dl/dl_step)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n, k, d, pol);
% [Xam,Aam(:,tl/tl_step,dl/dl_step),Ram(:,tl/tl_step,dl/dl_step),Tam(:,tl/tl_step,dl/dl_step),Phiam(:,tl/tl_step,dl/dl_step)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n2, k2, d, pol);

% This block is combined to the block 6
%%  6. Subplot
% Pease find more details in Carlos's Ph.D. thesis Section 6.1 or keywords: CIE 1931 color space

tiledlayout(2,2) 
figure(1)
   
addpath '.\Color_map\github_repo\tbx'
addpath '.\Color_map\github_repo\misc'
  
 
% Reflection
 nexttile
 c=[];
 tl=[tl_max:-tl_step:tl_min]; 
 dl=[dl_min:dl_step:dl_max]; 
 for i=1:a; 
 for j=1:b;

 thickPCM=tl(i);
 thickITO=dl(j);
 
 d=1e-9*[thickPCM,30,thickITO,25,1000]; % layer thickness
 theta0=theta_inc*pi/180;  % incident angle in radians
 [X,A(:,i,j),R(:,i,j),T(:,i,j),Phi(:,i,j)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n, k, d, pol);
 [Xam,Aam(:,i,j),Ram(:,i,j),Tam(:,i,j),Phiam(:,i,j)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n2, k2, d, pol);
 
 xyzAMO = rspd2xyz(X,Ram(:,i,j));
 
 XYcieAMO(1) = xyzAMO(1) / (xyzAMO(1)+xyzAMO(2)+xyzAMO(3)); % x = X / (X+ Y+ Z)
 XYcieAMO(2) = xyzAMO(2) / (xyzAMO(1)+xyzAMO(2)+xyzAMO(3)); % y = Y / (X+ Y+ Z)
 % X: A linear combination of cone-response curves chosen to be non-negative
 % Y: Defined as the luminance
 % Z: Approximatedvalue to the blue stimulation
 
 hold on 
 title('Am Reflection color')
 c=[c;xyzAMO/max(xyzAMO)]; 
 
 end
 end
 grid_size = [a b]; % https://blogs.mathworks.com/steve/2020/03/10/how-to-display-color-swatches/?doing_wp_cron=1641536777.3834769725799560546875
 daspect([1 1 1]); % Making the color swatches square
 gap=0;
 xlabel('ITO (nm)');
 ylabel('PCM (nm)');
 colorSwatches(c,grid_size,gap=0);


  
nexttile
 c=[];
 tl=[tl_max:-tl_step:tl_min]; 
 dl=[dl_min:dl_step:dl_max]; 
 for i=1:a; 
 for j=1:b;
 
 xyzCRY = rspd2xyz(X,R(:,i,j));
 
 XYcieCRY(1) = xyzCRY(1) / (xyzCRY(1)+xyzCRY(2)+xyzCRY(3)); % x = X / (X+ Y+ Z)
 XYcieCRY(2) = xyzCRY(2) / (xyzCRY(1)+xyzCRY(2)+xyzCRY(3)); % y = Y / (X+ Y+ Z) 
 
 hold on
 title('Cry Reflection color')
 c=[c;xyzCRY/max(xyzCRY)]; 
 
 end
 end
 grid_size = [a b];
 daspect([1 1 1]); % Making the color swatches square
 gap=0;
 xlabel('ITO (nm)');
 ylabel('PCM (nm)');
 colorSwatches(c,grid_size,gap=0);
 
 
% Transmission
 nexttile
 c=[];
 tl=[tl_max:-tl_step:tl_min]; 
 dl=[dl_min:dl_step:dl_max];
 for i=1:a; 
 for j=1:b;
     
 xyzAMO = rspd2xyz(X,Tam(:,i,j));
  
 XYcieAMO(1) = xyzAMO(1) / (xyzAMO(1)+xyzAMO(2)+xyzAMO(3)); % x = X / (X+ Y+ Z)
 XYcieAMO(2) = xyzAMO(2) / (xyzAMO(1)+xyzAMO(2)+xyzAMO(3)); % y = Y / (X+ Y+ Z)

 hold on
 title('Am Transmission color')
 c=[c;xyzAMO/max(xyzAMO)];  

 end
 end
 grid_size = [a b];
 daspect([1 1 1]); % Making the color swatches square
 gap=0;
 xlabel('ITO (nm)');
 ylabel('PCM (nm)');
 colorSwatches(c,grid_size,gap=0);
 
 
 nexttile
 c=[];
 tl=[tl_max:-tl_step:tl_min]; 
 dl=[dl_min:dl_step:dl_max]; 
 for i=1:a; 
 for j=1:b;
     
 xyzCRY = rspd2xyz(X,T(:,i,j));
 
 XYcieCRY(1) = xyzCRY(1) / (xyzCRY(1)+xyzCRY(2)+xyzCRY(3)); % x = X / (X+ Y+ Z)
 XYcieCRY(2) = xyzCRY(2) / (xyzCRY(1)+xyzCRY(2)+xyzCRY(3)); % y = Y / (X+ Y+ Z)
 
 hold on
 title('Cry Transmission color')
 c=[c;xyzCRY/max(xyzCRY)];

 end
 end
 grid_size = [a b];
 daspect([1 1 1]); % Making the color swatches square
 gap=0;
 xlabel('ITO (nm)');
 ylabel('PCM (nm)');
 colorSwatches(c,grid_size,gap=0);