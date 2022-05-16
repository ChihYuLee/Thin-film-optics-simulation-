%% THIS PROGRAM CALCULATES A, R, T, AND COLOR GAMUT FOR BOTH PCM PHASES OF A STACK OF LAYERS IN SWEEP THICKNESS
%% HOW TO USE:
% 1. Input the parameters in the block 1.
% 2. Be sure that the refractive index of the material you want to use is
%    included in the file import_refr_in. 
% 3. Enter the refractive index of each layer in order (in block #3). First entry
%    corresponds to first layer.
% 4. Input the thickness of each layer in block 4. The number of entries must be
%    consistent with the number of layers in parameters.
% 5. Run.
%% By Carlos Rios, August 2015

clear all
%close all
clc

%% Revised By Yi-Siou Huang, January 2022: Subplot
% Figure(1) Both phases A, R and T in sweep thickness 
% Figure(2) Both phases color gamut in sweep thickness

%% 1. PARAMETERS
theta_inc=0;        %% Incident angle in degrees
number_layers=5;    %% Total number of layers (without counting air)
                    %ATTENTION: number_layers must be consistent with the # of entries in d and n,k.
global min_lambda   % Just to be sure min_lambda is taken inside import_refr_in
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

A=zeros(500,5);
R=zeros(500,5);
T=zeros(500,5);
Phi=zeros(500,5);

dl_min=50;
dl_max=150;
dl_step=10;


%% 4. LAYER THICKNESSES first corresponds to upper layer
% thickITO=dl;
% thickPCM=15;

% d=1e-9*[thickPCM,30,thickITO,100,1000];

% This section is combined into section 6

%% 5. TRANSFER MATRIX CALCULATIONS
% theta0=theta_inc*pi/180;  % incident angle in radians
% [X,A(:,dl/dl_step),R(:,dl/dl_step),T(:,dl/dl_step),Phi(:,dl/dl_step)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n, k, d, pol);;
% [X,Aam(:,dl/dl_step),Ram(:,dl/dl_step),Tam(:,dl/dl_step),Phiam(:,dl/dl_step)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n2, k2, d, pol);

% This section is combined into section 6

%% 6. SUBPLOT
figure(1)
tiledlayout('flow')         


nexttile
for dl=dl_min:dl_step:dl_max
thickITO=dl;
thickPCM=15;
%d=[thickPCM, Ag, thickITO, Ag, ones(500)]
d=1e-9*[thickPCM,30,thickITO,100,1000];
theta0=theta_inc*pi/180;  % incident angle in radians
[X,A(:,dl/dl_step),R(:,dl/dl_step),T(:,dl/dl_step),Phi(:,dl/dl_step)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n, k, d, pol);;
[X,Aam(:,dl/dl_step),Ram(:,dl/dl_step),Tam(:,dl/dl_step),Phiam(:,dl/dl_step)]= transfer_matrix(theta0, n0, max_lambda, min_lambda, number_layers, n2, k2, d, pol);
plot(X,R(:,dl/dl_step),'-','DisplayName',['R=',mat2str(dl),'nm'])
title('Cry Reflection vs Wavelength')
hold on
xlabel('Wavelength (nm)') 
%legend(['R=',mat2str(dl),'nm']) 
axis([min(X) max(X) 0 1])
legend show
end
         
nexttile
for dl=dl_min:dl_step:dl_max
plot(X,T(:,dl/dl_step),'-','DisplayName',['T=',mat2str(dl),'nm'])
title('Cry Transimission vs Wavelength')
hold on
xlabel('Wavelength nm') 
% legend(['T=',mat2str(dl),'nm']) 
axis([min(X) max(X) 0 1])
legend show
 end
         
nexttile
for dl=dl_min:dl_step:dl_max
plot(X,A(:,dl/dl_step),'-','DisplayName',['A=',mat2str(dl),'nm'])
title('Cry Absorption vs Wavelength')
hold on
xlabel('Wavelength nm') 
% legend(['A=',mat2str(dl),'nm']) 
axis([min(X) max(X) 0 1])
legend show
end
         
         
nexttile
for dl=dl_min:dl_step:dl_max
plot(X,Ram(:,dl/dl_step),'-','DisplayName',['Ram=',mat2str(dl),'nm'])
title('Am Reflection vs Wavelength')
hold on
xlabel('Wavelength nm') 
%legend(['Ram=',mat2str(dl),'nm']) 
axis([min(X) max(X) 0 1])
legend show
end
         
nexttile
for dl=dl_min:dl_step:dl_max
plot(X,Tam(:,dl/dl_step),'-','DisplayName',['Tam=',mat2str(dl),'nm'])
title('Am Transimission vs Wavelength')
hold on
xlabel('Wavelength nm') 
% legend(['Tam=',mat2str(dl),'nm']) 
axis([min(X) max(X) 0 1])
legend show
end
         
nexttile
for dl=dl_min:dl_step:dl_max
plot(X,Aam(:,dl/dl_step),'-','DisplayName',['Aam=',mat2str(dl),'nm'])
title('Am Absorption vs Wavelength')
hold on
xlabel('Wavelength nm') 
% legend(['Aam=',mat2str(dl),'nm']) 
axis([min(X) max(X) 0 1])
legend show
end
         
numITO=num2str(thickITO);
numPCM=num2str(thickPCM);
mkdir('results_FROC');
save(['results_FROC/ITO_',numITO,'_PCM',numPCM,'nm_AM.txt'],'Ram','-ascii');
save(['results_FROC/ITO_',numITO,'_PCM',numPCM,'nm_CRY.txt'],'R','-ascii');

%--------------------------------------------------------------------------------------------------------------------

% COLOR 
% This section contains the color plots continued from Section 6: SUBPLOT
figure(2)
tiledlayout(1,2)

addpath '.\Color_map\github_repo\tbx'
addpath '.\Color_map\github_repo\misc'

nexttile
plotChromaticity
for dl=dl_min:dl_step:dl_max       

xyzAMO = rspd2xyz(X,Ram(:,dl/dl_step));
  
XYcieAMO(1) = xyzAMO(1) / (xyzAMO(1)+xyzAMO(2)+xyzAMO(3)); %this is X
XYcieAMO(2) = xyzAMO(2) / (xyzAMO(1)+xyzAMO(2)+xyzAMO(3)); %this is Y

% colorTable=[( XYcieAMO(1) XYcieAMO(2)),( XYcieCRY(1) XYcieCRY(2))];
hold on
title('Amorphous')
scatter( XYcieAMO(1), XYcieAMO(2),10,'black');
% text (XYcieAMO(1), XYcieAMO(2),'Am')
 
end

nexttile
plotChromaticity
for dl=dl_min:dl_step:dl_max       
 
xyzCRY = rspd2xyz(X,R(:,dl/dl_step));
 
XYcieCRY(1) = xyzCRY(1) / (xyzCRY(1)+xyzCRY(2)+xyzCRY(3));
XYcieCRY(2) = xyzCRY(2) / (xyzCRY(1)+xyzCRY(2)+xyzCRY(3));

% colorTable=[( XYcieAMO(1) XYcieAMO(2)),( XYcieCRY(1) XYcieCRY(2))];
hold on
title('Crystal')
scatter( XYcieCRY(1), XYcieCRY(2),10,'red'); 
% text (XYcieCRY(1), XYcieCRY(2),'Cry')
 
end