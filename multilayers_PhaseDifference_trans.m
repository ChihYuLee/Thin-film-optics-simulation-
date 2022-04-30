%% THIS PROGRAM CALCULATES A, R and T FOR BOTH SbSe PHASES, and PHASE DIFFERNCE of R(Phi) and T(Psi) OF A STACK OF LAYERS
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
%% Revised By Yi-Siou Huang, April 2022
clear all
close all
clc

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
n=struct('n_val',{n_SbSecry, n_Ag, n_ITO, n_Ag, ones(500)});
k=struct('k_val',{k_SbSecry, k_Ag, k_ITO, k_Ag, zeros(500)});

n2=struct('n_val',{n_SbSeam, n_Ag, n_ITO, n_Ag, ones(500)});
k2=struct('k_val',{k_SbSeam, k_Ag, k_ITO, k_Ag, zeros(500)});

%% 4. LAYER THICKNESSES first corresponds to upper layer
d=1e-9*[15,30,150,100,1000];
% d=1e-9*[20,30,120,25,1000]; good for T
% d=1e-9*[10,12,80,100,1000]; good for R

   
%% 5. TRANSFER MATRIX CALCULATIONS
% Phi is phase difference in reflection
% Psi is phase differene in transmission
theta0=theta_inc*pi/180;  % incident angle in radians
[X,A,R,T,Phi,Psi]= transfer_matrix_trans(theta0, n0, max_lambda, min_lambda, number_layers, n, k, d, pol);
[Xam,Aam,Ram,Tam,Phiam,Psiam]= transfer_matrix_trans(theta0, n0, max_lambda, min_lambda, number_layers, n2, k2, d, pol);

%% 6. PLOT
    % ALL DOWN HERE IS TO GENERATE THE PLOT, IT DOES NOTHING TO DO WITH THE
    % CALCULATION ITSELF
         figure, 
         p=plot(X,R,'-',X,T,'--', X,Ram,'-',X,Tam,'--',X,(Phi-Phiam)/(2*pi),'--',X, (Psi-Psiam)/(2*pi),'-')
         set(p,'LineWidth',2);
         xlabel('Wavelength (nm)') 
         legend('Rcry','Tcry','Ram','Tam','Phase diff of R','Phase diff of T') 
%          legend('R','T','A', 'Ram','Tam','Aam') 
         axis([min(X) max(X) -1 1])
         
         
        imax = find(max(R) == R);
        text(min_lambda,0,['Maximum R=  ',num2str(R(imax)),' @ ',num2str(X(imax)),'nm',],...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','left',...
        'FontSize',8)
        text(min_lambda,0,['Layers=',mat2str(d*1e9),' nm   ','Angle=', mat2str(theta0*180/pi)],...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left',...
        'FontSize',8);
        title(strcat(pol,' spectrum',' Inc: ', num2str(theta_inc),'^{o}'));

        
        %num=num2str(150*10^9);
        %saveas(gcf,['ITO_',num,'nm'],'bmp');
