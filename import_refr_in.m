%function import_refr_in(min_lambda)

global n_GSTcry n_GSTam k_GSTcry k_GSTam n_AISTcry n_AISTam k_AISTcry k_AISTam n_SbSecry n_SbSeam k_SbSecry k_SbSeam
global n_ITO k_ITO n_Pt k_Pt n_Ag n_Ag2 k_Ag k_Ag2 n_SiO k_SiO n_Ge k_Ge n_TiO2 k_TiO2

%% GST REFRACTIVE INDEX
ref_ind_GST=dlmread('refractive_index_GST_pei.txt');

% This is a 5 column matrix. 1. Wavelength (from 250:10:1750). 2. n_cry 3. k_cry 4. n_am 5.k_am
lambda_GST=ref_ind_GST(:,1);
n_GSTcry=ref_ind_GST(min_lambda-299:500,2);
k_GSTcry=ref_ind_GST(min_lambda-299:500,3);
n_GSTam=ref_ind_GST(min_lambda-299:500,4);
k_GSTam=ref_ind_GST(min_lambda-299:500,5);
%% If a small step of wavelengths is required use the following
% interpolations for getting the refractive index per 1nm
% n_am=interp1(lambda_GST,n_am1,[min_lambda:1:1745],'linear','extrap');
% n_cry=interp1(lambda_GST,n_cry1,[min_lambda:1:1745],'linear','extrap');
% k_am=interp1(lambda_GST,k_am1,[min_lambda:1:1745],'linear','extrap');
% k_cry=interp1(lambda_GST,k_cry1,[min_lambda:1:1745],'linear','extrap');
% figure, plot((300:10:1750),n_am1,'o',(300:10:1750),k_am1,'o',(300:10:1750),n_cry1,'o',(300:10:1750),k_cry1,'o',(300:1:1745),n_am,'*',(300:1:1745),k_am,'*',(300:1:1745),n_cry,'*',(300:1:1745),k_cry,'*');

clear ref_ind_GST

%% AIST REFRACTIVE INDEX
ref_ind_AIST=dlmread('refractive_index_AIST.txt');
% This is a 5 column matrix. 1. Wavelength (from 250:10:1750). 2. n_cry 3. k_cry 4. n_am 5.k_am
lambda_AIST=ref_ind_AIST(:,1);
n_am1=ref_ind_AIST(:,2);
k_am1=ref_ind_AIST(:,3);
n_cry1=ref_ind_AIST(:,4);
k_cry1=ref_ind_AIST(:,5);
%If a small step of wavelengths is required use the following
% interpolations for getting the refractive index per 1nm
n_AISTam=interp1(lambda_AIST,n_am1,min_lambda:1:1745,'linear','extrap');
n_AISTcry=interp1(lambda_AIST,n_cry1,min_lambda:1:1745,'linear','extrap');
k_AISTam=interp1(lambda_AIST,k_am1,min_lambda:1:1745,'linear','extrap');
k_AISTcry=interp1(lambda_AIST,k_cry1,min_lambda:1:1745,'linear','extrap');
% figure, plot((300:10:1750),n_am1,'o',(300:10:1750),k_am1,'o',(300:10:1750),n_cry1,'o',(300:10:1750),k_cry1,'o',(300:1:1745),n_am,'*',(300:1:1745),k_am,'*',(300:1:1745),n_cry,'*',(300:1:1745),k_cry,'*');
clear n_cry1 k_cry1 n_am1 k_am1 ref_ind_AIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ITO REFRACTIVE INDEX
ref_ind_ITO=dlmread('refractive_index_ITOam.txt');
lambda_ITO=ref_ind_ITO(:,1);
n_ITO=ref_ind_ITO(min_lambda-299:501,2);
k_ITO=ref_ind_ITO(min_lambda-299:501,3);
% figure, plot(lambda_ITO,n_ITO, lambda_ITO,k_ITO)

%% Pt REFRACTIVE INDEX (Just to compare with Pei)
ref_ind_Pt=dlmread('Pt.txt');
lambda_Pt=ref_ind_Pt(:,1);
n_Pt1=ref_ind_Pt(:,2);
k_Pt1=ref_ind_Pt(:,3);

n_Pt=interp1(lambda_Pt,n_Pt1,[min_lambda:1:800],'linear','extrap');
k_Pt=interp1(lambda_Pt,k_Pt1,[min_lambda:1:800],'linear','extrap');
% figure, plot((300:1:800),n_Pt,'o',(300:1:800),k_Pt,'o');

%% Ag REFRACTIVE INDEX (Just to compare with Pei)
ref_ind_Ag=dlmread('refractive_index_Ag_Ellipsometry_meausred.txt');
lambda_Ag=ref_ind_Ag(:,1);
n_Ag=ref_ind_Ag(min_lambda-299:501,2);
k_Ag=ref_ind_Ag(min_lambda-299:501,3);

%% Ag REFRACTIVE INDEX in a longer spectrum
ref_ind_Ag=dlmread('refractive_index_Ag_long.txt');
lambda_Ag2=ref_ind_Ag(:,1);
n_Ag2=ref_ind_Ag(:,2);
k_Ag2=ref_ind_Ag(:,3);

n_Ag2=interp1(lambda_Ag2,n_Ag2,[min_lambda:1:10000],'linear','extrap');
k_Ag2=interp1(lambda_Ag2,k_Ag2,[min_lambda:1:10000],'linear','extrap');

%% SiO REFRACTIVE INDEX (Just to compare with Pei)
ref_ind_SiO=dlmread('refractive_index_SiO.txt');
lambda_SiO=ref_ind_SiO(:,1);
n_SiO=ref_ind_SiO(min_lambda-299:501,2);
k_SiO=ref_ind_SiO(min_lambda-299:501,3);

% n_Pt=interp1(lambda_Pt,n_Pt1,[min_lambda:1:800],'linear','extrap');
% k_Pt=interp1(lambda_Pt,k_Pt1,[min_lambda:1:800],'linear','extrap');
% figure, plot((300:1:800),n_Pt,'o',(300:1:800),k_Pt,'o');
%end

%% Ge REFRACTIVE INDEX
ref_ind_Ge=dlmread('refractive_index_Ge.txt');
lambda_Ge=ref_ind_Ge(:,1);
n_Ge=ref_ind_Ge(min_lambda-299:501,2);
k_Ge=ref_ind_Ge(min_lambda-299:501,3);

%% TiO2 REFRACTIVE INDEX
ref_ind_TiO2=dlmread('refractive_index_TiO2.txt');
lambda_TiO2=ref_ind_TiO2(:,1);
n_TiO2=ref_ind_TiO2(min_lambda-299:501,2);
k_TiO2=ref_ind_TiO2(min_lambda-299:501,3);

%% SbSe REFRACTIVE INDEX
ref_ind_SbSe=dlmread('refractive_index_SbSe_Ellipsometry_meausred_10nm.txt');

% This is a 5 column matrix. 1. Wavelength (from 300:10:800). 2. n_cry 3. k_cry 4. n_am 5. k_am
lambda_SbSe=ref_ind_SbSe(:,1);
n_SbSecry=ref_ind_SbSe(min_lambda-299:500,2);
k_SbSecry=ref_ind_SbSe(min_lambda-299:500,3);
n_SbSeam=ref_ind_SbSe(min_lambda-299:500,4);
k_SbSeam=ref_ind_SbSe(min_lambda-299:500,5);

