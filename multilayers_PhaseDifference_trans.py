#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 19:55:58 2022

@author: chihyulee

Translate "multilayers_PhaseDifference_trans.m" from Carlos group

Calculate A,R,T for both states; 
Calculate phase difference of R(Phi) and A(Psi) of different layers 
"""

import numpy as np 
import pandas as pd 
import warnings
import import_refr_in_transfer_matrix_trans as f
import matplotlib.pyplot as plt 
import timeit
import time
%matplotlib widget
import ipywidgets as widgets

 

############### Parameters ################
theta_inc=0       ##Incident angle in degrees
number_layers=5   ##Total number of layers (without counting air)
min_lambda=300    ##Minimum lambda to run simulation
max_lambda=799    ##Maximum lambda to run simulation
n0=1              ##refraction index of incident medium
pol='TE'          ##Choose either TE or TM polarization (relevant under oblique incidence)
PCM='SbSe'        ##Chose the PCM materials 
d_PCM= float(input('PCM thickness:')) #15    ##Thickness of PCM
d_Ag1= float(input('Ag1 thickness:')) #30     ##Thickness of Ag1
d_ITO= float(input('ITO thickness:')) #150
d_Ag2= float(input('Ag2 thickness:')) #100
d_air=1000        ##Thickness of Air 



############### Import all the refractive idexes from the files #############
start = timeit.default_timer()

# Order from top to bottom 
stack_layers=[PCM,'Ag','ITO','Ag','Air']
stack_thickness=[d_PCM, d_Ag1, d_ITO, d_Ag2, d_air] 
stack_thickness_nm=[i *1e-9 for i in stack_thickness]
n_crystalline=np.zeros((max_lambda-min_lambda+1, number_layers))  
k_crystalline=np.zeros((max_lambda-min_lambda+1, number_layers))  
n_amorphous=np.zeros((max_lambda-min_lambda+1, number_layers))    
k_amorphous=np.zeros((max_lambda-min_lambda+1, number_layers))    

for i,l in enumerate(stack_layers):
    
    #materials layer 
    if i is not len(stack_layers)-1:
        n_crystalline[:,i]= f.read_refractive_index_file(l,min_lambda)[1][:max_lambda-min_lambda+1,0]
        k_crystalline[:,i]= f.read_refractive_index_file(l,min_lambda)[1][:max_lambda-min_lambda+1,1]
        n_amorphous[:,i]= f.read_refractive_index_file(l,min_lambda)[1][:max_lambda-min_lambda+1,2]
        k_amorphous[:,i]= f.read_refractive_index_file(l,min_lambda)[1][:max_lambda-min_lambda+1,3]
   #Air layer
    else:
        n_crystalline[:,i]= np.ones(max_lambda-min_lambda+1)
        k_crystalline[:,i]= np.zeros(max_lambda-min_lambda+1)
        n_amorphous[:,i]= np.ones(max_lambda-min_lambda+1)
        k_amorphous[:,i]= np.zeros(max_lambda-min_lambda+1)



############## Transfer Matrix Calculations ##################

#incident angle in radians
theta0= theta_inc*np.pi/180 


#array containing X, A, R, T, Phi, Psi
results_c= f.transfer_matrix_trans(theta0, n0, max_lambda, min_lambda, number_layers, n_crystalline, k_crystalline, stack_thickness_nm, pol)
results_a= f.transfer_matrix_trans(theta0, n0, max_lambda, min_lambda, number_layers, n_amorphous, k_amorphous, stack_thickness_nm, pol)

############## Plot Response after Calculation ##############################

x= results_a[0, :]/1e-9
plt.figure(figsize = (10,6))
plt.xlabel('Wavelength(nm)')
plt.plot(x, results_c[2,:],'-', c='blue', label='$R_{cry}$')
plt.plot(x, results_c[3,:],'--',c= 'red', label='$T_{cry}$')
plt.plot(x, results_a[2,:],'-', c='yellow', label='$R_{amo}$')
plt.plot(x, results_a[3,:],'--',c= 'purple', label='$T_{amo}$')
plt.plot(x, (results_c[-2,:]-results_a[-2,:])/(np.pi*2),'-', c='green', label='$\phi_{R}$' )
plt.plot(x, (results_c[-1,:]-results_a[-1,:])/(np.pi*2),'--', c='c', label='$\psi_{T}$' )
plt.legend(loc=1)

stop = timeit.default_timer()
time= stop-start
materials= zip(stack_layers, stack_thickness)
max_r= np.real(np.max(results_c[2,:]))
max_r_index= np.where(results_c[2,:]==max_r)
plt.title(list(materials))
plt.suptitle(pol +' spectrum Inc: '+ str(theta_inc)+' ${^\circ}$'+
             ','+'Max R: {:.3f}@ {}nm '.format(max_r,np.real(x[max_r_index])[0])+
             ','+ 'Time:{:.3f}s'.format(time))



