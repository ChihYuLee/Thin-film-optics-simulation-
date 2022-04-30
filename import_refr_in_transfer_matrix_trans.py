#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 13:22:02 2022

@author: chihyulee

Translate "import_refr_in.m" file and "transfer_matrix_trans.m" from Carlos group
"""
import numpy as np 
import pandas as pd 
import warnings
from scipy import interpolate 

###### GST refractive index #######


def read_refractive_index_file(materials,min_lambda=300):
    """

    Parameters
    ----------
    materials : Layer of materials are used 
    min_lambda : Minimum wavelength of interest 

    Returns
    -------
    wavelengths: Wavelenths contained in txt file 
    indexes: Different Indexes of n k as a function of wavelengths in amorphous and crystalline states
    (n_crystalline, k_crystalline, n_amorphous, k_amorphous) or (n, k)
    extrapolated: Extrapolated values 

    """
    
    global wavelengths, indexes, extrapolated
    wavelengths=None
    indexes=None
    extrapolated= None
    
    
    if materials== "GST":
        file_data = np.loadtxt("refractive_index_GST_pei.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]
    
    #with max_lambda_interp=1745
    elif materials=="AIST":
        file_data = np.loadtxt("refractive_index_AIST.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:] # 4 columns 
        
        max_lambda_interp=1745
        extrapolated=np.zeros((max_lambda_interp-min_lambda+1,indexes.shape[1]))
        for i in range(extrapolated.shape[1]):
            f=interpolate.interp1d(wavelengths.flatten(), indexes[:,i], fill_value='extrapolate')
            extrapolated[:,i]= f(np.arange(min_lambda, max_lambda_interp+1))
    
    elif materials== "ITO": 
        file_data = np.loadtxt("refractive_index_ITOam.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:] #only 2 columns 
        indexes=np.tile(indexes, (1, 2))
    
    #Compared with Pei
    elif materials== "Pt":
        file_data = np.loadtxt("refractive_index_ITOam.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]  #only 2 columns
        indexes=np.tile(indexes, (1, 2))
        
    #Compared with Pei
    elif materials== "Ag":
        file_data = np.loadtxt("refractive_index_Ag_Ellipsometry_meausred.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]  #only 2 columns
        indexes=np.tile(indexes, (1, 2))
        
    #Longer Spectrum with max_lambda_interp=10000
    elif materials== "Ag2":
        file_data = np.loadtxt("refractive_index_Ag_long.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[:,1:]  #only 2 columns
        indexes=np.tile(indexes, (1, 2))
        
        max_lambda_interp=10000
        extrapolated=np.zeros((max_lambda_interp-min_lambda+1,indexes.shape[1]))
        for i in range(extrapolated.shape[1]):
            f=interpolate.interp1d(wavelengths.flatten(), indexes[:,i], fill_value='extrapolate')
            extrapolated[:,i]= f(np.arange(min_lambda, max_lambda_interp+1))
        
    #Compared with Pei 
    elif materials== "SiO2":
        file_data = np.loadtxt("refractive_index_SiO.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]  #only 2 columns
        indexes=np.tile(indexes, (1, 2))
        
    elif materials== "Ge":
        file_data = np.loadtxt("refractive_index_Ge.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]  #only 2 columns
        indexes=np.tile(indexes, (1, 2))
        
    elif materials== "TiO2":
        file_data = np.loadtxt("refractive_index_TiO2.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]  #only 2 columns
        indexes=np.tile(indexes, (1, 2))
        
    elif materials== "SbSe":
        file_data = np.loadtxt("refractive_index_SbSe_Ellipsometry_meausred_10nm.txt", dtype=np.float32)
        wavelengths= file_data[:,0][:,None]
        indexes=file_data[int(min_lambda-wavelengths[0][0]):,1:]  #4 columns
        
    return wavelengths, indexes, extrapolated
 
def transfer_matrix_trans(theta0, n0, max_lambda, min_lambda, number_layers, n, k, d, pol):
    
    """
    Parameters
    ------
    theta0 : Incident angle in degrees
    n0 : Refraction index of incident medium 
    max_lambda : Maximum lambda to run simulation / Beware: some materials n and k are given up to 800nm
    min_lambda : Minimum lambda to run simulation
    number_layers :  Total number of layers (without counting air)
    n : Refractive index n of crystalline or amorphous states, shape= (wavelengths,layers) #(500,5)
    k : Refractive index k of crystalline or amorphous states = (wavelengths, layers)#(500,5)
    d : Layers of thickness
    pol : Choose either TE or TM polarization (relevant under oblique incidence)


    Returns
    -------
    #array containing X, A, R, T, Phi, Psi
    X: wavelengths in nm 
    A: Apsorption 
    R: Reflection
    T: Transmittance 
    Phi: Phase difference of R in stack of layers 
    Psi: Phase difference of A in stack of layers 
    
    """
    # Parameters for calculation 
    matrix_lambda= max_lambda-min_lambda
    admitanceinvacuum=2.6544*10**(-3)
    #m= np.dstack([np.identity(2)]* (matrix_lambda+1))
    m= np.array([[1,0],[0,1]],dtype=complex)
    y=np.empty((matrix_lambda+1, number_layers),dtype=complex)
    eta= np.empty((matrix_lambda+1, number_layers),dtype=complex)
    theta= np.empty((matrix_lambda+1, number_layers),dtype=complex)
    delta= np.empty((matrix_lambda+1, number_layers-1),dtype=complex)

    
    
    # First layer (PCM)---> theta1, eta1, y1
    theta1= np.arcsin(np.divide(n0*np.sin(theta0),(n[:,0]-1j*k[:,0])))
    y0=n0*admitanceinvacuum 
    y1= (n[:,0]-1j*k[:,0])*admitanceinvacuum 
    if pol=='TE':
        eta0= y0*np.cos(theta0)
        eta1= y1*np.cos(theta1)
    elif pol=='TM':
        eta0= y0/np.cos(theta0)
        eta1= y1/np.cos(theta1)
    else:
        print('Ther is an error')
 
    y[:,0]=y1
    eta[:,0]=eta1
    theta[:,0]=theta1
    M=np.empty((2,2,matrix_lambda+1), dtype=complex)
    for j in range(matrix_lambda+1):
        M[:,:,j]= m
    
    
    
    # Rest of layers 
    lambda_complete= np.arange(min_lambda, max_lambda+1)*10**(-9)
    A=np.empty(matrix_lambda+1)
    R=np.empty(matrix_lambda+1)
    T=np.empty(matrix_lambda+1)
    Phi=np.empty(matrix_lambda+1)
    Psi=np.empty(matrix_lambda+1)
  
        

    
    for i in range(number_layers-1): #different layers
        
        theta[:,i+1]= np.arcsin(np.divide(n[:,i]*np.sin(theta[:,i]), (n[:,i+1]-1j*k[:,i+1])))
        y[:,i+1]=(n[:,i+1]-1j*k[:,i+1])*admitanceinvacuum
        
        if pol=='TE':
            eta[:,i+1]= y[:,i+1]*np.cos(theta[:,i+1])
        elif pol=='TM':
            eta[:,i+1]= y[:,i+1]/np.cos(theta[:,i+1])
        else:
            print('Ther is an error')
        

        delta[:,i]=2*np.pi*d[i]*np.sqrt(n[:,i]**2
                                        -k[:,i]**2
                                        -n0*(np.sin(theta0))**2
                                        -2j*(n[:,i]*k[:,i]))/lambda_complete
        
        s= np.array([[np.cos(delta[:,i]),1j*np.sin(delta[:,i])/eta[:,i]]
                     ,[1j*eta[:,i]*np.sin(delta[:,i]),np.cos(delta[:,i])]], 
                    dtype= complex)


        
        for j in range(s.shape[-1]): #different wavelengths
            M[:,:,j]= M[:,:,j]@s[:,:,j]
            

    #Outputs       
    MM= np.empty((s.shape[0],s.shape[-1]), complex)
    for j in range(s.shape[-1]):  #different wavelengths
        MM[:,j]= M[:,:,j]@ np.array([1,eta[j][-1]], dtype=complex)
        
    B=MM[0,:]
    C=MM[1,:]
    Y=C/B
    r0= np.divide((eta0-Y),(eta0+Y))
    
    R= r0*np.conj(r0)
    T= np.divide(4*eta0*np.real(eta[:,-1]),(eta0*B+C)*np.conj(eta0*B+C))
    A= np.divide(4*eta0*np.real(B*np.conj(C)-eta[:,-1]),(eta0*B+C)*np.conj(eta0*B+C))
    Phi= np.arctan(np.divide(np.imag(eta[:,-1]*(B*np.conj(C)-C*np.conj(B))),
                   np.real((eta[:,-1])**2)*B*np.conj(B)-C*np.conj(C)))
    Psi= np.arctan(np.divide( -np.imag(eta0*B+C),np.real(eta0*B+C)))
    X= lambda_complete
    output= np.array([X, A, R, T, Phi, Psi])
    
    return output
    
    
    
    
       
 
            
  
            
        
        
        
