#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 13:22:02 2022

@author: chihyulee

Translate "import_refr_in.m" file from Carlos group
"""
import numpy as np 
import pandas as pd 
import warnings

###### GST refractive index #######


def read_refractive_index_file(materials):
    """

    Parameters
    ----------
    materials : Layer of materials are used 

    Returns
    -------
    Different Indexes of n k as a function of wavelengths in amorphous and crystalline states 

    """
    if materials== "GST":
        file_data = np.loadtxt("refractive_index_GST_pei.txt", dtype=float)
        print(file_data.shape)
        
read_refractive_index_file('GST')