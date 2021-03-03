#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 18:45:03 2020

@author: konrad
"""

# get linear fao values

import os
import pandas as pd
import scipy.stats
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import numpy as np

def SimpleLinearFunction(x,a,b):
    return a * x + b

def get_linear_fao(EC_100,EC_50):
    '''

    :return:
        EC 90 wert herausgeben
        
    '''
    
    # Create corresponding yield reduction values for EC_100 nd EC_50
    
    ec_values = (EC_100,EC_50)
    yield_reduction = (1.0,0.5)
    
    # fit linear function to EC_100 and EC_50 values
    
    coeffs, matcov = curve_fit(SimpleLinearFunction, ec_values, yield_reduction)

    a = coeffs[0]
    b = coeffs[1]

    # value for target rel. yield    
    Yr = 0.9
    
    # use fittetd parameters to calculate EC90
    x =  ((Yr -b ) / a)

    return x

EC_100_EC_50_vals_MINHAS = pd.read_csv('/home/konrad/Nextcloud/salzpaper/data/EC100_EC50_vals_MINHAS.csv',index_col=0)
EC_100_EC_50_vals_MINHAS['EC_90'] = None

# calulate EC90 value for each crop

for crop in EC_100_EC_50_vals_MINHAS.index:
    EC_100=EC_100_EC_50_vals_MINHAS.loc[crop,'EC_100']
    EC_50=EC_100_EC_50_vals_MINHAS.loc[crop,'EC_50']
    
    EC90 = get_linear_fao(EC_100,EC_50)
    
    EC_100_EC_50_vals_MINHAS.loc[crop,'EC_90'] = EC90 

EC_100_EC_50_vals_MINHAS.to_csv('/home/konrad/Nextcloud/salzpaper/data/EC_100_EC_90_EC_50_vals_MINHAS.csv')