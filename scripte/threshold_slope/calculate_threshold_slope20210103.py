#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 08:04:04 2020

@author: konrad
"""

import os
import pandas as pd
import scipy.stats
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt


import numpy as np

# calculate threshold slope data for experiments
def load_df(fn, list_cols=['yieldreduction_list','salinitylevel_list']):
    '''
    Load a csv file to a pandas DataFrame. Calculate a new column containing
    values for max yield reduction based on salinity and yield reduction list columns
    also clears and split columnnames to get the names right
    and drop nan values
    
    input:
        fn = filename of csv file
        list_cols= definition of two columns containg 1st lists of 
        yield reduction values and 2nd list of salinity levels (must have same n)
    output:
        pd.DataFrame 
    
    '''
    
    
    # load excelcsv file to a pandas dataframe
    df = pd.read_csv(fn,delimiter=';',encoding='latin1')
    
    
    # covert liststring columns to real python lists
    for colname in list_cols:
        df[colname] = df[colname].str.replace(']',"")
        df[colname] = df[colname].str.replace('[',"")
        df[colname] = df[colname].str.split(',')

    # take the "list od list" including str and float values and access the values of the list
    # to get the last object --> 'max_yieldreduction'

    df['max_yieldreduction'] = [float(i[-1]) for i in df.yieldreduction_list]
    
    # data is a list of strings, we have to convert to list of floats
    for colname in list_cols:
        l=[]
        for i in df[colname]:
            l.append([float(x) for x in i])
        df[colname] = l
                # drop nan alues out of the lists
        df[colname] = df[colname].apply(lambda x: [i for i in x if str(i) != "nan"])

        
    df['crop_variety']= rename(df['crop_variety'])   
#    df=df.set_index(df['crop_variety'])
    return df 

 
def rename(column_name):
    '''
    check column and if a names is duplicated add a number to it
    this way we cann get unique index by having duoplicated names
    
    takes:
        a pandas dataframe column, checks for duplicate strings and if found
    adds a number to it to make it unique
    
    returns: a list containing renamed items in previous given order
    
    '''
    result = []
    for fname in column_name.values:
        orig = fname
        i=1
        while fname in result:
            fname = orig + str(i)
            i += 1
        result.append(fname)
    return result    
    

def piecewise_calibration(x,x0,y0):
    '''
    piecewise function according to Stepphun et al 2005. Root zone salinity
    called Threshold slope function. This function is used for calibrating the
    piecewise regressions. It ignores the third part of the stepphun equation to
    take all datapoints into account when optimizing. 
    
    
    input: 
        x: list of x values
        x0,y0 parameters of function
    
    output:
        np.array containing results of calculation
    '''
#    also applies little 'crack' to avoid sudden drops of the graph because of 
#    the problem if x1 > x where it should not be
#    if x[-1] > x1:
#        x1=x[-1]
        
    pw = np.piecewise(x , [x <= x0, x0<x], [lambda x:1, lambda x: 1-1*y0*(x-x0)])
#    print(pw)
    return pw




def piecewise(x,x0,x1,y0):
    '''
    piecewise function according to Stepphun et al 2005. Root zone salinity
    called Threshold slope function
    
    input: 
        x: list of x values
        x0,x1,y0 parameters of function
    
    output:
        np.array containing results of calculation
    '''
#    also applies little 'crack' to avoid sudden drops of the graph because of 
#    the problem if x1 > x where it should not be
    if x[-1] > x1:
        x1=x[-1]
        
    pw = np.piecewise(x , [x <= x0, np.logical_and(x0<x, x<=x1), x>=x1] , 
                             [lambda x:1, lambda x: 1-1*y0*(x-x0), lambda x: -5])
#    print(pw)
    return pw


def parameter_guess(repetitions,salinity,rel_yield,crop_type,crop):
    '''
    
    '''
    
    # satring variable for optimization process. Storing the minimal estimatetd error. Starting in inf to make it bigger the maximum possible error.
    perr_min = np.inf
    # variable for the best estimated parameterset
    p_best = None
    for iteration in range(repetitions):
            
        #crete a tuple of two values between 0 and 1 and multply it by 50 to get random parameterguesses between 0 an 25
        x0_initial = np.random.rand(1)*50
        # same for 
#            y0_initial = np.random.rand(1)*2000 dont know why this is times 2000 as y is a fraction and should be betwee 0 and 1
        
        y0_initial = np.random.rand(1)
        # fitting a piceswise function on the raw data. Using salinity on x
        # and rel_yield as y data
        # 
        
        # not sure if maxfev is necceary at such high value?
        
        p , e = curve_fit(piecewise_calibration, salinity, rel_yield,
                          bounds=((0,0), (np.inf,np.inf)),p0=(x0_initial[0],y0_initial[0]))#,
#                              maxfev=repetitions)
        
        # sum squared error for differenzce between measured and calibrated values
        perr = np.nanmean((rel_yield - piecewise_calibration(salinity, *p))**2)
#        perr = np.sum((np.abs(rel_yield-piecewise_calibration(salinity, *p)))**2)
        
        # check if error of actual parameter_guess is lower then previous ones
        # if so set error as new minimal error value an store the parameter set 
        # in p_best vairable
        
        if(perr < perr_min):
            perr_min = perr
            p_best = p
        
            print('crop: ',crop_type,' ',crop,' guess: ',iteration,' perr_min: ',perr_min,' y0_initial: ',y0_initial,' x0_initial: ',x0_initial)                
  

    return p_best

    
def get_threshold_slope_salipa(raw_data,repetitions = 40000):
    '''
    create linear regressions for inserted datzabase values based on the threshold slope model
    
    returns pd.dataframe
    '''
    
#    https://datascience.stackexchange.com/questions/8457/python-library-for-segmented-regression-a-k-a-piecewise-regression
    
    
    # empty list to append results later on
    res=[]
    
    
    # iterate every entry in selected raw_data

    
    for entry in raw_data.index:
    # ger x and y data
    
        # x axis data is a list of values for irrigation water salinitiy in ds/m
        salinity = raw_data.salinitylevel_list[raw_data.index == entry].values[0]
        
        # convert list to np.array
        salinity = np.asarray(salinity)
        print(salinity)
        # y axis data is a list of values for corresponding relative yield    
        rel_yield = raw_data.yieldreduction_list[raw_data.index == entry].values[0] #convert from percetage to fraction
        
        # 
        # convert list to np.array and from percentage to fraction 
        rel_yield=np.asarray(rel_yield)*0.01
        print(rel_yield)

        crop = raw_data.crop_variety.loc[raw_data.index == entry][entry]
        crop_type = raw_data.crop_type.loc[raw_data.index == entry][entry]

        # start parameter estimation by repitively guessing random parametersets 
        # and optimizing for minimal error
        
        p_best =parameter_guess(repetitions,salinity,rel_yield,crop_type,crop)              
        x0 = p_best[0]
        print('X0 ',x0)
        y0 = p_best[1]
    

        #  we calculate x at y=0 to find starting point of the third piece of the regression function (stepphun 2005)
        x1=(1-y0*x0)/y0
        


        # append x1 to p_best parametereset
        
        p_best=np.asarray([p_best[0],p_best[1],x1])
        print(p_best)
        
        
        r_square = np.corrcoef(rel_yield, piecewise(salinity, *p_best))[0, 1] ** 2
        rmse = scipy.sqrt(sum((rel_yield-piecewise(salinity, *p_best))**2)/len(piecewise(salinity, *p_best)))
        
        
        name = crop_type+'_'+str(repetitions)+'_'+str(df.loc[entry].crop_variety)
        
        print(("Threshold Slope  linear r2:%.2f   rmse:%.2f") % (r_square,rmse))
        res.append([name,crop,crop_type,r_square,rmse,p_best,x0,x1,y0,salinity,rel_yield])
#        ECe.append(y_fao_list)


     
    res=pd.DataFrame(res,columns=['name','crop','crop_type','r²','rmse','p_best','x0','x1','y0','salinity','rel_yield'])
#    df= pd.DataFrame(ECe,columns=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10','EC_0'])

    return res
# plotting

def plot(x0,x1,salinity,rel_yield,crop,crop_type,r_square,rmse,name):
    
    
    ax = plt.gca()
    # plot 1 first segment
    ax.plot([0,x0],[1,1],linestyle='--',color='red')
    # plot 2 second segment
    ax.plot([x0,x1],[1,0],linestyle='--',color='red')

    # plot 3 third segment
    ax.plot([x1,250],[0,0],linestyle='--',color='red',
            label=("Threshold Slope $r²$:%.2f, $rmse$:%.2f") % (r_square,rmse))

    # plot 4 observation
    ax.scatter(salinity,rel_yield,color='black',label='Observed', marker='o')
    
    ax.set_xlim(0,x1+x0)
    ax.set_ylim(0,1.3)
    ax.set_xlabel('Salinity $dS \ m^{-1}$')
    ax.set_ylabel('Yield reduction [-]')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title(crop)
    ax.legend(ncol=1,fontsize=7,loc='lower left',bbox_to_anchor=[.01, 0.01]).draw_frame(False)
    plt.savefig(name+'.png',dpi=150,bbox_inches='tight')
    
#    plt.close()

    
    
    
        

# problem ids:
#problems = [13,60,65,66,67,72,73,80,85,88,100,110,121]
problems = [133]
#
df = load_df('/home/konrad/Nextcloud/salzpaper/data/salipa2020.csv')

# help in finding errors
# df2=pd.DataFrame(df['salinitylevel_list'].apply(len)-df['yieldreduction_list'].apply(len))



df = df[df.index.isin(problems)]
#

result= get_threshold_slope_salipa(df,repetitions=15000)

result.to_csv('ts_result.csv',sep=';')
#
#    
for entry in result.index:
    
    x0 = result.loc[entry,'x0']
    x1 = result.loc[entry,'x1']
    salinity = result.loc[entry,'salinity']
    rel_yield = result.loc[entry,'rel_yield']
    crop = result.loc[entry,'crop']
    crop_type = result.loc[entry,'crop_type']
    r_square = result.loc[entry,'r²']
    rmse = result.loc[entry,'rmse']
    name = result.loc[entry,'name']
    
    
    
    plot(x0,x1,salinity,rel_yield,crop,crop_type,r_square,rmse,name)
    