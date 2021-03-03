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

        
    # drop all rows in df containing nans in columns selected by list cols
    df = df.dropna(subset=list_cols)
    # take the "list od list" including str and float values and access the values of the list
    # to get the last object --> 'max_yieldreduction'

    df['max_yieldreduction'] = [float(i[-1]) for i in df.yieldreduction_list]
    
    # data is a list of strings, we have to convert to list of floats
    for colname in list_cols:
        l=[]
        for i in df[colname]:
            l.append([float(x) for x in i])
        df[colname] = l
        
    df['crop_variety']= rename(df['crop_variety'])   
    df=df.set_index(df['crop_variety'])
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


def get_threshold_slope_salipa(df,savefig=True,calculate_params=True,ts_parameter = 40000):
    '''
    create linear regressions for inserted datzabase values based on the threshold slope model
    
    returns pd.dataframe
    '''
    
#    https://datascience.stackexchange.com/questions/8457/python-library-for-segmented-regression-a-k-a-piecewise-regression
    
    crop_type = str(df.crop_type.values)
    
    ECe=[]
    res=[]
    
    for i in df.index:
    # ger x and y data
        x = df.salinitylevel_list[df.index == i].values[0]
        y = df.yieldreduction_list[df.index == i].values[0] #convert from percetage to fraction
        
        x=np.asarray(x)
        print(x)
#        print(np.asarray(y))
        y=np.asarray(y)*0.01
        print(y)
        crop = df.crop_variety.loc[df.index == i][i]
        crop_type = df.crop_type.loc[df.index == i][i]

        perr_min = np.inf
        p_best = None
        
        if calculate_params:
            for n in range(ts_parameter):
                #crete a tuple of two values between 0 and 1 and multply it by 50 to get random parameterguesses between 0 an 25
                x0_initial = np.random.rand(1)*50
                y0_initial = np.random.rand(1)*2000
                # 
                p , e = curve_fit(piecewise_calibration, x, y,bounds=((0,0), (1000,1000)),p0=(x0_initial,y0_initial),maxfev=ts_parameter)
                perr = np.sum((np.abs(y-piecewise_calibration(x, *p)))**2)
                if(perr < perr_min):
                    perr_min = perr
                    p_best = p
    #                print(p_best,'n: ',n)
            xd = np.linspace(0, max(x)+6, 10000)
            print(str(crop)+str(p_best))
            x0 = p_best[0]
            y0  = p_best[1]
        
        else:
            x0 = df['x0']
            y0 = df['y0']
            
#        y_fao_list = get_piecewise_ec_by_yr(fao_ec_list,x0,x1,y0)

        #  we calculate x at y=0 to find starting point of the third piece of the regression function (stepphun 2005)
        x1=(1+y0*x0)/y0
        
        # append x1 to p_best parametereset
        
        p_best=np.asarray([p_best[0],p_best[1],x1])
        print(p_best)
        
        
        r_square = np.corrcoef(y, piecewise(x, *p_best))[0, 1] ** 2
        rmse = scipy.sqrt(sum((y-piecewise(x, *p_best))**2)/len(piecewise(x, *p_best)))
        
        print(("Threshold Slope  linear r2:%.2f   rmse:%.2f") % (r_square,rmse))
        res.append([crop,crop_type,r_square,rmse,x0,x1,y0])
#        ECe.append(y_fao_list)

        
        if savefig:
            plot=plot_piecewise(p_best,xd,x,y,r_square,rmse,savefig=savefig,Name=crop_type+'_'+str(ts_parameter)+'_'+str(df.loc[i].crop_variety))
    res=pd.DataFrame(res,columns=['crop','crop_type','r²','rmse','x0','x1','y0'])
#    df= pd.DataFrame(ECe,columns=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10','EC_0'])
    
#    df = pd.concat([res,df],axis=1)
    df=df.set_index(['crop'])
    return res    


# plotting

def makeFig(w=10,h=4,fontsize=8.):
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = 'Times'
    plt.rcParams['font.cursive'] = 'Apple Chancery'
    plt.rcParams['axes.labelsize'] = fontsize
    plt.rcParams['xtick.labelsize'] =fontsize
    plt.rcParams['ytick.labelsize'] = fontsize  
    fig = plt.figure(frameon=False)
    fig.set_facecolor('white') 
    #DefaultSize = fig.get_size_inches()
    #fig.set_size_inches( (DefaultSize[0]*.7, DefaultSize[1]*1) )
    fig.set_size_inches(w,h)
    return fig  
    
def plot_piecewise(p_best,xd,x,y,r_square,rmse,savefig=True,Name='Name'):
    '''
    create a plot for threshold slope calculations
    '''
    y_out = piecewise(xd, *p_best)
    xlen = int(round(x.max()) + 6)
    fig = makeFig(w=5,h=5,fontsize=12.)
    ax1 = fig.add_subplot(111)
    ax1.set_xlim(0, xlen)
    ax1.set_ylim(0, 1.05)
    ax1.set_xticks(range(0, xlen))
    ax1.set_xticklabels([str(i) for i in range(0,xlen)])
    ax1.set_yticks(np.arange(0,1.01,0.1))
    ax1.set_yticklabels([str(round(i,2)) for i in np.arange(0,1.01,0.1)])
    
    ax1.plot(x,y,label='Observed', marker='o',color="k",linestyle="none")
    ax1.plot(xd, y_out, 'r-',ls='--', 
             label=("Threshold Slope $r²$:%.2f, $rmse$:%.2f") % (r_square,rmse))
    #    ax1.plot(EC_list,Yr_list,c=color,linestyle=ls)
    ax1.scatter(x,y)
    
    # desingn
    ax1.set_xlabel('Salinity $dS \ m^{-1}$')
    ax1.set_ylabel('Yield reduction [-]')
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_title(Name)
    ax1.legend(ncol=1,fontsize=7,loc='lower left',bbox_to_anchor=[.01, 0.01]).draw_frame(False)
    plt.figure(figsize=(5,5))
    if savefig:
        fig.savefig("/home/konrad/Nextcloud/salzpaper/scripte/threshold_slope/testcase_problems/Fig/"+Name+'.png')
    plt.close('all')
    return df


# problem ids:
#problems = [13,60,65,66,67,72,73,80,85,88,100,110,121]
#problems = [121]
#
df = load_df('/home/konrad/Nextcloud/salzpaper/data/salipa_db2018_KB.csv')

#df = df.loc[df['index'].isin(problems)]


#df_ts= get_threshold_slope_salipa(df)