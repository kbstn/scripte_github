#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 08:06:54 2020

@author: konrad
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
# create barplots

# paths
ts_filter ='/home/konrad/Nextcloud/salzpaper/data/ts_filter.csv'
mod_filter ='/home/konrad/Nextcloud/salzpaper/data/mod_filter.csv'
fao_interpolated = '/home/konrad/Nextcloud/salzpaper/data/FAO_values_linear_interploated.csv'


#dataframes
df_ts = pd.read_csv(ts_filter)
df_ts.columns= ['Crop', 'Crop.1', 'Crop Type', 'r²', 'rmse', 'x0', 'x1', 'y0', 'EC 100',
       'EC 97', 'EC 90', 'EC 75', 'EC 50', 'EC 10', 'EC 0', 'Crop Variety',
       'Reference', 'Experiment Type Original', 'Experiment Type',
       'Experiment Numbers', 'Salinitylevel List', 'Yieldreduction List',
       'Yield', 'Yield_type', 'Max. Yieldreduction']

df_ts=df_ts.replace({'Crop Type': {'wheat': 'Wheat','sorghum': 'Sorghum','potato':'Potato','cucumber':'Cucumber',
                             'corn':'Corn','alfalfa':'Alfalfa','date_palm':'Date Palm','tomato':'Tomato'}})





df_mod = pd.read_csv(mod_filter)
df_mod.columns=['Unnamed: 0', 'Crop Type', 'C50', 'p', 'r²', 'rmse', 'EC 100', 'EC 97',
       'EC 90', 'EC 75', 'EC 50', 'EC 10', 'Crop Variety', 'Reference',
       'Experiment Type Original', 'Experiment Type', 'Experiment Numbers',
       'Salinitylevel List', 'Yieldreduction List',
       'Yield', 'Yield_type', 'Max. Yieldreduction']

df_mod=df_mod.replace({'Crop Type': {'wheat': 'Wheat','sorghum': 'Sorghum','potato':'Potato','cucumber':'Cucumber',
                             'corn':'Corn','alfalfa':'Alfalfa','date_palm':'Date Palm','tomato':'Tomato'}})




df_fao = pd.read_csv(fao_interpolated,index_col='crops')
df_fao.columns=['EC 100', 'EC 97', 'EC 90', 'EC 75', 'EC 50', 'EC 10', 'EC 0']

#df_fao=df_fao.replace({'Crop Type': {'wheat': 'Wheat','sorghum': 'Sorghum','potato':'Potato','cucumber':'Cucumber',
#                             'corn':'Corn','alfalfa':'Alfalfa','date_palm':'Date Palm','tomato':'Tomato'}})

fig, axes = plt.subplots(nrows = 2,sharex=True, sharey=True)
fig.set_size_inches(15,15,forward=True)
#ec_vals= ['EC_90','EC_50']
ec_vals= ['EC 90','EC 50']




def plotting(ec_value, ax,df_ts,df_mod,df_fao):
    # hier dann deine plotting function
    # die plotting function verwendet einfach das ax das wir ihr mitgebegn ha ben
    
    crop_types=['Alfalfa','Corn','Cucumber','Date Palm','Potato','Sorghum','Tomato','Wheat']
    sns.set_style("white")
    sns.set_context("paper")
    
    
    df = pd.DataFrame(index=crop_types)
    
    df[ec_value+'_linremmax']=df_ts.groupby('Crop Type')[ec_value].max()
    df[ec_value+'_linremmin']=df_ts.groupby('Crop Type')[ec_value].min()
    df[ec_value+'_linremmean']=df_ts.groupby('Crop Type')[ec_value].mean()
    df[ec_value+'_modremmax']=df_mod.groupby('Crop Type')[ec_value].max()
    df[ec_value+'_modremmin']=df_mod.groupby('Crop Type')[ec_value].min()
    df[ec_value+'_modremmean']=df_mod.groupby('Crop Type')[ec_value].mean()
    df[str(ec_value)+'_fao']=df_fao[ec_value]
    print(df_fao[ec_value])
#    df[ec_value+'_fao']=df_fao[ec_value]
    y =  np.arange(len(crop_types))


    
    # danger: normal bar is drawn from zero to max value, while using "left"
    # the draw of bar starts at "left value" so we need du substract left from x value
    # to get things right
    ax.plot(df[ec_value+'_fao'],y,marker='o', color='r', ls='',label='Maas (1977)')
    
    ax.plot(df[ec_value+'_modremmean'],y,marker='|', color='darkgreen', ls='',label='',mew=1,markersize=28)
    ax.plot(df[ec_value+'_linremmean'],y,marker='|', color='darkblue', ls='',label='',mew=1,markersize=28)
    ax.barh(y,(df[ec_value+'_linremmax'].values-df[ec_value+'_linremmin'].values),
            left=df[ec_value+'_linremmin'],alpha=0.65,label='Threshold-Slope (| = mean)', color ='#4B6CBB')
    
    ax.barh(y,(df[ec_value+'_modremmax'].values-df[ec_value+'_modremmin'].values),
            left=df[ec_value+'_modremmin'],alpha=0.65, label='Modified-Discount, (| = mean)',color='#42A559')
    
    sns.despine(left=True,bottom=True)
    ax.set_yticks(np.arange(len(crop_types)))
    ax.set_yticklabels(crop_types, fontsize=fontsize)
    
    plt.xticks(fontsize=fontsize)
    ax.set_xlim(0,df.max().max()+1)
    ax.grid(False)
    ax.xaxis.grid(color='lightgrey')
    ax.set_xlabel('Salinity $dS \ m^{-1}$', fontsize=fontsize)
    ax.set_title("Yield potential ranges and values at "+ec_value, fontsize=fontsize+4)
    

for ec_val,ax in zip(ec_vals,axes.flatten()):
    
    plotting(ec_val,ax,df_ts,df_mod,df_fao)
    
    
plt.legend(loc="lower right", fancybox=True)
plt.savefig('/home/konrad/Nextcloud/salzpaper/scripte/ec_barplots/ec_barplots.png', dpi=200)