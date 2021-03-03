#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 08:15:53 2020

get a total dataset 
apply selectetd filter patterns
return and save a filtered dataset



@author: konrad
"""
import pandas as pd
#
#filter data

def filter_rsquare(df,threshold=0.75):
    '''
    filter given dataframe by rsquare values, droppiung every entry below 
    threshold AND every entry containging r²=1 (because ist only fitted between two points)    
    '''
    
    return df[df["r²"] >= threshold]

def filter_YR(df,threshold=45):
    '''
    take dataframe created with "load_df" and return an new one ontaining only
    data of experiments  with an maximum yield reduction of at least "threshold" percent
    
    '''
    
    return df[df.max_yieldreduction <= threshold]
#    
    
def filter_N(df,threshold=4):
    '''
    take dataframe created with "load_df" and return an new one ontaining only
    data of experiments  with minimum number of threshold N experimewnts
    '''
    
    return df[df.experiment_numbers >= threshold]

def filter_yield_type(df,exclude=['Other','root length']):
    '''
    take dataframe created with "load_df" and return an new one ontaining only
    data of experiments  do not have the yield type in "exclude"
    '''
    return df[~df['yield_type'].isin(exclude)]    


mod=pd.read_csv('/home/konrad/Nextcloud/salzpaper/data/mod_total.csv',index_col=0)   
 
dfn= filter_N(mod,4)
# filter outr all with max yield reduction higher than 45%
dfy=filter_YR(dfn,threshold=45)

#filter out all wit r² worse than 0.8
#dfr=filter_rsquare(dfy,threshold=0.8)
#filter out all containing non specific yield type
mod_accepted=filter_yield_type(dfy,exclude=['Other','root length'])


ts= pd.read_csv('/home/konrad/Nextcloud/salzpaper/data/ts_40000_total.csv',index_col=0)

    #filter out all experiments less than 4 irrigation steps

dfn= filter_N(ts,4)

# filter outr all with max yield reduction higher than 45%
dfy=filter_YR(dfn,threshold=45)

#filter out all wit r² worse than 0.8
#dfr=filter_rsquare(dfy,threshold=0.8)
#filter out all containing non specific yield type
ts_accepted=filter_yield_type(dfy,exclude=['Other','root length'])

# save filtered datasets to csv files
path = '/home/konrad/Nextcloud/salzpaper/data/'

mod_accepted.to_csv(path+'mod_filter.csv')
ts_accepted.to_csv(path+'ts_filter.csv')
