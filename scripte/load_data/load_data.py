#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 17:15:04 2021

Create a handsome csv file out of xlsx database document


@author: konrad
"""

import pandas as pd

filename = '/home/konrad/Nextcloud/salzpaper/data/Database.xlsx'


xlf = pd.ExcelFile(filename)

sheetnames= xlf.sheet_names


salipa = pd.DataFrame(columns=['crop_variety','reference',
                               'experiment_type_original','experiment_numbers',
                               'yield','year','crop_type','salinitylevel_list',
                               'yieldreduction_list','max_yieldreduction'])


for name in sheetnames:
    if name != 'Sheet1':
        print(name)
        tmp_df1 = pd.read_excel(filename,header=1,skiprows=0,sheet_name=name)
        
        tmp_df2 = tmp_df1[['Crop variety','Reference','pot, field or greenhouse',
                         'Number of Tests', 'Variable_Cropyield','year of planting']]
        
        tmp_df2['yieldreduction_list'] = tmp_df1[['% yield_1','% yield_2','% yield_3',
      '% yield_4','% yield_5','% yield_6','% yield_7','% yield_8','% yield_9',
      '% yield_10','% yield_11','% yield_12']].values.tolist()

        tmp_df2['salinitylevel_list'] = tmp_df1[['salt_1 [dS/m]','salt_2 [dS/m]',
                                      'salt_3 [dS/m]','salt_4 [dS/m]',
                                      'salt_5 [dS/m]','salt_6 [dS/m]',
                                      'salt_7 [dS/m]','salt_8 [dS/m]',
                                      'salt_9 [dS/m]','salt_10 [dS/m]',
                                      'salt_11 [dS/m]','salt_12 [dS/m]']].values.tolist()
        tmp_df2['crop_type'] = name
        tmp_df2=tmp_df2.rename(columns={'Crop variety':'crop_variety',
                               'Reference':'reference',
                               'pot, field or greenhouse':'experiment_type_original',
                               'Number of Tests':'experiment_numbers',
                               'Variable_Cropyield':'yield',
                               'year of planting':'year'})
        
        salipa=salipa.append(tmp_df2,ignore_index=True)

salipa.to_csv('salipa2020.csv',sep=';')