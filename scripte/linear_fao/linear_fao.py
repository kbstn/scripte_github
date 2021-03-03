#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 08:37:09 2020

@author: konrad
"""
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/konrad/Nextcloud/salzpaper/scripte/original')



from salinity_script_final4 import get_linear_fao

fn='/home/konrad/Nextcloud/salzpaper/data/FAO_EC_0_EC_100.csv'

df,xy = get_linear_fao(fn)

df.to_csv('/home/konrad/Nextcloud/salzpaper/data/FAO_values_linear_interploated.csv')
    