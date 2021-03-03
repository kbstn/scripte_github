#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 17:23:22 2021

@author: konrad
"""
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from ast import literal_eval



def ModifiedDiscountFunction(C,C50):
    return (1. / (1. + (C/C50)**3.))



def get_moddisc(raw_data):
    # empty list to append results later on
    res=[]
    # iterate every entry in selected raw_data

    
    for entry in raw_data.index:
    # ger x and y data
    
        # x axis data is a list of values for irrigation water salinitiy in ds/m
        salinity = raw_data.salinity[raw_data.index == entry].values[0]
        
        # convert list to np.array
#        salinity = np.asarray(salinity)
        print(salinity)
        # y axis data is a list of values for corresponding relative yield    
        rel_yield = raw_data.rel_yield[raw_data.index == entry].values[0] #convert from percetage to fraction
        
        # 
        # convert list to np.array and from percentage to fraction 
#        rel_yield=np.asarray(rel_yield)*0.01
        print(rel_yield)

        crop = raw_data.crop.loc[raw_data.index == entry][entry]
        crop_type = raw_data.crop_type.loc[raw_data.index == entry][entry]


        popt, pcov = curve_fit(ModifiedDiscountFunction, salinity, rel_yield)
        res.append(*popt)

    # return list of parameters for every entry in raw_data dataframe
    return res
    
    
def get_moddiscold(df,savefig=True):
    '''

    :return:
        
    '''
    ### read datatable
    
    ### variable for results
    
    res = []
    xy = pd.DataFrame(index=np.linspace(0,40,50))
    
    cropstr= df.crop_variety.values
    
    crop_type = df.crop_type.values
    
    ECe = []
    ### make analysis for all crops in df
    for i in df.index:
        #ger x and y data
        xdata = df.salinitylevel_list[df.index == i].values[0]
        ydata = df.yieldreduction_list[df.index==i].values[0] #convert from percetage to fraction
        
        xdata=np.asarray(xdata)
        ydata=np.asarray(ydata)/100
        xdata = xdata[np.logical_not(np.isnan(xdata))]
        ydata = ydata[np.logical_not(np.isnan(ydata))]
        crop = df.crop_variety.loc[df.index == i][i]
        crop_type = df.crop_type.loc[df.index == i][i]
        
        trialX = np.linspace(0,40,50)
        # fit simple linear function
      #  popt, pcov = curve_fit(SimpleLinearFunction, xdata, ydata)
     #   y_SimpleLinearFunction = SimpleLinearFunction(trialX, *popt)
        # fit ModifiedDiscountFunction
        popt, pcov = curve_fit(ModifiedDiscountFunction, xdata, ydata)
        y_ModifiedDiscountFunction = ModifiedDiscountFunction(trialX, *popt)
        print(y_ModifiedDiscountFunction)
        # define wanted ec values ( same as FAO uses)
        fao_ec_list = np.asarray([1.00,0.97,0.90,0.75,0.5,0.1])
        par_ModifiedDiscountFunction = popt

        # calculate ECx values for definde ec values in fao_ec_list
        c50 = par_ModifiedDiscountFunction[0]
        ECtemp = []
        for i in fao_ec_list:
            ECtemp.append(c50 * (((1/(float(i)))-1)**(1/3)))
            # alöways gives -1????? -1000.0

        
        ECe.append(ECtemp)
        
        
        # xy[crop.decode('UTF-8')+'_x'] = trialX
        #xy[crop.decode('UTF-8') +fn[:-4]+'moddisc'] = y_ModifiedDiscountFunction
        xy=y_ModifiedDiscountFunction    
  

        p = 3#par_ModifiedDiscountFunction[1]
        r_square = np.corrcoef(ydata, ModifiedDiscountFunction(xdata, *popt))[0,1]**2
        rmse = scipy.sqrt(sum((ydata-ModifiedDiscountFunction(xdata, *popt))**2)/len(ModifiedDiscountFunction(xdata, *popt)))
        print(("Modified Discount Function C50:%.2f   p:%.2f   r2:%.2f   rmse:%.2f") % (c50,p,r_square,rmse))


        # fit ModifiedWeilbullFunction
     #   popt, pcov = curve_fit(ModifiedWeilbullFunction, xdata, ydata)
      #  y_ModifiedWeilbullFunction = ModifiedWeilbullFunction(trialX, *popt)
        # fit BiExponetioalResponseFunction
      #  popt, pcov = curve_fit(BiExponetioalResponseFunction, xdata, ydata)
     #   y_BiExponetioalResponseFunction = BiExponetioalResponseFunction(trialX, *popt)
        ## fit ModifiedGompertzFunction --> not working
        #popt, pcov = curve_fit(ModifiedGompertzFunction, xdata, ydata)
        #y_ModifiedGompertzFunction = ModifiedGompertzFunction(trialX, *popt)
        # create figure and set figure properties

        xlen = int(round(xdata.max())+6)

        fig = makeFig(w=5,h=5,fontsize=8.)
        ax1 = fig.add_subplot(111)
        ax1.set_xlim(0,xlen)
        ax1.set_ylim(0,1.01)
        ax1.set_xticks(range(0,xlen))
        ax1.set_xticklabels([str(i) for i in range(0,xlen)])
        ax1.set_yticks(np.arange(0,1.01,0.1))
        ax1.set_yticklabels([str(round(i,2)) for i in np.arange(0,1.01,0.1)])
        ax1.set_xlabel('Salinity $dS \ m^{-1}$')
        ax1.set_ylabel('Yield reduction [-]')
        ax1.get_xaxis().tick_bottom()
        ax1.get_yaxis().tick_left()
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)


        #print(df.crop_variety.loc[df.index == i][i])
       # crop = 'crop'
        
        print(crop)

        ax1.set_title(crop_type+" "+crop)
        # plot observed
        ax1.plot(xdata, ydata, label='Observed', marker='o',color="k",linestyle="none")
        ax1.plot(trialX, y_ModifiedDiscountFunction, 'r-',ls='--', 
                 label=("Modified Discount $r²$:%.2f, $rmse$:%.2f") % (r_square,rmse))
        # plot estimated data
        ax1.legend(ncol=1,fontsize=7,loc='lower left',bbox_to_anchor=[.01, 0.01]).draw_frame(False)
#        bbox_to_anchor=[.95, 1.01]
        plt.figure(figsize=(5,5))
        plt.close('all')

        # get rid of numpybite sting 'b' stuff

        if savefig:
            fig.savefig("Fig/"+"moddisc_"+crop_type+crop+".png",transparent=True,dpi=300)
        res.append([crop_type,c50,p,r_square,rmse])

    import pandas as pd

    res_mod = pd.DataFrame(res,index=cropstr,columns=['crop_type','C50','p','r²','rmse'] )
    EC0_100 = pd.DataFrame(ECe,index=cropstr,columns=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10'])

    res = pd.concat([res_mod,EC0_100],axis=1)
    res.to_csv('res_'+fn[:-4]+'.csv')
    return res,xy,crop




import ast

def from_np_array(array_string):
    array_string = ','.join(array_string.replace('[ ', '[').split())
    return np.array(ast.literal_eval(array_string))



result_ts = pd.read_csv('/home/konrad/Nextcloud/salzpaper/data/ts_result.csv', encoding='latin1',converters={'salinity':from_np_array,'rel_yield':from_np_array})

#list_cols=['salinity','rel_yield']

result_ts['mod_popt'] = get_moddisc(result_ts)