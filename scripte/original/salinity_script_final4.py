# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 09:39:27 2018

@author: bestian-k
"""

import pandas as pd
import scipy.stats
from scipy.optimize import curve_fit

from numba import jit

import numpy as np
from matplotlib import pyplot as plt
from IPython.display import display, HTML
import seaborn as sns

# dfine filenames
fn = 'FAO_EC_0_EC_100.csv'
fn2 ='salipa_db2018_KB.csv'
fn3 = 'MaasECvalues.csv'


'''
           fao_EC_100  fao_EC_0
crop_type                      
alfalfa           2.0        16
corn              1.7        10
cucumber          2.5        10
date_palm         4.0        32
potato            1.7        10
sorghum           6.8        13
tomato            2.5        13
wheat             6.0        20

'''

# define YR steps for EC calculations 
fao_ec_list = np.asarray([1.,0.97,0.90,0.75,0.5,0.1,0.])

# parameter defining how often piecwise function get applied to check for 
# parameter with minimum error
# define functions
#to load data:

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
    
def load_FAO(fn):
    '''
    load ec values from fao database
    '''
    return pd.read_csv(fn,delimiter=';',index_col=0)
    
# filter functions 

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
#def filer_min_EC(df,threshold=)
    
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
    

def drop_fliers_from_df(df, col_names=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10']):
    '''
    Take a dataframe containing "crop_type" and "EC_XX" columns. calculate for 
    every crop type and every EC value given in the"col_names list the q25 and q75
    quantlies. Calculate the fences to identify "extreme outlieres". drop 
    identified outliers from origin dataframe and return new dataframe
    '''     
    
    # group data by crop type to enable crop specific calcualtions
    grouped =df.groupby(by='crop_type')
    
    # calculate for every crop type
    for name,group in grouped:
        
        # and for every defined EC value
        for ecvalue in col_names:
            
           # calculate quantilesdf
            q1 = group[ecvalue].quantile(0.25)
            q3 = group[ecvalue].quantile(0.75)
            iqr = q3-q1 #Interquartile range
            # define the fences for extreme outlier identification
            fence_low  = q1-1.5*iqr
            fence_high = q3+1.5*iqr
            
            # get index for every value lower than the low fence 
            toolow = group[group[ecvalue] < fence_low].index
            # drop those indices from ther original dataframe
            for index in toolow:
                if index in df.index:
                    df= df.drop(index)
            # get index for every value higher than high fence
            toohigh = group[group[ecvalue] > fence_high].index
            # drop those indices from the original dataframe
            for index in toohigh:
                if index in df.index:
                    df= df.drop(index)
    # return new dataframe containing only values inbetween the fences               
    return df
# define functions used in the script,


def SimpleLinearFunction(C,a,b):
    return a - b * C

   
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
                             [lambda x:1, lambda x: 1-1*y0*(x-x0), lambda x: 0])
#    print(pw)
    return pw
    
def ModifiedDiscountFunction(C,C50):
    return (1. / (1. + (C/C50)**3.))

def Get_EC_for_Moddisc(Yr,c50,p=3):

    return c50 * (((1/Yr)-1)**(1/p))
    

def get_ec_by_yr(Yr,a,b):
    return ((Yr-a)/-b)
@jit
def get_piecewise_ec_by_yr(Yr,x0,x1,y0):
    

    x = 1/y0+x0-Yr/y0
    
    return x    
    
def create_mean_ec_values(df,ec=['EC_97','EC_10'],db='lin'):
    '''
    Take a dataframe containing 'EC_XX' and 'crop_type' columns and create 
    a mean value for every crop_type, return a new dataframe with crop_type and 
    selected EC values in mean
    also takes db keyword to add it to columnnames to avoid confusions in further work
    '''
    
    mean_df = df.groupby(['crop_type']).mean()
    mean_df= mean_df[ec]
    columns = []
    for i in mean_df.columns:
        columns.append(db+'_'+i)
    mean_df.columns=columns
    return mean_df

    
def get_linear_fao(fn):
    '''

    :return:
    '''
    df = pd.read_csv(fn,index_col=0)
    
    
    
    df['YR_100'] = 100
    df['YR_0'] = 0
    
    
    df = df.T
    # ec 100 + ec 0
    ec=df[:2]
    # relyield 100 + 0
    yr=df[2:4]
    

    crops = df.columns
    print(crops)
    ### variable for results
    ECe = []
    xy = pd.DataFrame(index=np.linspace(1, 40, 50))
    rsqr=[]
    rmselist=[]
    
    
    ### make analysis for all crops
    for crop in df.columns:
        print(crop)
        # ger x and y data
        xdata = ec[crop].values
        ydata = yr[crop].values / 100.  # convert from percetage to fraction
        # create x data for fitting
        trialX = np.linspace(1, 40, 50)
        # fit simple linear function
        coeffs, matcov = curve_fit(SimpleLinearFunction, xdata, ydata )
        a = coeffs[0]
        b = coeffs[1]
        print('a: ',a,' b: ', b)
        y_fao_list = get_ec_by_yr(fao_ec_list, *coeffs)
        
        print(fao_ec_list) # [1.   0.97 0.9  0.75 0.5  0.1  0.  ]
        print(y_fao_list) # [ 2.5    2.725  3.25   4.375  6.25   9.25  10.   ]
        ECe.append(y_fao_list)
        
#        plot=plot_linear(savefig=True,Name=crop,Yr_list=fao_ec_list,EC_list=y_fao_list,x=xdata,y=ydata)
        
        

#        fig.savefig("linear_" + crop + ".png", transparent=True, dpi=300)
    df= pd.DataFrame(ECe,columns=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10','EC_0'])
    df['crops']= crops
    df=df.set_index(['crops'])




    return df,xy
    
    
    
    
## @jit shouild increas calculation time!!! 

def get_threshold_slope_salipa(df,savefig=True,calculate_params=True,ts_parameter = 60000):
    '''
    create linear regressions for inserted datzabase values based on the threshold slope model
    
    returns pd.dataframe
    '''
    
    crop_type = str(df.crop_type.values)
    
    ECe=[]
    res=[]
    
    for i in df.index:
    # ger x and y data
        x = df.salinitylevel_list[df.index == i].values[0]
        y = df.yieldreduction_list[df.index == i].values[0] #convert from percetage to fraction
        
        x=np.asarray(x)
#        print(np.asarray(y))
        y=np.asarray(y)
    
        crop = df.crop_variety.loc[df.index == i][i]
        crop_type = df.crop_type.loc[df.index == i][i]

        perr_min = np.inf
        p_best = None
        
        if calculate_params:
            for n in range(ts_parameter):
                k = np.random.rand(3)*50
                # 
                p , e = curve_fit(piecewise, x, y,bounds=((0,0,0), (np.inf,np.inf,np.inf)),p0=k,maxfev=ts_parameter)
                perr = np.sum(np.abs(y-piecewise(x, *p)))
                if(perr < perr_min):
                    perr_min = perr
                    p_best = p
    #                print(p_best,'n: ',n)
            xd = np.linspace(0, max(x)+6, 10000)
            print(str(crop)+str(p_best))
            x0 = p_best[0]
            x1 = p_best[1]
            y0 = p_best[2]
        
        else:
            x0 = df['x0']
            x1 = df['x1']
            y0 = df['y0']
        y_fao_list = get_piecewise_ec_by_yr(fao_ec_list,x0,x1,y0)
        r_square = np.corrcoef(y, piecewise(x, *p_best))[0, 1] ** 2
        rmse = scipy.sqrt(sum((y-piecewise(x, *p_best))**2)/len(piecewise(x, *p_best)))
        print(("Threshold Slope  linear r2:%.2f   rmse:%.2f") % (r_square,rmse))
        res.append([crop,crop_type,r_square,rmse,x0,x1,y0])
        ECe.append(y_fao_list)

        
        if savefig:
            plot=plot_piecewise(p_best,xd,x,y,r_square,rmse,savefig=savefig,Name=crop_type+'_'+str(ts_parameter)+'_'+str(df.ix[i].crop_variety))
    res=pd.DataFrame(res,columns=['crop','crop_type','r²','rmse','x0','x1','y0'])
    df= pd.DataFrame(ECe,columns=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10','EC_0'])
    
    df = pd.concat([res,df],axis=1)
    df=df.set_index(['crop'])
    return df    
    
def get_moddisc(df,savefig=True):
    '''

    :return:
    '''
    ### read datatable
    
    ### variable for results
    
    from matplotlib.ticker import FormatStrFormatter
    
    
    res = []
    import pandas as pd
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
        fig.savefig("Fig/"+Name+'.png')
    plt.close('all')
    return df

def make_barplot(lindfrem,moddfrem,faodf,ec_value='EC_90',crop_types=['wheat', 'potato',
                                                                'cucumber', 'corn',
                                                                'alfalfa', 'date_palm',
                                                                'tomato', 'sorghum'],
                                                                savefig=True, fontsize=12):
    '''
    Create a multi horizontal barplot:
        - for every crop type a row
        - each row contains 2 (overlapping) bars
        - each bar made of max and min for selected ec value
        - additionally plot a point for FAO values if avaliable
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    
    sns.set_style("white")
    sns.set_context("paper")
    
    
    df = pd.DataFrame(index=crop_types)

    df[ec_value+'_linremmax']=lindfrem.groupby('crop_type')[ec_value].max()
    df[ec_value+'_linremmin']=lindfrem.groupby('crop_type')[ec_value].min()
    df[ec_value+'_linremmean']=lindfrem.groupby('crop_type')[ec_value].mean()
    df[ec_value+'_modremmax']=moddfrem.groupby('crop_type')[ec_value].max()
    df[ec_value+'_modremmin']=moddfrem.groupby('crop_type')[ec_value].min()
    df[ec_value+'_modremmean']=moddfrem.groupby('crop_type')[ec_value].mean()
    df[ec_value+'_fao']=faodf.groupby(df.index)[ec_value].max()
    
    df = df.ix[crop_types]

    y =  np.arange(len(crop_types))
    
    fig, ax = plt.subplots()
    fig.set_size_inches(16,9,forward=True)

    # danger: normal bar is drawn from zero to max value, while using "left"
    # the draw of bar starts at "left value" so we need du substract left from x value
    # to get things right
    ax.plot(df[ec_value+'_fao'],y,marker='o', color='r', ls='',label='Maas (1977)')
    
    ax.plot(df[ec_value+'_modremmean'],y,marker='|', color='darkgreen', ls='',label='',mew=1,markersize=22)
    ax.plot(df[ec_value+'_linremmean'],y,marker='|', color='darkblue', ls='',label='',mew=1,markersize=22)
    ax.barh(y,(df[ec_value+'_linremmax'].values-df[ec_value+'_linremmin'].values),
            left=df[ec_value+'_linremmin'],alpha=0.75,label='Threshold-Slope (| = mean)')
    ax.barh(y,(df[ec_value+'_modremmax'].values-df[ec_value+'_modremmin'].values),
            left=df[ec_value+'_modremmin'],alpha=0.65, label='Modified-Discount, (| = mean)')
    print((df[ec_value+'_linremmax']-df[ec_value+'_linremmin']))
    print(df[ec_value+'_linremmin'])
    sns.despine(left=True,bottom=True)
    ax.set_yticks(np.arange(len(crop_types)))
    ax.set_yticklabels(crop_types, fontsize=fontsize)
    
    plt.xticks(fontsize=fontsize)
    ax.set_xlim(0,df.max().max()+1)
    ax.grid(False)
    ax.xaxis.grid(color='lightgrey')
    ax.set_xlabel('Salinity $dS \ m^{-1}$', fontsize=fontsize)
    ax.set_title("Yield potential ranges and values at "+ec_value, fontsize=fontsize+4)
    
#    # create custom legend
#    from matplotlib.patches import Patch
#    from matplotlib.lines import Line2D
#    # seaborn green= 55A868 blie = 4C72B0
#    maas =Line2D(color='red',marker='o', label='Maas (1977)')
#    ts = Patch(color='#4C72B0', label='Threshold-Slope')
#    mod = Patch(color='#55A868', label='Modified-Discount')
#    ts_mean = Line2D(color='#4C72B0',marker='|', label='Threshold-Slope Mean')
#    mod_mean = Line2D(color='55A868',marker='|', label='Modified-Discount')
#
#    legend = ax.legend(handles=[maas, ts, mod, ts_mean, mod_mean],bbox_to_anchor=(0.1, 5.75), frameon=True, fancybox=True, facecolor="white", edgecolor="grey")
##    for text in legend.get_texts():
##        text.set_color("grey")
    plt.legend(loc="lower right", fancybox=True,fontsize=fontsize-2)
    plt.show()
    if savefig:
            fig.savefig("Fig/final_figures/"+"barplots_"+ec_value+".png",transparent=True,dpi=300)
    df.to_csv('Fig/'+ec_value+'_final_values.csv')
    
def pieplot(df,kind='pie',my_color='darkgrey'):
    
    '''
    create pyplot of number of experiments
    '''
    import seaborn as sns
    
    sns.set_style("white")
    
    figsize=(4,4)
    df.experiment_numbers.value_counts().sort_values(ascending=True).plot(kind=kind, legend=False, figsize=figsize, color=my_color)
    # get axis
    ax = plt.gca()
    # get rid of spines
    sns.despine(ax=ax)
    #set label for x axis    
    ax.set_xlabel('Number of experiments')
    ax.set_ylabel('Irrigation steps per experiment')
    plt.show()
    
    df.experiment_type.value_counts().sort_values(ascending=True).plot(kind=kind, legend=False, figsize=figsize, color=my_color)
    ax = plt.gca()
    sns.despine(ax=ax)
    ax.set_xlabel('Number of experiments')
    ax.set_ylabel('Experiment types')
    plt.show()
    
    df.yield_type.value_counts().sort_values(ascending=True).plot(kind=kind, legend=False, figsize=figsize,color=my_color)
    ax = plt.gca()
    ax.set_xlabel('Number of experiments')
    ax.set_ylabel('Irrigation steps per experiment')
    sns.despine(ax=ax)
    plt.show()
    
    df['bins'] = pd.cut(df['max_yieldreduction'],bins=[0,20,40,60,80,150], labels=['0','20','40','60','80'])
    a=df.groupby('bins').size()
    a.plot(kind=kind,figsize=figsize, color = my_color)
    ax = plt.gca()
    ax.set_xlabel('Number of experiments')
    ax.set_ylabel('Maximum observed Yr decline per experiment')
    sns.despine(ax=ax)
    plt.show()

    

    
    
def make_boxplot(lindf, lindftotal,moddf,moddftotal):

    from matplotlib import pyplot as plt
    import pandas as pd
    import seaborn as sns
    sns.set(color_codes=True)
    
    ticks =['threshold slope','modified discount']

    # plotlayout 2x2 plot row and column sharing
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
    # r2 salipa min mod vs lin
    data1 = pd.concat([lindf['r²'].dropna(),moddf['r²'].dropna()],axis=1, names=['threshold slope','modified discount'])
    ax1 = data1.boxplot(ax=ax1)
    ax1.set_title('$Accepted Data, n={}$'.format(len(lindf)))
    ax1.set_ylabel("$R^2$")
    ax1.grid(False)
    #ax1.set_axis_bgcolor("whitesmoke")

    # r2 salipa maax mod vs lin
    data2 = pd.concat([lindftotal['r²'].dropna(),moddftotal['r²'].dropna()],axis=1, names=['threshold slope','modified discount'])
    ax2 = data2.boxplot(ax=ax2)

    ax2.set_title('$Total Data, n={}$'.format(len(lindftotal)))
    # ax2.set_ylabel("$R^2$")
    ax2.grid(False)
    #ax2.set_axis_bgcolor("whitesmoke")


    # rmse salipa min mod vs lin
    data3 = pd.concat([lindf['rmse'].dropna(),moddf['rmse'].dropna()],axis=1, names=['threshold slope','modified discount'])
    ax3 = data3.boxplot(ax=ax3)
    
    # ax3.set_title('$SaliPa_min$')
    ax3.set_ylabel("$rmse$")
    ax3.grid(False)
    ax3.set_xticklabels(ticks,minor=False)
    ax3.set_ylim(-0.01,0.3)
    #ax3.set_axis_bgcolor("whitesmoke")
    # rmse salipa min mod vs lin
    data4 = pd.concat([lindftotal['rmse'].dropna(),moddftotal['rmse'].dropna()],axis=1, names=['threshold slope','modified discount'])
    ax4 = data4.boxplot(ax=ax4)
    ax4.set_xticklabels(ticks,minor=False)
    ax4.set_ylim(-0.01,0.3)
    # ax4.set_title('$SaliPa_max$')
    # ax4.set_ylabel("rmse")
    ax4.grid(False)
    #ax4.set_axis_bgcolor("whitesmoke")

    
    
    plt.show()
# make tables
def make_table(df,savename='table'):
    '''
    
    '''
    
    df2=pd.DataFrame()
    # EC columns
    ec_cols= df.columns[df.columns.str.startswith('EC')]
        
    
    for i in ec_cols:                    
    # get maximum ec values for every crop type
   # maxcol= [i+'_max' for i in ec_cols]
        df2[i+'_max']=df.groupby(by='crop_type').max()[i]
    # get minimun ec value for every crop type
    #mincol= [i+'_min' for i in ec_cols]

        df2[i+'_min']=df.groupby(by='crop_type').min()[i]
    
        #save as excel
        
    writer = pd.ExcelWriter(savename+'.xlsx', engine='xlsxwriter')
    df2.to_excel(writer, sheet_name='Sheet1')

    return df2


def reload_dataframes(fn2):
    '''
    
    '''
    df=load_df(fn2)
    df_total=df
    df = filter_YR(df,threshold=45)
    df= filter_N(df,threshold=4)
    df = filter_yield_type(df,exclude=['Other','root length'])
    
    lindftotal = pd.read_csv('lindftotal.csv')
    moddftotal = pd.read_csv('moddftotal.csv')

    return df,df_total,lindftotal,moddftotal
    



def recalculate_threshold_slope_salipa(df,lindftotal,savefig=True):
    '''
    create linear regressions for inserted datzabase values based on the threshold slope model
    '''
    
    
    ECe=[]
    res=[]
    
    # viterate the filtered df
    for i in df.index:
    # ger x and y data
        x = df.salinitylevel_list[df.index == i].values[0]
        y = df.yieldreduction_list[df.index == i].values[0] #convert from percetage to fraction
        
        x=np.asarray(x)
        y=np.asarray(y)/100
    
        xd = np.linspace(0, max(x)+6, 10000)

        crop = df.crop_variety.loc[df.index == i][i]
        crop_type = df.crop_type.loc[df.index == i][i]
        print(crop)
        r2 = lindftotal['r²'].loc[lindftotal.index == i][i]
        rmse = lindftotal['rmse'].loc[lindftotal.index == i][i]
        x0 = lindftotal['x0'].loc[lindftotal.index == i][i]
        x1 = lindftotal['x1'].loc[lindftotal.index == i][i]
        y0 = lindftotal['y0'].loc[lindftotal.index == i][i]
        p_best = (x0,x1,y0)
  
        print(x0,x1,y0)
        y_fao_list = get_piecewise_ec_by_yr(fao_ec_list,x0,x1,y0)
        r_square = np.corrcoef(y, piecewise(x, x0,x1,y0))[0, 1] ** 2
        rmse = scipy.sqrt(sum((y-piecewise(x, x0,x1,y0))**2)/len(piecewise(x, *p_best)))
        print(("Threshold Slope  linear r2:%.2f   rmse:%.2f") % (r_square,rmse))
        res.append([crop,crop_type,r_square,rmse,x0,x1,y0])
        ECe.append(y_fao_list)

        
        if savefig:
            plot=plot_piecewise(p_best,xd,x,y,r_square,rmse,savefig=savefig,Name=crop_type+'_'+crop)
    res=pd.DataFrame(res,columns=['crop','crop_type','r²','rmse','x0','x1','y0'])
    df= pd.DataFrame(ECe,columns=['EC_100','EC_97','EC_90','EC_75','EC_50','EC_10','EC_0'])
    
    df = pd.concat([res,df],axis=1)
    df=df.set_index(['crop'])
    return df    
    
def make_boxplot_ec(ts,ts_accepted,mod,mod_accepted,EC1='EC_50',EC2='EC_90',savefig=True,font_size=20):
    
    from matplotlib import pyplot as plt
    import pandas as pd
    import seaborn as sns
    sns.set(color_codes=True)
    
    
    for i in pd.unique(ts.crop_type):
        ticks =['threshold slope','modified discount']
    # plotlayout 2x2 plot row and column sharing
        f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        f.set_size_inches(8,9,forward=True)
        f.suptitle(i, fontsize=font_size+4)
    # r2 salipa min mod vs lin
        #first subplot accepted ec50
        data1=pd.concat([mod_accepted[mod_accepted.crop_type==i][EC1],
                     ts_accepted[ts_accepted.crop_type==i][EC1]],axis=1,names=ticks)
        
        ax1 = data1.boxplot(ax=ax1,showfliers=False)
        ax1.set_title('$Accepted Data, n={}$'.format(len(data1)),fontsize=font_size)
        ax1.set_ylabel(EC1,fontsize=font_size)
        ax1.grid(False)
        ax1.yaxis.grid(color='lightgrey')
        ax1.set_facecolor('white')
        ax1.tick_params(axis='y',labelsize=font_size-4)

#        ax1.set_ylim(0.05,)
#        ax1.margins(y=0.05)
        
                #third subplot total ec50 
        data2=pd.concat([mod[mod.crop_type==i][EC1],
                     ts[ts.crop_type==i][EC1]],axis=1,names=ticks)

        ax2 = data2.boxplot(ax=ax2,showfliers=False)
#        ax2.margins(y=0.05)
        # ax3.set_title('$SaliPa_min$')
        ax2.grid(False)
        ax2.set_ylim(0.05,)
        ax2.yaxis.grid(color='lightgrey')
        ax2.set_facecolor('white')
#        ax2.set_xticklabels(ticks,minor=False)
        ax2.set_title('$Total Data, n={}$'.format(len(data2)),fontsize=font_size)
         #second subplot accepted ec90
        data3=pd.concat([mod_accepted[mod_accepted.crop_type==i][EC2],
                     ts_accepted[ts_accepted.crop_type==i][EC2]],axis=1, names=ticks)
        ax3 = data3.boxplot(ax=ax3,showfliers=False)
#        ax3.margins(y=0.05)
        ax3.set_xticklabels(ticks,minor=False,fontsize=font_size,rotation=45)
        ax3.set_ylabel(EC2,fontsize=font_size)
        ax3.grid(False)
        ax3.yaxis.grid(color='lightgrey')
        ax3.set_facecolor('white')
        ax3.tick_params(axis='y',labelsize=font_size-4)
#       ax3.set_ylim(0.05,)
#       fourth total ec90

        data4=pd.concat([mod[mod.crop_type==i][EC2],
                     ts[ts.crop_type==i][EC2]],axis=1,names=ticks)
        ax4 = data4.boxplot(ax=ax4,showfliers=False)
        ax4.set_xticklabels(ticks,minor=False,fontsize=font_size,rotation=45)
        ax4.set_ylim(0.05,)
        ax4.grid(False)
        ax4.yaxis.grid(color='lightgrey')
        ax4.set_facecolor('white')
# ax4.set_title('$SaliPa_max$')
# ax4.set_ylabel("rmse")
        
        # make room on bottom so labels arnt cut off
        plt.gcf().subplots_adjust(bottom=0.25)
        if savefig:
            f.savefig("Fig/final_figures/"+i+'_'+EC1+'_'+EC2+'.png')

        plt.show()

 
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
    

if __name__ == '__main__':
        
    #    
    #def test_significance(column1,column2):
    #    return 
    #    ts['r²]
    
    
    
    df= load_df(fn2)
#    
#    
#    #mod,xy,crop = get_moddisc(df)
#    #df2=load_df(fn2)
#    #columns=['crop_variety', 'reference','experiment_type_original', 'experiment_type', 'experiment_numbers','salinitylevel_list', 'yieldreduction_list', 'yield', 'yield_type',
#    #       'max_yieldreduction']
#    # mod[columns]=df2[columns]
#    
#    
#    mod=pd.read_csv('mod_total.csv',index_col=0)
#    
#    
#    
#    dfn= filter_N(mod,4)
#    # filter outr all with max yield reduction higher than 45%
#    dfy=filter_YR(dfn,threshold=45)
#    
#    #filter out all wit r² worse than 0.8
#    #dfr=filter_rsquare(dfy,threshold=0.8)
#    #filter out all containing non specific yield type
#    dfyt=filter_yield_type(dfy,exclude=['Other','root length'])
#    mod_accepted=dfyt
#    #numer of references cited in accepted df:
#    print(len(mod_accepted.groupby(mod_accepted.reference)))
#    
#    # ts parameter 10000 seems good enough
#    ##ts= get_threshold_slope_salipa(df,savefig=True,calculate_params=True,ts_parameter=40000)
#    #df2=load_df(fn2)
#    #columns=['crop_variety', 'reference','experiment_type_original', 'experiment_type', 'experiment_numbers','salinitylevel_list', 'yieldreduction_list', 'yield', 'yield_type',
#    #       'max_yieldreduction']
#    #df=pd.read_csv('ts_60000.csv')
#    #df=df.set_index(df.crop) 
#    #df[columns]=df2[columns]
#    
#    
#    
#    
#    
#    ts= pd.read_csv('ts_40000_total.csv',index_col=0)
#    # filter procedure:
#        
#        #filter out all experiments less than 4 irrigation steps
#    
#    dfn= filter_N(ts,4)
#    
#    # filter outr all with max yield reduction higher than 45%
#    dfy=filter_YR(dfn,threshold=45)
#    
#    #filter out all wit r² worse than 0.8
#    #dfr=filter_rsquare(dfy,threshold=0.8)
#    #filter out all containing non specific yield type
#    dfyt=filter_yield_type(dfy,exclude=['Other','root length'])
#    ts_accepted=dfyt
#    ##numer of references cited in accepted df:
#    #print(len(ts_accepted.groupby(mod_accepted.reference)))
#    ##
#    #
#    ##daten stimmen nicht!!!! müssen überprüft werden!!! sind nämlich beide dieselben?!!!!!
#    #print(scipy.stats.ttest_ind(mod_accepted['r²'].dropna(),ts_accepted['r²'].dropna(),equal_var=False))
#    ##Out[81]: Ttest_indResult(statistic=0.0, pvalue=1.0)
#    #
#    #print(scipy.stats.ttest_ind(mod['r²'].dropna(),ts['r²'].dropna(),equal_var=False))
#    ##Out[93]: Ttest_indResult(statistic=-3.5772907988389986, pvalue=0.000385363368351187)
#    #
#    #
#    #print(scipy.stats.ttest_ind(mod_accepted['r²'].dropna(),mod['r²'].dropna(),equal_var=False))
#    #print(scipy.stats.ttest_ind(ts_accepted['r²'].dropna(),ts['r²'].dropna(),equal_var=False))
#    
#    #drop fliers
#    
##    
##    ## works
##    ts_accepted = drop_fliers_from_df(ts_accepted,col_names=['EC_50','EC_90'])
##    mod_accepted = drop_fliers_from_df(mod_accepted,col_names=['EC_50','EC_90'])
##    make_boxplot_ec(ts,ts_accepted,mod,mod_accepted)
##    #
##    #
##    
##    
##    
##    get_threshold_slope_salipa(ts,savefig=True,calculate_params=False,ts_parameter = 60000)
##    get_moddisc(df,savefig=True) 
##        
#        
#    #
#    #
#    #ts = drop_fliers_from_df(ts,col_names=['EC_50','EC_90'])
#    #ts_accepted = drop_fliers_from_df(ts_accepted,col_names=['EC_50','EC_90'])
#    #mod = drop_fliers_from_df(mod,col_names=['EC_50','EC_90'])
#    #mod_accepted = drop_fliers_from_df(mod_accepted,col_names=['EC_50','EC_90'])
#    #
#    #
#    ##table_maas
#    fao,xy=get_linear_fao(fn)
#    #writer = pd.ExcelWriter('fao.xlsx', engine='xlsxwriter')
#    #fao.to_excel(writer, sheet_name='Sheet1')
#    #
#    #table_mod= make_table(mod_accepted,savename='table_mod')
#    #table_ts = make_table(ts_accepted,savename='table_ts')
#    #
#    #
#    ## create barplots
#    for i in fao.columns[:-1]:
#        make_barplot(ts_accepted,mod_accepted,fao,ec_value=i,fontsize=20)
    #    
    #    
    #    