#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 09:56:35 2021

@author: konrad
"""
import pandas as pd
from matplotlib import pyplot as plt
import ast
from scipy.optimize import curve_fit
import numpy as np
from sklearn.linear_model import LinearRegression





# objective function
def objective(x, a, b, c):
	return a * x + b

def threshold_slope(x,a,b):
    '''
    orginal formula from vangenuchten:
        
    Yr= 1 - b*(C-Ct) at Ct < C < C0
    
    Yr= relative yield
    b = absolute value of declining slope in Yr with C
    Ct = maximimum value of salinity without a yield reduction
    C0 = lowest value of C where Yr = 0
    

    Parameters
    ----------
    x : irrigation water salinity (C).
    a : slope of declining yield (b)
    b : breakpoint point of highest salinty where yield is still 100%. (Ct)

    Returns
    -------
    relative yield in dependency of irrigation water salinty.



    '''
    return 1 -a*(x-b)



# fit curve

def simple_fit(function,xvals,yvals):
    '''
    given linear function and x and y valus this function returns the otimal
    parameters (popt) and ther estimated covariance (pcov) on the base of
    scipy.optimize.curve_fit
    
    Parameters
    ----------
    function : objectiv function. python function expected
    xvals : array like series of x values
    yvals : array like series of y values

    Returns
    -------
    popt,pcov

    '''
    
    return curve_fit(threshold_slope,xvals,yvals)




def get_simple_fit_values(function,x_array, popt):
    '''
    return array like series of y values based on a function and a series of 
    x values

    Parameters
    ----------
    function :  objectiv function. python function expected
    
    x_array: gereated array of data for x
    
    popt :  the optimal
    parameters of a curve_fit function (popt).

    Returns
    -------
    None.

    '''
    a,b = popt
    
    return function(x_array,a,b)


def linearfit_uncentered(x_values,y_values):
    '''
    account for the problem that the breakpoint b is zero or negative. 
    In such cases we have to adjust the calibration to ensure that every fitted 
    regression starts at least at (1|0) following this stackoverflow post:
    https://stackoverflow.com/questions/58531601/linear-fit-constrained-to-go-through-the-first-point-in-python
    
    
    if first y value is less than 1 we dont aplly piecewise but linear fit centered on (1|0)'''
    
    x_values = np.array(x_values)
    
    y_values = np.array(y_values)
    
    
    lm = LinearRegression(fit_intercept = False)
    
    # center data on x_values[0], y_values[0]
    y_values2 = y_values - y_values[0]
    x_values2 = x_values - x_values[0]
    
    # fit model
    lm.fit(x_values2.reshape(-1, 1), y_values2)
    # predict line
    preds = lm.predict(np.arange(0, 40, 0.1).reshape(-1,1))
    

    # add x_0 and y_0 back to the predictions
    new_x_values= np.arange(0,40, 0.1) + x_values[0]
    new_y_values= preds  + y_values[0]
    
    
    
    return new_x_values[new_y_values>0],new_y_values[new_y_values>0]
    


def piecewise(x,a,b):
    '''
    

    Parameters
    ----------
    x : irrigation water salinity (C).
    
    a : slope of declining yield (b)
    
    b : breakpoint point of highest salinty where yield is still 100%. (Ct)

    Returns
    -------
    series of y values (values for relative yield) based on given x values

    '''

#    print(x,a,b)
#    if b < 0:
#        print(type(x),'to low!')
        
    # return np.piecewise(x , [x <= b, b<x], [lambda x:1, lambda x: 1-1*a*(x-b)])
    return np.piecewise(x , [x <= b, b<x], [lambda x:1, lambda x: 1-1*a*(x-b)])

def estimate_breakpoint(function,x_values,y_values,repetitions=5000):
    '''
    

    Parameters
    ----------
    function : TYPE
        DESCRIPTION.
    x_values : TYPE
        DESCRIPTION.
    y_values : TYPE
        DESCRIPTION.
    repetitions : Number of repitition for breakpoint estimation. The default is 40000.

    Returns
    -------
    p_best : TYPE
        DESCRIPTION.
    sqerror : TYPE
        DESCRIPTION.

    '''

    
    # satring variable for optimization process. Storing the minimal 
    # estimatetd error. Starting in inf to make it bigger the maximum 
    # error.
#    sqerror_min = np.inf
#    # variable for the best estimated parameterset
#    p_best = None
#    for iteration in range(repetitions):
#            
    #
#    a_initial = (np.random.rand(1)*10)[0]
#    
#    #initialk parametergues for breakpoint
#    b_initial = (np.random.rand(1)*25)[0]
#    
    
    # fit curve piecwise using initial gusses
    popt,pcov = curve_fit(piecewise,x_values,y_values)#,
                     #      p0=(a_initial,b_initial))
    
    # claculate the sum squared error between calculated and measured data
    sqerror = np.nanmean((y_values - piecewise(x_values, *popt))**2)
    
    # check if guessed error is below the best rum. if so replace best run 
    # with actual parameters
#        if(sqerror < sqerror_min):
#            sqerror_min = sqerror
#            p_best = popt
#        
#            print(' guess: ',iteration,
#                   ' perr_min: ',sqerror_min,' a_initial: ',a_initial,
#                   ' b_initial: ',b_initial)                
#            print('pbest_tmp: ',p_best)

        
    return popt,sqerror
        
def calc_xy_values(x_values,y_values,x_out=np.linspace(0,15,num=100),
                      function=threshold_slope,method=simple_fit,linearfit_uncentered=False):
    '''
    take two lists of x and y values. apply the simple_fit using the given 
    objective function (threshold_slope)
    also generate series of new x values for claculation (x_out, can be changed)

    Parameters
    ----------
    x_values : list of irrigation water salinity values.
    
    y_values : list of relative yield values.
        DESCRIPTION.
    x_out : np array of values. The default is np.linspace(0,15,num=250).
    
    function : objective function for the fitting can be changed here.
    The default is threshold_slope.
    
    method: method used for fitting(simple fit or estimate_breakpoint)
    

    Returns
    -------
    two arrays of same length, first x_out and second for every value in x_out 
    the corresponding calculated value for y_out

    '''
    popt_simple,pcov = method(function, x_values, y_values)

    
        
        
    return x_out,get_simple_fit_values(function,x_out,popt_simple)




def fit_and_plot(x_values,y_values,crop_variety):
    
       #fit for simple linear model
    x_simple,y_simple= calc_xy_values(x_values,y_values,function=threshold_slope,
                                      method=simple_fit)
    # print(x_simple,y_simple)
    
    #fit piecewise
    x_piecewise,y_piecewise = calc_xy_values(x_values,y_values,function= piecewise, 
                                             method=simple_fit)
    # print(x_piecewise,y_piecewise)
    
    #fit breakpoint
    x_bp,y_bp = calc_xy_values(x_values,y_values,function= piecewise,
                               method=estimate_breakpoint)
    
    if y_bp[0] < 1:
        x_bp,y_bp=linearfit_uncentered(x_values,y_values)
        
    # plot the data
    ax = plt.gca()
    ax.scatter(x_values,y_values,color='r')
#    ax.plot(x_simple,y_simple,color='b')
#    ax.plot(x_piecewise,y_piecewise,color='k')
    ax.plot(x_bp,y_bp,color='g')
    plt.title(crop_variety)
    plt.show()
    plt.close()

def xy_crop_from_data(fn='problems.csv'):
    
    data = pd.read_csv(fn,sep=',')
    
    for i in data.index:
        x_values= ast.literal_eval(data.salinitylevel_list.loc[i])
        y_values= ast.literal_eval(data.yieldreduction_list.loc[i])
        y_values = [value*.01 for value in y_values]
        crop_variety = data.crop_variety.loc[i]
        
        
        fit_and_plot(x_values,y_values,crop_variety)
        
        






#     get crop_variations xy data from data
    
#     get crop_variations simple linear model
    
#     get crop_variations threshold slope model

#     get crop_variations breakpoint estimation model
    
#     plot all variations in one plots
    
#     create giant subplot. caluclate gridsize based on lenght of data
    
    
xy_crop_from_data(fn='problems.csv')

# # load pbest data
# data = pd.read_csv('problems.csv', sep=',')

# # X data. (independent) lists of irrigation water salinity values in ds/M 
# x_values=data.salinitylevel_list.apply(lambda x: ast.literal_eval(x))


# x_values=x_values[1]

# #Y Data (dependend). lists of values for relative yield
# y_values=data.yieldreduction_list.apply(lambda x: ast.literal_eval(x))
# # every item in every list gets converted from % to fraction
# y_values = y_values.apply(lambda x: [i*.01 for i in x])
# y_values=y_values[1]



# #fit for simple linear model
# x_simple,y_simple= calc_xy_values(x_values,y_values,function=threshold_slope,
#                                   method=simple_fit)


# #fit piecewise
# x_piecewise,y_piecewise = calc_xy_values(x_values,y_values,function= piecewise, 
#                                          method=simple_fit)


# #fit breakpoint
# x_bp,y_bp = calc_xy_values(x_values,y_values,function= piecewise,
#                            method=estimate_breakpoint)









    
    

#piecewise example

# https://stackoverflow.com/questions/41641880/using-scipy-curve-fit-with-piecewise-function



# curvefit initial guess ist p0. 
# sollte immer guessen im bereich zwischen dem letzten x über 100 und dem ersten unter 100 und
# einfach zwischen beiden werten in 0.1 schritten durch enumeraten und jedesmal den fehler messen und dann den besten nehmen



# p0=array_like, optional

#     Initial guess for the parameters (length N). If None, then the initial values will all be 1 (if the number of parameters for the function can be determined using introspection, otherwise a ValueError is raised).


# do this tutorial:
    # https://machinelearningmastery.com/curve-fitting-with-python/


# plot simple regression

# linear function
# f = -1*m*x

# scipy.optimize.curve_fit(f,X[0],Y[0])
# check NNLS
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.lsq_linear.html#scipy.optimize.lsq_linear




