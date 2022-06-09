#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:48:48 2019

@author: koenig.g
"""
#*************************#
# External function for   #
# Plotting stuff          #
# During the cost function#
# Computation, by         #
# G. Koenig the 06/12/2019#
#*************************#

#****Required packages******#
import matplotlib.pyplot as plt
import pandas as pd


def plotting_stuff(ts_mod,ts_IS,STAT,VAR, cost_loc,fig_loc,n_iter):
    """ Function to plot velocity and overelevation and save them in some plots
    
    
    INPUTS : 
    ts_mod,ts_IS : Overelevation and velocities for the model and in-situ data 
                   (m and m.s-1), in pandas timeserie so that also include
                   the dates
    STAT : Name of the station
    VAR : List with the name of the variables
    cost_loc : Local cost function, to serve as a rms error
    fig_loc : Place to save the figure
    n_iter = Number to add to the figure so that we can see how it evolved
    OUTPUTS :
    None"""
    
    # We are going to loop on the variables
    # We declare the figure and the axis
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    # We plot the stuff
    ax.plot(ts_IS,label='In situ data')
    ax.plot(ts_mod,label = 'Model outputs')
    ax.tick_params(axis='x', rotation=-75)
    # Here it will vary with the variable
    if VAR=='Zeta' :    
        ax.set_ylabel('Overelevation (m)')
        ax.set_title(VAR+str(' at ')+STAT+'\n'+' RMS error : ' \
                     + str(cost_loc) +'m^2')
    else :
        ax.set_ylabel('Velocity (m.s-1)')
        ax.set_title(VAR+str(' at ')+STAT+'\n'+' RMS error : ' \
                     + str(cost_loc) +'m^2.s^-2')
    
    ax.legend(loc='best')
    # And now we save it in png
    fig.savefig('{}{}_{}_{}.png'.format(fig_loc,STAT,VAR,n_iter),bbox_inches='tight')
    # And we close it
    plt.close()
    return None
