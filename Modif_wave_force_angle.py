#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 13:27:01 2021

@author: koenig.g
"""

def main(GRD_FILE='../GRID/GRD_reduced_smoothed_2m_trench.nc', WAVE_FILE='../GRID/WF_init.nc', 
         WAVE_COEFF=0.083,WAVE_COEFF_NORTH=0.083,ANGLE=-15.) :
    """ Here we modify some grid files in order to determine the good wave-breaking parameter
    
    INPUTS : 
    GRD_file : The grid file, we need it for the coral reef
    WAVE_FILE : The wave file, where we put our parameters
    WAVE_COEFF: The wave coefficient of interested
    WAVE_COEFF_NORTH : The wave coefficient in the north of the domain
    ANGLE : The direction of the waves in degrees from the east direction 
            of the domain
    
    OUTPUTS: Nothing, a modified file"""
    
    #********Packages import***********#
    import numpy as np
    import netCDF4 
    
    #***We set a function to put the wave_force*#
    
    def set_wave_force(start,end,wave_coeff,angle,reef_mask,wave_nc,
                       var_wave_1,var_wave_2):
        
        for j in range(start,end):
            for k in range(0,60):
                if reef_mask[j,k]==1 and reef_mask[j,k-1]==0:
                    # And now we can use it to estimate the direction of waves
                    wave_nc.variables[var_wave_1][j,k]=wave_coeff*np.cos(angle)
                    # And the angle in the y direction
                    grad_y = reef_mask[j+1,k] - reef_mask[j-1,k]
                    # We have to set the radiation stress at different places
                    # Depending on the incidence angle of waves on the barreer
                    wave_nc.variables[var_wave_2][j,k]=wave_coeff*np.sin(angle)                         
                            
    #********Now the files we want*****#
    #********To read*******************#
    file_grd = GRD_FILE
    file_wave = WAVE_FILE
    
    var_grd='mask_reef'
    var_land='mask_rho'
    var_wave_1='FXwave'
    var_wave_2='FYwave'
    angle = np.deg2rad(ANGLE)
    
    #****We are not interested in the**#
    #***Entire domain, so we restrict**#
    #***It*****************************#
    y_st_1 = 2
    y_end_1 = 13
    
    y_st_2 = 24
    y_end_2= 121
    
    y_st_3 = 124
    y_end_3 = 144
    
    y_st_4 = 150
    y_end_4 = 180
    
    # And the value of wave-force we want to introduce
    wave_coeff = WAVE_COEFF
    wave_coeff_north = WAVE_COEFF_NORTH
    
    #****Now we have to import*********#
    #****The reef mask*****************#
    grid_nc=netCDF4.Dataset(file_grd,'r')
    reef_mask = grid_nc.variables[var_grd][:,:]
    mask_land = grid_nc.variables[var_land][:,:]
    grid_nc.close()
    
    #******And now we use it to modify***#
    #******The other mask****************#
    wave_nc=netCDF4.Dataset(file_wave,'r+')
    
    # I am trying a very little something #
    
    wave_nc.variables[var_wave_2][:,:] = 0.
    
    #*** And now we call the function four times
    set_wave_force(y_st_1,y_end_1,wave_coeff,angle,reef_mask,wave_nc,
                   var_wave_1,var_wave_2)
    set_wave_force(y_st_2,y_end_2,wave_coeff,angle,reef_mask,wave_nc,
                   var_wave_1,var_wave_2)
    set_wave_force(y_st_3,y_end_3,wave_coeff,angle,reef_mask,wave_nc,
                   var_wave_1,var_wave_2)
    set_wave_force(y_st_4,y_end_4,wave_coeff_north,angle,reef_mask,wave_nc,
                   var_wave_1,var_wave_2)
    
    #****And then I run the final loop
#    for j in range(0,194):
#        for k in range(0,110):
#            if mask_land[j,k] == 1 :
#                if wave_nc.variables[var_wave_2][j,k] < -wave_coeff :
#                    wave_nc.variables[var_wave_2][j,k] =0.
                    

    wave_nc.close()

    return None

if __name__ == '__main__':
    main()
