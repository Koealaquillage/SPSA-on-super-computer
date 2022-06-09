def main(GRD_FILE='../GRID/GRD_reduced.nc', DRAG_FILE='../GRID/DRAG_increased.nc', 
		COR_FRIC = 0.1, OTHER_FRIC=0.01,COR_FRIC_NORTH=1.) :
    """ Here we modify some grid files in order to determine the good wave-breaking parameter
    
    INPUTS : 
    GRD_FILE : The grid file, we need it for the coral reef
    DRAG_FILE : The drag file, with the friction in it
    COR_FRIC : The drag coefficient on the reef
    OTHER_FRIC : The drag coefficient elsewhere
    COR_FRIC_NORTH : The drag coefficient in the northern part of the reef

    OUTPUTS: Nothing, a modified file"""
    
    #********Packages import***********#
    import numpy as np
    import netCDF4 
                            
                            
    #********Now the files we want*****#
    #********To read*******************#
    file_grd = GRD_FILE
    file_drag = DRAG_FILE
  
    
    # And the value of wave-force we want to introduce
    drag_cor = COR_FRIC
    drag_other = OTHER_FRIC
    drag_cor_north = COR_FRIC_NORTH
    #****Now we have to import*********#
    #****The reef mask*****************#
    grid_nc=netCDF4.Dataset(file_grd,'r')
    reef_mask = grid_nc.variables['mask_reef'][:,:]
    grid_nc.close()
    
    #******And now we use it to modify***#
    #******The other mask****************#
    drag_nc=netCDF4.Dataset(file_drag,'r+')
    
    # I am trying a very little something #
    
    drag_nc.variables['bottom_drag_coef'][:,:]=drag_other
    # And for the southern part
    for i in range(0,155):
        for j in range(0,110):
            if reef_mask[i,j] == 1 :
                drag_nc.variables['bottom_drag_coef'][i,j] = drag_cor
    # And now I run something for the northern part
    for i in range(155,194):
        for j in range(0,110):
            if reef_mask[i,j] == 1 :
                drag_nc.variables['bottom_drag_coef'][i,j] = drag_cor_north

    return None

if __name__ == '__main__':
    main()
