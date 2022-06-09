#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 14:02:15 2019

@author: koenig.g
"""


#****************************#
# Code to compute the SPSA   #
# With an external code in   #
# Which we wanna store the   #
# Variables. By G. Koenig,   #
# The (08/07/2019)           #
# Modified to take spin-up in#
# -to account (G.Koenig)     #
# (14/08/2019)               #
# Modify to save some of the #
# Data in a json file        #
# (G. Koenig, 27/08/2019)    #
# Modified to use accelerated#
# Versions of CROCO (OPENMP) #
# (G. Koenig, 28/08/2019)    #
# Modify to have a tidal sign#
# -al as input               #
# By G.Koenig, the 30/08/2019#
# Modified with a cleaner    #
# Cost function (G.Koenig,   #
# 03/09/2019)                #
# Modified to use another ver#
# -sion of CROCO (03/10/2019)#
# Modified to use another    #
# External library to modify #
# The wave force (G.Koenig   #
# 26/08/2020 )               #
# Now we modify both friction#
# And the wave-breaking      #
# (G. Koenig, 15/10/2020)    #
# I modified the cost functio#
# -n and the breaking script #
#****************************#

#****Packages import*********#
import numpy as np
import random
import os
import shutil # To import the copy
import Cost_func_dataframe_modified_3 as Cost_func
# Because it is heavy to recopy
import Modif_wave_force_angle as Modif_wave_force # For the modification of the wave force
import Modif_fric_2_parts as Modif_fric# For the modification of the friction
import netCDF4
import subprocess
import threading # To parallelize computations
import time
import scipy.interpolate # For the spline approximations
import copy # To have the module deepcopy to copy dictionnaries
import json 
#***Variables declaration****#

#SPSA variables
Y = [] # An empty list to store the modifications
J,J_plus,J_minus = [],[],[]
grad = [] 
pert_vec = []
n_max = 50 # Number of iterations
L = [4] # Sizes of the perturbation vectors
c = np.array([1e-3,1e-3,1e-3,2e-4]) # Perturbation coefficient
a = np.array([1e-5,1e-5,1e-5,1e-5]) # Gain coefficient
dic_data_saver={'Wave_coeff':[],# Dictionnary to store the data
                'Wave_coeff_north':[],
		'Coral_fric':[],
		'Other_fric':[],
                'cost_func':[]}# Of the gradient descent
ctr = 0 # We need a counter for the for loop
ctr_json = 0 # And another counter for the json

#Outputs of model and data
data_fold = '../../DONNEES/DonneesPropres_2016/'
fold_van = '/home/koenig.g/SPSA_CROCO/CROCO_WF_CAMILLE/SPSA_CROCO_WF_FRIC_WF/ESSAI_OUANO_WF_ZETA' 
# Name of the original file
fold_0 = '/home/koenig.g/SPSA_CROCO/CROCO_WF_CAMILLE/SPSA_CROCO_WF_FRIC_WF/ESSAI_OUANO_WF_ZETA_0'
# Name of the principal_file
fold_plus = '/home/koenig.g/SPSA_CROCO/CROCO_WF_CAMILLE/SPSA_CROCO_WF_FRIC_WF/ESSAI_OUANO_WF_ZETA_+'
fold_minus = '/home/koenig.g/SPSA_CROCO/CROCO_WF_CAMILLE/SPSA_CROCO_WF_FRIC_WF/ESSAI_OUANO_WF_ZETA_-'
# And something for the outputs
outpath = '/SORTIE/his_CROCO_2D.00000.nc'
# Localization of the plots
fig_loc = '/home/koenig.g/SPSA_CROCO/CROCO_WF_CAMILLE/SPSA_CROCO_WF_FRIC_WF/PLOTS/'

# And now the locations of the files to be modified
WF_FILE = '/GRID/WF_init.nc'
GRD_FILE = '/GRID/GRD_reduced_smoothed_2m_trench.nc'
DRAG_FILE = '/GRID/DRAG_increased.nc'
# And the initial parameter for the wave force
Wave_coeff = 0.008
Wave_coeff_north = 0.008
# And for the friction
Fric_cor = 0.1
Fric_rest = 0.05
#Specifics for our problem
tide_comp = ['M2']#['M2','S2'] # The tidal component that we consider
nb_tide = 0 # An index to access correctly the tides in the forcing file

#***Functions declaration********#

def perturb_vector(c,L):
    """ We use this function to perturb vector in a symmetric fashion.
    
    INPUTS :
    c : Amplitude of the perturbations
    L : Size of the vector
    
    OUTPUTS :
    pert_vec : Vector of perturbations
    """
    
    pert_vec = c*(np.array([random.randint(0,1) for p in range(L)])-0.5)*2. 

         
    # We return the symetrically perturbed vector and the perturbation
    return pert_vec

def comp_grad(J_plus,J_minus,pert_vec):
    """ We use this functions to compute the gradient approximations.
    
    INPUTS :
    J_plus : First value (upward) of cost function
    J_minus : Second value (downward )of cost function
    pert_vec : Vector of perturbations
    
    OUTPUTS :
    grad : Gradient estimation """
    
    grad = (J_plus-J_minus)/(2.*pert_vec) # I take advantage of the fact that 
    # Python authorizes divisions by vectors
    
    return grad


def launch_model(fold) :
    """Launch the CROCO model in a particular sub folder
    
    INPUTS:
    fold : folder where we wann launch the CROCO instance
    
    OUTPUTS:
    None"""
    # Flag for execution
    out = False
    # We launch it
    cmd = "/home/koenig.g/SPSA_CROCO/CROCO_WF_CAMILLE/SPSA_CROCO_WF_FRIC_WF/launch_model.sh"
    subprocess.call([cmd,fold])
    # And we wait for it to execute
    while out==False :
        try :
            # We check if we are on the list of the cluster jobs
            cmd = "ssh koenig.g@cluster.osupytheas.fr squeue| awk '{print $4}' | grep koenig.g"
            name = subprocess.check_output(cmd, shell=True)
            # We go on waiting
            time.sleep(30)
            # And we print something to check
            print("Still runnin' bro' ")
        except subprocess.CalledProcessError as e :
            out=True
   
    return None

 # Function to store the output in json files    
def store_dic_data_saver(dic_data_saver,json_file,nb_it_max):
    """ Function to store the data into json file. First for compatibility's 
    sake we need to remove the numpy part of most of dictionnaries
    
    INPUTS :
    dic_data_saver : A dictionnary to store data
    json_file : The json file in which we are going to store the data
    
    OUTPUT : 
    A nice json file """
    
    # First thing we are going to create a copy of the dictionnary
    dic_data_saver_cp = copy.deepcopy(dic_data_saver)
    
    # Now we are going to denumpy some of the part of the array
    for i in range(nb_it_max):
        dic_data_saver_cp['Wave_coeff'][i]=\
            dic_data_saver_cp['Wave_coeff'][i].tolist()
        dic_data_saver_cp['Wave_coeff_north'][i]=\
            dic_data_saver_cp['Wave_coeff_north'][i].tolist()
        # And for the friction
        dic_data_saver_cp['Coral_fric'][i]=\
            dic_data_saver_cp['Coral_fric'][i].tolist()
        dic_data_saver_cp['Other_fric'][i]=\
            dic_data_saver_cp['Other_fric'][i].tolist()
    # And now we save it 
    with open(json_file, 'w') as fp:
        json.dump(dic_data_saver_cp, fp)

# Start of the script
print('Local time at the beginning is '+str(time.asctime()))

#  We initialize the vectors to store the perturbations
#Y.append(np.zeros(L[]))
Y.append(np.array([Wave_coeff,Wave_coeff_north,Fric_cor,Fric_rest]))
# Now we loop on all the sizes of perturbation vectors that we have
for k in L :
    
    # We reinitialize the gradient and the perturbations
    grad=np.zeros(k) 
    pert_vec=np.zeros(k)
    grad_estim = np.zeros(k)
    # We interpolate the latest solution
    Y_new=Y[-1]
    Y.append(Y_new)
    
    # Now I save it into my dictionnary
    dic_data_saver['Wave_coeff'].append(Y[-1][0])
    dic_data_saver['Wave_coeff_north'].append(Y[-1][1])
    dic_data_saver['Coral_fric'].append(Y[-1][2])
    dic_data_saver['Other_fric'].append(Y[-1][3])
    
    # We begin by estimating the first cost function 
    #and initializing perturbations
    print('Forcing vector is {}'.format(Y[-1]))
    Modif_wave_force.main(fold_0+GRD_FILE,fold_0+WF_FILE,Y[-1][0],Y[-1][1])
    Modif_fric.main(fold_0+GRD_FILE,fold_0+DRAG_FILE,Y[-1][2],Y[-1][3])
    launch_model(fold_0)
    cost_func = Cost_func.main(path_IS=data_fold,path_mod=fold_0+\
                               outpath,fig_loc=fig_loc,
                               tide_comp=tide_comp,n_iter=ctr_json)
    J.append(cost_func)
    
    # Saving the outputs
    dic_data_saver['cost_func'].append(J[-1])
    
    # We update the counter for the json file
    ctr_json+=1
    
    # And now, we can let the music do the looping
    for i in range(n_max):
	# We update the gain and the perturbation coefficient
        c_i = c/np.power((1 + i),0.101)
        a_i = a/np.power((50+i),0.602)
        # We generate the perturbations
        pert_vec=perturb_vector(c_i,k)
        print('The perturbed vector is {}'.format(Y[ctr*(n_max+1)+i+1]+pert_vec))
        
	# We add them
        Modif_wave_force.main(fold_plus+GRD_FILE,fold_plus+WF_FILE,Y[ctr*(n_max+1)+i+1][0]+pert_vec[0],Y[ctr*(n_max+1)+i+1][1]+pert_vec[1])
        Modif_wave_force.main(fold_minus+GRD_FILE,fold_minus+WF_FILE,Y[ctr*(n_max+1)+i+1][0]-pert_vec[0],Y[ctr*(n_max+1)+i+1][1]+pert_vec[1])
        # And we add the rest
        Modif_fric.main(fold_plus+GRD_FILE,fold_plus+DRAG_FILE,Y[ctr*(n_max+1)+i+1][2]+pert_vec[2],Y[ctr*(n_max+1)+i+1][3]+pert_vec[3])
        Modif_fric.main(fold_minus+GRD_FILE,fold_minus+DRAG_FILE,Y[ctr*(n_max+1)+i+1][2]-pert_vec[2],Y[ctr*(n_max+1)+i+1][3]-pert_vec[3])
        
	# And launch the two instances
        # First we define the threads
        threads = [threading.Thread(target=launch_model,args=(x,)) 
                   for x in [fold_plus,fold_minus]]
        # We open the processing space for them
        [thread.start() for thread in threads]
        # And we launch them in this processing space
        [thread.join() for thread in threads]
        
	# Once done, we compute the cost functions
        cost_func_plus = Cost_func.main(path_IS=data_fold,path_mod=fold_plus+\
                                        outpath,fig_loc=fig_loc,tide_comp=tide_comp,
					n_iter=ctr_json)
        cost_func_minus = Cost_func.main(path_IS=data_fold,path_mod=fold_minus+\
                                         outpath,fig_loc=fig_loc,tide_comp=tide_comp,
					 n_iter=ctr_json)
        # And we append them
        print('The perturbed cost functions are {}(+) and {}(-)'.format(cost_func_plus,cost_func_minus))
        J_plus.append(cost_func_plus)
        J_minus.append(cost_func_minus)
        
	# Now we can compute the gradient
        grad=comp_grad(J_plus[ctr*(n_max)+i],J_minus[ctr*(n_max)+i],
                              pert_vec)
        print('The gradient is worth {}'.format(grad))
	# And we update the vector of perturbations
        sign_grad = np.sign(J_plus[ctr*(n_max)+i] - J_minus[ctr*(n_max)+i])
        # And for that we use the previous data perturbations
        grad_estim = 0.8*grad_estim + sign_grad*pert_vec
        # And now we add the grad estimation
        Y.append(Y[ctr*(n_max+1)+i+1]-grad_estim)
        print('The use vector is {} and the final vector is {}'.format(Y[ctr*(n_max+1)+i+1],Y[-1]))
        # We then compute the cost function associated with those new perturbations
        Modif_wave_force.main(fold_0+GRD_FILE,fold_0+WF_FILE,Y[-1][0],Y[-1][1])
        Modif_fric.main(fold_0+GRD_FILE,fold_0+DRAG_FILE,Y[-1][2],Y[-1][3])
        launch_model(fold_0)
        cost_func = Cost_func.main(path_IS=data_fold,path_mod=fold_0+\
                                   outpath, fig_loc=fig_loc,tide_comp=tide_comp,
                                   n_iter=ctr_json)
        J.append(cost_func)
        
	# Now we save the outputs
        dic_data_saver['cost_func'].append(J[-1])
        dic_data_saver['Wave_coeff'].append(Y[-1][0])
        dic_data_saver['Wave_coeff_north'].append(Y[-1][1])
        dic_data_saver['Coral_fric'].append(Y[-1][2])
        dic_data_saver['Other_fric'].append(Y[-1][3])
        
	# We update the counter for the json file 
        ctr_json+=1
    
    # We increment the main counter
    ctr+=1    

# End of time computation
print('Local time at the end is is '+str(time.asctime()))

# And saving the stuff into a json file
store_dic_data_saver(dic_data_saver,'data.json',ctr_json)
