#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:06:18 2019

@author: koenig.g
"""

#*********************************************#
# Here we do a cleaner script to compute the  #
# Cost function, because the last one present-#
# Ed some bugs. By G. Koenig, the 04/12/2019  #
# Here we added possibility of choosing the   #
# Tidal components of interest. By G. Koenig  #
# The 20/01/2020                              #
# Modified to take into account the rotated   #
# Domain of CROCO. By G. Koenig,              #
# the 17/02/2020                              #
# Modified to treat u and v as independant    #
# Component instead of a tidal ellipse by     #
# Koenig, the 10/08/2020                      #
# Remodified to take the ellipse component,   #
# But only in the semi-major axis direction   #
# By G. Koenig, the 30/09/2020                #
# I modified the sensor position for the new  #
# Bathymetry. By G. Koenig the 01/02/2020     #
#*********************************************#


def main(path_IS = '/media/koenig.g/LaCie/1ERE_ANNEE/FEVRIER_18/15_02'+\
     '/Donnees_2016/DonneesPropres_2016/' ,path_mod = '/home/koenig.g/These/'+\
     '3EME_ANNEE/CROCO_WF/TEST_WF_YANN/Sortie/his_WF.nc',fig_loc = 'PLOTS/',
     tide_comp = 'auto',n_iter = 0,save_json = True) :
    """ This script is to compute a standard RMS cost function between data and
    CROCO's output. Note that we automated most of the time-reading and used panda
    dataframes to handle different timespans for different datas.

    INPUTS :
    path_IS : Path for the in-situ data ( with a very specific format)
    path_mod : Path for the output of CROCO
    fig_loc : Path for savings the outputs of models and data with some plots
    tide_comp : Tidal component of interest. In format ['M2',...], and ['M2'] if only one.
                if 'auto' is chosen, all the components are considered ( this is Utide syntax)
    n_iter : Argument for plotting the results. So that we can differentiate at which iteration
             of the particle filter we were
    save_json : Boolean to know if we are going to save our coefficients in json files
    OUTPUTS :
    cost_func : The computed cost function"""
    #*****Packages import*********#
    import numpy as np 
    import pandas as pd
    import netCDF4
    import scipy.io
    import matplotlib.pyplot as plt
    import Plotting_stuff_for_cost_func_1 # A function to plot stuff
    import utide
    import matplotlib.dates as mdates
    import json
    #******Variables declaration*********#
    
    # Temporary dictionnary for data
    temp_dic = {'Date': [], 'Station': [], 'Umean' :[], 'Vmean': [], 
                'Zeta': []}
    # Dataframes for the data and the model
    columns = ['Date','Station','Umean','Vmean','Zeta']
    df_IS = pd.DataFrame(columns = columns)
    df_mod = pd.DataFrame(columns = columns)
    
    freq = 'H' # Hour frequency for the data recording
    time_var = [] # To store a time array
    # List of the stations that we are going to use
    list_stat=['GDigo','Isie','Nord','Ouarai','PDigo','PtFixe','Recif1',
               'Recif2','Tenia']
    # List of the variables that we are going to use 
    list_var=['Zeta','Umean','Vmean']
    # Dictionnary to have the name equivalence between the model and the script
    dic_var_mod={'Zeta':'zeta',
                 'Umean':'ubar',
                 'Vmean':'vbar',
                 'lat_rho' : 'lat_rho'}
    # Now we have a dictionnary for the excluded data
    dic_exc_data = {'Zeta' : ['Ouarai','Tenia','Nord','GDigo','Isie'],
                     'Umean' : [], 'Vmean' : []}
    # Starting date
    start_date='2016-04-21 18:00:00'
    end_date = []
    # And the dictionnary for the positions of stations
    stat_pos={'GDigo': np.array([131,42]),
              'Isie':np.array([123,34]),
              'Nord':np.array([159,43]),
              'Ouarai':np.array([147,36]),
              'PDigo':np.array([127,52]),
              'PtFixe':np.array([100,48]),
              'Recif1':np.array([87,44]),
              'Recif2':np.array([82,40]),
              'Tenia':np.array([28,54])}
    # And a dictionnary for the latitude
    temp_lat_dic = {'GDigo':[],
		    'Isie':[],
		    'Nord':[],
		    'Ouarai':[],
		    'PDigo':[],
		    'PtFixe':[],
		    'Recif1':[],
		    'Recif2':[],
		    'Tenia':[]}
    stat_lat={'GDigo': np.array([]),
              'Isie':np.array([]),
              'Nord':np.array([]),
              'Ouarai':np.array([]),
              'PDigo':np.array([]),
              'PtFixe':np.array([]),
              'Recif1':np.array([]),
              'Recif2':np.array([]),
              'Tenia':np.array([])}

    # And the cost function
    cost_func = 0
    # The weight for the zeta, so that we equilibrate everything
    W_Zeta = 1.
    # The variable for the rotation of the domain
    rot_dom = np.deg2rad(45.)
    # And we have a dictionnary for the tidal coefficients 
    dic_coeff_model = {}
    dic_coeff_data = {}
    # And now we fill it
    for STAT in list_stat :
        dic_coeff_model[STAT] = {}
        dic_coeff_data[STAT] = {}

    # Now we create the import functions
    def data_import_netcdf(netid,stat,stat_pos,list_var,dic_var_mod,temp_dic):
        """ Here we used an already open netcdf file to load data into a 
            dictionnary preformatted.
            INPUTS :
            netid : Id of the opened netcdf file
            stat : station we are interested in
            stat_pos : tuple having the position of the station on the grid
            list_var : List of the variables we want to use
            dic_var_mod : Dictionnary for the equivalence between the name in
            the model and in the script
            
            OUTPUTS :
            temp_dic : Formatted dictionnary containing the data we are 
            interested in"""
            
        for VAR in list_var :
            # We quickly give shorter name to the stat_pos
            id_1,id_2 = stat_pos
            # We use the name equivalence
            VAR_mod=dic_var_mod[VAR]
            # And now we import the data
            if VAR_mod == 'lat_rho':
                temp_dic[VAR]=netid[VAR_mod][id_1,id_2]
            else :
                temp_dic[VAR]=netid[VAR_mod][:,id_1,id_2]
    
        # And the name of the station
        temp_dic['Station']=STAT
    
        # And we're done
        return temp_dic
    
    def data_import_IS(path_IS,stat,temp_dic) :
        """ Here we import the data that we have for the different stations 
            into a preformatted dictionnary.
        
        INPUTS :
        path_IS : Path to in-situ data
        stat : The station we are interested in
        
        
        OUTPUTS :
        temp_dic : The preformatted dictionnary to store our data"""
        
        # Local name of the file
        file_IS=path_IS+'Vit_'+stat+'.mat'
        # And we open it
        data_IS=scipy.io.loadmat(file_IS)
        # We have to treat a lil' the data. We use a panda data frame so
        # That we have more freedom in choosing the correct timing
        df_time = pd.DataFrame({'year':data_IS['Temps']['year'][0][0][:,0],
                           'month':data_IS['Temps']['month'][0][0][:,0],
                           'day':data_IS['Temps']['day'][0][0][:,0],
                           'hour':data_IS['Temps']['hour'][0][0][:,0],
                           'minute':data_IS['Temps']['minute'][0][0][:,0],
                           'second':data_IS['Temps']['year'][0][0][:,0]})
        # And we convert it to datetime
        df_time = pd.to_datetime(df_time)
        # We put it in the dictionnary
        temp_dic['Date']=df_time
        # We set all the rest in the dictionnary as well
        temp_dic['Umean'] = np.nanmean(data_IS['vitesse'][0][0][0][:],axis=1)
        temp_dic['Vmean'] = np.nanmean(data_IS['vitesse'][0][0][1][:],axis=1)
        temp_dic['Zeta'] = data_IS['P'][0,0][2][:,0]
        # And the name of the station
        temp_dic['Station']= stat
        # And we return
        return temp_dic
    
    
    #*******Now we do the import itself*************************************#
    # First for the netcdf. We beginning by importing it
    ind_netcdf = netCDF4.Dataset(path_mod)
    # And by creating a time variable on it
    size_time = ind_netcdf.dimensions['time'].size
    print('The number of time steps is {}'.format(size_time))
    time_var = pd.date_range(start = start_date,periods = size_time,freq=freq)
    #******Now we can import the model data************#   "
    for STAT in list_stat :
        # We begin with the latitude
        stat_lat[STAT] =data_import_netcdf(ind_netcdf,STAT,stat_pos[STAT],['lat_rho'],
                                           dic_var_mod,temp_lat_dic)[STAT]
                                           
        temp_dic = data_import_netcdf(ind_netcdf,STAT,stat_pos[STAT],list_var,
                                      dic_var_mod,temp_dic)
        # We add the time
        temp_dic['Date'] = time_var
        # And we put it in a temporary dataframe
        new_df = pd.DataFrame.from_dict(temp_dic)
        # And we put it int the dataframe
        df_mod = df_mod.append(new_df,ignore_index=True)
    
    # Always closing the netcdf
    ind_netcdf.close()
    
    #****Now we order a little the stuff so that we have access to it in the 
    #*** Good order
    df_mod = df_mod.set_index(['Station'])
    # And I sort it
    df_mod.sort_index(inplace=True)

    #*****And now we filter it to keep only the tidal component*****#
    for STAT in list_stat :
        # We explicit the variables first
        time = mdates.date2num(df_mod.loc[STAT].Date.values)
        zeta = df_mod.loc[STAT,'Zeta'].values
        u,v = df_mod.loc[STAT,'Umean'].values,df_mod.loc[STAT,'Vmean'].values
	# We modify them slightly to take into account the rotated domain
        u = u*np.cos(rot_dom) - v*np.sin(rot_dom)
        v = u*np.sin(rot_dom) + v*np.cos(rot_dom)
        # And now we can tidal filter them, by first extracting the
        # Coefficients
        coeff_zeta = utide.solve(time,zeta,lat=stat_lat[STAT],nodal=False,trend=False,
                                 method='robust',conf_int='linear',
                                 Rayleigh_min=.95, constit=tide_comp)
        coeff_u_v = utide.solve(time,u,v,lat=stat_lat[STAT],nodal=False,trend=False,
                                method='robust',conf_int='linear',Rayleigh_min=.95,
                                constit=tide_comp)

        # And now we fill our favorites dictionnaries 
        dic_coeff_model[STAT]['Amplitude elevation'] = coeff_zeta['A'].tolist()
        dic_coeff_model[STAT]['Phase elevation'] = coeff_zeta['g'].tolist()
        dic_coeff_model[STAT]['Mean elevation'] = coeff_zeta['mean'].tolist()
        dic_coeff_model[STAT]['Mean u'] = coeff_u_v['umean'].tolist()
        dic_coeff_model[STAT]['Mean v'] = coeff_u_v['vmean'].tolist()
        dic_coeff_model[STAT]['Semi-major velocity'] = coeff_u_v['Lsmaj'].tolist()
        dic_coeff_model[STAT]['Semi-major phase'] = coeff_u_v['g'].tolist()
        dic_coeff_model[STAT]['Semi-major angle'] = coeff_u_v['theta'].tolist()
        dic_coeff_model[STAT]['Semi-minor velocity'] = coeff_u_v['Lsmin'].tolist()
        
	# And then by reconstructing it
        tidal_zeta = utide.reconstruct(time,coeff_zeta)
        tidal_u_v = utide.reconstruct(time,coeff_u_v)
        df_mod.loc[STAT,'Zeta'] =  tidal_zeta.h
        df_mod.loc[STAT,'Umean'] =  tidal_u_v.u 
        df_mod.loc[STAT,'Vmean'] = tidal_u_v.v
        # Now I do extract the mean current direction and intensity
        dic_coeff_model[STAT]['Mean velocity'] = dic_coeff_model[STAT]['Mean u']**2 \
                                              + dic_coeff_model[STAT]['Mean v']**2
        dic_coeff_model[STAT]['Mean velocity'] = np.sqrt(dic_coeff_model[STAT]['Mean velocity'])
        dic_coeff_model[STAT]['Mean angle'] = np.arctan2(dic_coeff_model[STAT]['Mean v'],
                                                         dic_coeff_model[STAT]['Mean u'])


    #****And we import as well the in-situ data********#
    for STAT in list_stat :
        # We import the data in a dictionnary
        temp_dic = data_import_IS(path_IS,STAT,temp_dic)
        # And we put it in a temporary dataframe
        new_df = pd.DataFrame.from_dict(temp_dic)
        # And we put it int the dataframe
        df_IS = df_IS.append(new_df,ignore_index=True)
        
    #****And we do the same for the ordering and stuff
    df_IS = df_IS.set_index(['Station'])
    # And I sort it
    df_IS.sort_index(inplace=True)
    
    
    
    #****Now we can clean up the data a little bit********#
    for STAT in list_stat :
        df_IS.loc[STAT,'Umean']/=1000.
        df_IS.loc[STAT,'Vmean']/=1000.
    
        
    #***And we treat the data to have only the tidal signal***#
    for STAT in list_stat :
        # First we need to get the time signal in a good shape
        time = mdates.date2num(df_IS.loc[STAT].Date.values)
    
        # And now we get the signals
        zeta = df_IS.loc[STAT,'Zeta'].values
        u,v = df_IS.loc[STAT,'Umean'].values, df_IS.loc[STAT,'Vmean'].values
        # First we get the coefficients
        coeff_zeta=utide.solve(time,zeta,lat=stat_lat[STAT],nodal=False,trend=False,
                              method='robust',conf_int='linear',
                              Rayleigh_min=.95, constit=tide_comp)
        coeff_u_v=utide.solve(time,u,v,lat=stat_lat[STAT],nodal=False,trend=False,
                              method='robust',conf_int='linear',
                              Rayleigh_min=.95,constit=tide_comp)
        # And now we use only the tidal component
        tidal_u_v = utide.reconstruct(time,coeff_u_v)
        tidal_zeta = utide.reconstruct(time,coeff_zeta)
        df_IS.loc[STAT,'Zeta']= tidal_zeta.h
        df_IS.loc[STAT,'Umean'] = tidal_u_v.u
        df_IS.loc[STAT,'Vmean'] = tidal_u_v.v
        # And now we fill our favorites dictionnaries 
        dic_coeff_data[STAT]['Amplitude elevation'] = coeff_zeta['A'].tolist()
        dic_coeff_data[STAT]['Phase elevation'] = coeff_zeta['g'].tolist()
        dic_coeff_data[STAT]['Mean elevation'] = coeff_zeta['mean'].tolist()       
        dic_coeff_data[STAT]['Mean u'] = coeff_u_v['umean'].tolist()
        dic_coeff_data[STAT]['Mean v'] = coeff_u_v['vmean'].tolist()
        dic_coeff_data[STAT]['Semi-major velocity'] = coeff_u_v['Lsmaj'].tolist()
        dic_coeff_data[STAT]['Semi-minor velocity'] = coeff_u_v['Lsmin'].tolist()
        dic_coeff_data[STAT]['Semi-major phase'] = coeff_u_v['g'].tolist()
        dic_coeff_data[STAT]['Semi-major angle'] = coeff_u_v['theta'].tolist()

        # And then by reconstructing it
        tidal_zeta = utide.reconstruct(time,coeff_zeta)
        tidal_u_v = utide.reconstruct(time,coeff_u_v)
        df_IS.loc[STAT,'Zeta'] =  tidal_zeta.h
        df_IS.loc[STAT,'Umean'] =  tidal_u_v.u
        df_IS.loc[STAT,'Vmean'] = tidal_u_v.v
        # Now I do extract the mean current direction and intensity
        dic_coeff_data[STAT]['Mean velocity'] = dic_coeff_data[STAT]['Mean u']**2 \
                                                + dic_coeff_data[STAT]['Mean v']**2
        dic_coeff_data[STAT]['Mean velocity'] = np.sqrt(dic_coeff_data[STAT]['Mean velocity'])
        dic_coeff_data[STAT]['Mean angle'] = np.arctan2(dic_coeff_data[STAT]['Mean v'],
                                                         dic_coeff_data[STAT]['Mean u'])
        # And now we want to remove the mean out of that 
        df_IS.loc[STAT,'Zeta'] -= coeff_zeta['mean']
    #***And we exclude the data we do not want***********#
    for VAR in list_var :
        for STAT in dic_exc_data[VAR] :
            df_IS.loc[STAT,VAR] = np.nan

    #****And now we can create the cost function in itself*****#
    # Before that we need the time for the termination of modelling
    
    for STAT in list_stat :
        # Time variables
        t_i,t_f = time_var[0],time_var[-1]
        
        # If there are less in-situ data, we shorten it
        # Note that we can do comparisons because they are all timestamps
        if df_IS.loc[STAT].Date.iloc[-1] < t_f :
            t_f = df_IS.loc[STAT].Date.iloc[-1]
        # And if there is more we add things in it
        if df_IS.loc[STAT].Date.iloc[0] >  t_i : 
            t_i = df_IS.loc[STAT].Date.iloc[0]
        # I'd like to see the dates
        print('Initial time is {}'.format(t_i))
        print('Final time is {}'.format(t_f))
	# Ok, now we can compare the cost functions 
        # First we start with zeta 
        zeta_cost_func = 0
        for i in range(len(tide_comp)):
            # And now the local cost function
            loc_cost_func = np.cos((dic_coeff_data[STAT]['Phase elevation'][i] - dic_coeff_model[STAT]['Phase elevation'][i])*2.*np.pi/360.)
            loc_cost_func*= -2*(dic_coeff_data[STAT]['Amplitude elevation'][i]*dic_coeff_model[STAT]['Amplitude elevation'][i])
            loc_cost_func += dic_coeff_data[STAT]['Amplitude elevation'][i]**2
            loc_cost_func += dic_coeff_model[STAT]['Amplitude elevation'][i]**2
            # We add it to the zeta cost func
            zeta_cost_func += loc_cost_func
            # And we add it to the total cost function
            cost_func+= loc_cost_func
        # Now we can do the plotting
        # First we perform the data extraction and comparison
        data_IS = df_IS.loc[STAT].set_index(['Date'])
        data_IS = data_IS['Zeta'].loc[t_i:t_f]
        #Ok, so we cannot directly interpolate it so we resample at higher
        # Frequency
        data_IS = data_IS.resample('H').pad()
        # We select the good dates
        data_mod = df_mod.loc[STAT].set_index(['Date'])
        data_mod = data_mod['Zeta'].loc[t_i:t_f]
        # Now the plotting part
        Plotting_stuff_for_cost_func_1.plotting_stuff(data_mod,data_IS,
                                                    STAT,'zeta',zeta_cost_func,
                                                    fig_loc,n_iter)

        # Ok, now we do all the same for u, but we must add the constant terms

        Lsmaj_cost_func = 0
        # So we begin with the constant terms
        if STAT == 'OUARAI' :
            # We have a strange reading in Ouarai so we exclude it
            loc_cost_func =0.
        else :
            loc_cost_func = -2.*dic_coeff_data[STAT]['Mean velocity']*dic_coeff_model[STAT]['Mean velocity']
            loc_cost_func *= np.cos((dic_coeff_data[STAT]['Mean angle'] - dic_coeff_model[STAT]['Mean angle'])*2.*np.pi/360.)
            loc_cost_func += dic_coeff_data[STAT]['Mean velocity']**2
            loc_cost_func += dic_coeff_model[STAT]['Mean velocity']**2
        # And we add it to the LSmaj cost func
        Lsmaj_cost_func += loc_cost_func
        # But also to the total cost func
        cost_func += loc_cost_func

        for i in range(len(tide_comp)):
            # And now the local cost function
            loc_cost_func = np.cos(2.*np.pi/360.*((dic_coeff_data[STAT]['Semi-major angle'][i] - dic_coeff_model[STAT]['Semi-major angle'][i])-(dic_coeff_data[STAT]['Semi-major phase'][i] - dic_coeff_model[STAT]['Semi-major phase'][i])))
            loc_cost_func*= -2*(dic_coeff_data[STAT]['Semi-major velocity'][i]*dic_coeff_model[STAT]['Semi-major velocity'][i])
            loc_cost_func += dic_coeff_data[STAT]['Semi-major velocity'][i]**2
            loc_cost_func += dic_coeff_model[STAT]['Semi-major velocity'][i]**2
            # We add it to the u cost func
            Lsmaj_cost_func += loc_cost_func
            # And we add it to the total cost function
            cost_func += loc_cost_func
            # Now we can do the plotting
            # First we perform the data extraction and comparison
            data_IS = df_IS.loc[STAT].set_index(['Date'])
            data_IS = data_IS['Umean'].loc[t_i:t_f]
            #Ok, so we cannot directly interpolate it so we resample at higher
            # Frequency
            data_IS = data_IS.resample('H').pad()
            # We select the good dates
            data_mod = df_mod.loc[STAT].set_index(['Date'])
            data_mod = data_mod['Umean'].loc[t_i:t_f]
            # Now the plotting part
            Plotting_stuff_for_cost_func_1.plotting_stuff(data_mod,data_IS,
                                                        STAT,'ubar',Lsmaj_cost_func,
                                                        fig_loc,n_iter)
            # And we're going to do the same for vbar
            # First we perform the data extraction and comparison
            data_IS = df_IS.loc[STAT].set_index(['Date'])
            data_IS = data_IS['Vmean'].loc[t_i:t_f]
            #Ok, so we cannot directly interpolate it so we resample at higher
            # Frequency
            data_IS = data_IS.resample('H').pad()
            # We select the good dates
            data_mod = df_mod.loc[STAT].set_index(['Date'])
            data_mod = data_mod['Vmean'].loc[t_i:t_f]
            # Now the plotting part
            Plotting_stuff_for_cost_func_1.plotting_stuff(data_mod,data_IS,
                                                        STAT,'vbar',Lsmaj_cost_func,
                                                        fig_loc,n_iter)
        print(cost_func)
    # And we save the coefficients if we decided so
    if save_json :
        # First for the model outputs
        with open('M2_coeff_model.json', 'a') as dump_file :
            json.dump(dic_coeff_model, dump_file)
        # Then for the IS data
        with open('M2_coeff_data.json', 'a') as dump_file :
            json.dump(dic_coeff_data, dump_file)
    
    return cost_func

if __name__ == '__main__' :
    main()
