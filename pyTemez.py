# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 17:16:41 2024

@author: Hector Macian-Sorribes
"""


# Temez model function

def Temez_single_stage(meteorology, parameters, initial_conditions):
    
    # This function performs the Temez hydrological model for given meteorological inputs (P,ET0)       
    # Meteorological inputs must be a list with [P,ETP] values, both in mm
    # Parameter inputs should be a 5-element list in the form [Hmax,C,Imax,alpha,area]
        # Hmax in mm
        # C adimensional
        # Imax in mm
        # alpha in months^-1
        # (if more than one discharge branch is used, it should be a list with the alpha coefficients
        # per branch and the distribution coefficient of each branch)
        # area in Km^2
        
    # initial_conditions inputs should be a 2-elements vector in the form [H0,V0]
        # H0 in mm
        # Vo in Mm3 (one value per branch)
    # The code returns the hydrological response in Mm3 and the final conditions in the same format as initial ones
    
    # Separation of meteorological inputs
    
    P = meteorology[0]
    PET = meteorology[1]

    # Separation of input parameters  
    
    Hmax = parameters[0]
    C = parameters[1]
    Imax = parameters[2]
    area = parameters[4]
    
    # The alpha is sepparated in a different way: if there is one branch, then the distribution coefficient
    # is given a value of 1. If there is more than 1 branch, then the fourth element is spllited into
    # alphas and distribution coefficients
    
    aquifer_params = parameters[3]
    
    if type(aquifer_params) == float:
        
        alpha = aquifer_params
        distribution_coefficients = [1]
        
    elif type(aquifer_params) == int:
        
        alpha = aquifer_params
        distribution_coefficients = [1]
        
    else:
        
        split = int(len(aquifer_params)/2)
    
        alpha = aquifer_params[:split]
        distribution_coefficients = aquifer_params[split:]    
    
#    Separation of initial conditions

    if len(initial_conditions) == 2:

        Hini = initial_conditions[0]
    
        Vini = initial_conditions[1]
        
    else:
        
        Hini = initial_conditions[0]
    
        Vini = initial_conditions[1:]        
    
    # Calculus of P0    
    
    P0 = C*(Hmax - Hini)
    delta = Hmax - Hini + PET
    
    # Calculation of rainfall excess
    
    if P <= P0:
        
        T = 0
            
    else:
        
        T = ((P-P0)**2)/(P+delta-2*P0)

    # Calculation of evapotranspiration
        
    ET = min(PET , Hini + P - T)
    
    # Calculation of final soil moisture
    
    H = max(0 , Hini + P - T- ET)
        
    # Calculation of infiltration rates
    
    I = Imax * (T)/(T + Imax)
        
    # Calculation of direct surface runoff
    
    Qsurface = (T - I) * area / 1000  
    
    # Calculation of aquifer recharge
    
    if type(aquifer_params) == float:
        
        R = I * area / 1000
        
    elif type(aquifer_params) == int:
        
        R = I * area / 1000
        
    else:
        
        R = list(range(len(alpha)))
        V = list(range(len(alpha)))
        Qsubt = list(range(len(alpha)))
        
        for cell in range(len(alpha)):
            
            R[cell] = I * distribution_coefficients[cell] * area / 1000
        
    # Calculation of final storage and discharge from the aquifer cells
    
    if type(aquifer_params) == float:
        
        V = Vini * np.exp(-alpha) + R / alpha * (1 - np.exp(-alpha)) 
        Qsubt = Vini - V + R
        Qtotal = Qsurface + Qsubt
        
    elif type(aquifer_params) == int:
        
        V = Vini * np.exp(-alpha) + R / alpha * (1 - np.exp(-alpha))  
        Qsubt = Vini - V + R
        Qtotal = Qsurface + Qsubt
        
    else:  
        
        for cell in range(len(alpha)): 
            
            V[cell] = Vini[cell] * np.exp(-alpha[cell]) + (R[cell] / alpha[cell]) * (1 - np.exp(-alpha[cell])) 
            Qsubt[cell] = Vini[cell] - V[cell] + R[cell]
            Qtotal = Qsurface + np.sum(Qsubt)
   
    return(Qtotal, H, V)


# Function solving the Temez model on a given time series

    # It requires:
        
        # (1) The precipitation time series in mm
        # (2) The temperature time series in degree Celsius
        # (3) The temperatures of the previous twelve months, to calculate the PET at the first time stages
        # (4) The model parameters in the following way:
        #       - Hmax in mm
        #       - C adimensional
        #       - Imax in mm
        #       - alpha in months^-1
        #         (if more than one discharge branch is used, it should be a list with the alpha coefficients
        #         per branch and the distribution coefficient of each branch)
        #       - area in Km^2            
        # (5) The Initial conditions in the following way:
        #       - H0 in mm
        #       - Vo in Mm3 (one value per branch)
        # (6) PET correction factor
        
    # The code returns the hydrological response in Mm3            
        

def Temez_time_series(P_time_series, ET0_time_series, parameters, initial_conditions):
            
    time_steps = list(range(len(P_time_series)))
    
    Q = list(range(len(P_time_series)))
    
    # Start time loop
    
    for time_step in time_steps:
        
    # Calculation of dicharges
        
        [Q[time_step], initial_conditions[0], initial_conditions[1]] = Temez_single_stage([P_time_series[time_step], ET0_time_series[time_step]] , parameters, initial_conditions)
             
    return(Q)


############################
### SERPIS RIVER EXAMPLE ###
############################

parameters_Serpis = [   [692.97, 0.04, 231.70, 0.73, 468.60],
                        [791.22, 0.01, 41.77, 0.38, 108.24],
                        [1186.65, 0.05, 35.84, 0.51, 149.37]  ]

initial_conditions_Serpis = [    [0.00, 6.94],
                                 [16.41, 2.90],
                                 [124.21, 5.93]   ]

models = ['GFDL_ESM4','IPSL_CM6A-LR', 'MPI_ESM1_2_HR', 'MRI_ESM2_0', 'UKESM1_0_11']
scenarios = ['ssp126', 'ssp585']
subbasins = ['Beniarres', 'Encantada', 'Vernissa']

for model in models:
    
    for scenario in scenarios:

        name_pr = 'PATH TO MASTER RAINFALL FILE' + model + '_' + scenario + '_pr_2015_2100_Serpis.csv'
        name_ET0 = 'PATH TO MASTER ET0 FILE' + model + '_' + scenario + '_ET0_2015_2100_Serpis.csv'
        
        pr_dataframe = pd.read_csv(name_pr, sep = ';', index_col = 0) 
        ET0_dataframe = pd.read_csv(name_pr, sep = ';', index_col = 0) 
        
        for subbasin_number in range(len(subbasins)):
            
            if subbasin_number == 0:
                
                output_dataframe = pr_dataframe[pr_dataframe.keys().values[0 : 2]]
        
            pr_time_series = pr_dataframe[subbasins[subbasin_number]].values
            ET0_time_series = ET0_dataframe[subbasins[subbasin_number]].values
            
            output_dataframe[subbasins[subbasin_number]] = Temez_time_series(pr_time_series, ET0_time_series, parameters_Serpis[subbasin_number], initial_conditions_Serpis[subbasin_number])
            
        name_output = 'PATH TO OUTPUT FOLDER' + model + '_' + scenario + '_streamflow_2015_2100_Serpis.csv'
        
        output_dataframe.to_csv(name_output, sep=';')  
        
