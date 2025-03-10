#!/usr/bin/env python
# coding: utf-8

# In[1]:


# AMG 510
from pysb import *
def load_model_G12C_NEA():
    Model()
    Monomer('GEF', ['b'])# Active Guanine nucleotide exchange factor ex SOS,
    Monomer('GTP', ['b']) # Guanosine triphosphate
    Monomer('GDP', ['b']) # Guanosine diphosphate
    Monomer('muRASn' ,['b'])# mutated free RAS
    Monomer('muRAS' ,['b', 'state'], {'state': ['T', 'D']})
    Monomer('I', ['b']) # Inhibitor of the KRAS protein
    Monomer('muRASI', ['b']) # the complex of mutated G12C binded with inhibitor
    Monomer('muRAS_GDP_GEF_I_inter', ['b']) # the intermedia complex of mutated G12C binded with GEF
    Monomer('muRAS_GDP_GEF_I', ['b']) # the complex of mutated G12C binded with GEF and inhibitor

    #input rate constans-----------------------MUTATION----------------(need to recalculate the GEF GAP association and dissociation rate)
    # RAS intrinsic association and dissociation to GDP and GTP
    Parameter('mukd_GDP', 1.1e-4)# intrinsic dissociation of RAS_GDP
    Parameter('muka_GDP', 2.3e6)# intrinsic association of RAS and GDP (the value shown already times the conc. of GDP 1.8e-5) (2.3e6)
    Parameter('mukd_GTP', 2.5e-4)# intrinsic dissociation of RAS_GTP
    Parameter('muka_GTP', 2.2e6)# intrinsic association of RAS and GTP (the value shown already times the conc. of GTP 1.8e-4) (2.2e6)
    
    # RAS interchanging between the states binding GDP and GTP through GEF
    #Km_GEF_GDP = 3.86e-4
    Parameter('mukDa_GEF', 3400) #(kDd_GEF + kcat_GDP)/Km_GEF_GDP = (10+3.9)/3.86e-4=36010.36269 # GEF associate with RAS_GDP# 3400*250 = 8500000
    Parameter('mukDd_GEF', 9.224) # We first assign 10 # Dissociation of GEF_RAS_GDP
    Parameter('mukcat_GDP', 21667) # RASGDP to RASGTP catalyzed by GEF 3.9/1.8e-4=21667
    #Km_GEF_GTP = 3e-4
    Parameter('mukTa_GEF', 3000) # (kTd_GEF+kcat_GTP GEF)/Km_GEF_GTP= 3e-1+7.2e-1/3e-4=3400associate with RAS_GTP# 3000*250 = 750000
    Parameter('mukTd_GEF', 1.8e-1) # We first assign3.1e-1# Dissociation of GEF_RAS_GTP
    Parameter('mukcat_GTP', 40000) # RASGTP to RASGDP catalyzed by GEF 7.2e-1/1.8e-5=40000
    #GAP constant from: Ahmadian, Mohammad Reza, et al. "Individual rate constants for the interaction of Ras proteins with GTPase-activating proteins determined by fluorescence spectroscopy." Biochemistry 36.15 (1997): 4535-4541.
    # RAS intrinsic GTPase activity
    Parameter('mukhyd',  2.522058824e-4)# intrinsic GTPase activity of RAS# 2.522058824e-4

    Parameter('ka_RAS_I', 1e5)# here we first set the value to zero, but we will asign a value accrodingly later
    Parameter('kd_RAS_I', 7.73)# here we first set the value to zero, but we will asign a value accordingly later
    # RAS form a irreversible covalent bond with inhibitor
    Parameter('kon_RAS_I', 0.85) #/s This kiact value is from the literature[2] is kinact. #0.85
    Parameter('kon_RASn_I', 0.6) #/s drug forming bond with free RAS with 10 time lower rate than the GDP bound G12C

    # input the rules------------------MUTATION--------------------------
    # RAS intrinsic association and dissociation to GDP and GTP
    Rule('muRAS_bind_unbind_GDP', muRASn(b=None) + GDP(b=None)| muRAS(b=None, state='D'), muka_GDP, mukd_GDP) # RAS bind and unbind GDP
    Rule('muRAS_bind_unbind_GTP', muRASn(b=None) + GTP(b=None)| muRAS(b=None, state='T'), muka_GTP, mukd_GTP) # RAS bind and unbind GTP
    # RAS intrinsic GTPase activity
    Rule('muRAS_GTPase_activity', muRAS(b=None, state='T') >> muRAS(b=None, state='D'), mukhyd) # RAS intrinsic GTPase activity
    # RAS interchanging between the states binding GDP and GTP through GEF
    Rule('muRAS_GDP_bind_unbind_GEF', muRAS(b=None, state='D') + GEF(b=None) | muRAS(b=1, state='D') % GEF(b=1), mukDa_GEF, mukDd_GEF) # RAS_GDP bind and unbind GEF
    Rule('muRAS_GDP_to_GTP_by_GEF', muRAS(b=1, state='D') % GEF(b=1) + GTP(b=None) >> muRAS(b=None, state='T') + GEF(b=None) + GDP(b=None), mukcat_GDP) # RAS_GDP to RAS_GTP catalyzed by GEF
    Rule('muRAS_GTP_bind_unbind_GEF', muRAS(b=None, state='T') + GEF(b=None) | muRAS(b=1, state='T') % GEF(b=1), mukTa_GEF, mukTd_GEF) # RAS_GTP bind and unbind GEF
    Rule('muRAS_GTP_to_GDP_by_GEF', muRAS(b=1, state='T') % GEF(b=1) + GDP(b=None) >> muRAS(b=None, state='D') + GEF(b=None) + GTP(b=None), mukcat_GTP) # RAS_GTP to RAS_GDP catalyzed by GEF
    
    # inhibitor bind muRAS-GDP (here we simulate the two-step covalent binding)
    Rule('muRAS_GDP_associate_dissociate_I', muRAS(b=None, state='D') + I(b=None) | muRAS(b=1, state='D') % I(b=1), ka_RAS_I, kd_RAS_I) # muRAS:GDP form a reversible binding with inhibitor
    Rule('muRAS_GDP_bind_I', muRAS(b=1, state='D') % I(b=1) >> muRASI(b=None), kon_RAS_I) # muRAS:GDP bind with inhibitor

    # inibitor can also target free muRAS G12C with a lower rate
    #Rule('muRASn_associate_I', muRASn(b=None) + I(b=None) | muRASn(b=1) % I(b=1), ka_RAS_I, kd_RAS_I ) # muRASn bind with inhibitor
    #Rule('muRASn_bind_I', muRASn(b=1) % I(b=1) >> muRASI(b=None), kon_RASn_I) # muRASn bind with inhibitor

    # inhibitor can also target the complex of muRAS G12C and GEF with a lower rate
    #Rule('muRAS_GDP_GEF_associate_dissociate_I', muRAS(b=None, state='D') % GEF(b=None) + I(b=None) | muRAS_GDP_GEF_I_inter(b=None), ka_RAS_I, kd_RAS_I) # muRAS:GDP:GEF form a reversible binding with inhibitor
    #Rule('muRAS_GDP_GEF_bind_I', muRAS_GDP_GEF_I_inter(b=None) >> muRAS_GDP_GEF_I(b=None), kon_RASn_I) # muRAS:GDP:GEF bind with inhibitor
    
    # muRAS=9uL, GEF=3uL, 

    # Steps for NEA exchage assay (final volume is 15 uL)
    # 1. Pre-incubate drug with WT(5 mins) or G12C(20 mins) in 9uL. G12C stock concentration 25nM (final conc. 15nM). drug conc. 0~16.7uM(final conc. 0~10uM)
    # 2. read out the all the conc. and dilute to final volumen(15uL)
    # 3. Eu3+-GTP(final conc. 10 nM)
    # 4. SOScat(final conc. 5 nM)
    # 5. Measurement time: 10 min after mixing all the solutions.
    #input the initial conditions (here note that we only simulate the mutated G12C RAS)
    Parameter('muRAS_GDP_0', 2.5e-8) # concentration of muRAS:GDP added # 2.5e-8
    #Parameter('muRAS_0', 5e-9) # concentration of free muRAS added # 0
    Parameter('GDP_0', 5e-7) # concentration of GDP added # 0
    Parameter('GEF_0', 0) # concentration of GEF added(5e-9)
    Parameter('I_0', 0)# concentration of inhibitor added(1e-5)
    Parameter('GTP_0', 0) # concentration of GTP added(1e-8)

    Initial(muRAS(b=None, state='D'), muRAS_GDP_0)
    Initial(GEF(b=None), GEF_0)
    Initial(I(b=None), I_0)
    Initial(GTP(b=None), GTP_0)
    Initial(GDP(b=None), GDP_0)
    #Initial(muRASn(b=None), muRAS_0)

    # Observables (outputs)
    Observable('obsmuRAS_GTP', muRAS(b=None, state='T'))
    Observable('obsmuRAS_GDP', muRAS(b=None, state='D'))
    Observable('obsGEF', GEF(b=None))
    Observable('obsmuRAS_GDP_GEF', muRAS(b=1, state='D') % GEF(b=1))
    Observable('obsmuRAS_GTP_GEF', muRAS(b=1, state='T') % GEF(b=1))
    Observable('obsGDP', GDP(b=None))
    Observable('obsGTP', GTP(b=None))
    Observable('obsI', I(b=None))
    Observable('obsmuRASI', muRASI(b=None))
    Observable('obsmuRASn', muRASn(b=None))
    Observable('obsmuRASI_inter', muRAS(b=1, state='D') % I(b=1))

# main function
from pysb.simulator import ScipyOdeSimulator
#import pylab as pl
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


# In[6]:


## function to simulate the experiment protocol for NEA test
def NEA_AMG510(i_conc, gdp_conc, gtp_conc):
    # First we simulate the solution of protein KRAS G12C after the protein purification: theoretically, the KRASG12C should all bound with GDP. 
    # but the GDP can naturally dissociate from the protein so we need to run a period of time to let the system reach the equilibrium
    load_model_G12C_NEA()#
    I_0.value = i_conc
    GDP_0.value = gdp_conc
    t = np.linspace(0, 1200, 1201)
    simres = ScipyOdeSimulator(model, tspan = t, compiler= 'python').run()
    df = simres.dataframe
    df_clean = pd.DataFrame()
    for obs_name in df.keys():
        if obs_name.startswith('obs'):
            df_clean[obs_name] = df[obs_name]
    my_list = [str(i) for i in model.species]

    # Second the protein-inhibitor mixture further mix with GEF and GTP for 60 minutes
    restart_value = simres.species[-1]*(9/15) # take the last value from the last run
    restart_value[my_list.index('GEF(b=None)')] = 5e-9 # 5nM SOS
    restart_value[my_list.index('GTP(b=None)')] = gtp_conc # 10nM GTP
    t2 = np.linspace(1200, 4800, 3601)
    simres2 = ScipyOdeSimulator(model, tspan = t2, initials=restart_value, compiler="python").run()
    df = simres2.dataframe
    df_clean2 = pd.DataFrame()
    for obs_name in df.keys():
        if obs_name.startswith('obs'):
            df_clean2[obs_name] = df[obs_name]
    my_list = [str(i) for i in model.species]
    combined_df = pd.concat([df_clean, df_clean2.iloc[1:]], axis=0, ignore_index=True)
    # concatenate the two dataframes
    #return restart_value, simres.species[-1]
    return df_clean2['obsmuRAS_GTP'].iloc[-1]


# In[7]:


def search_amg_ic50(gdp_conc, gtp_conc):
    max_signal = NEA_AMG510(0*(15/9), gdp_conc*(15/9), gtp_conc)
    target_value = max_signal*0.5

    # binary search
    low = 0
    high = 1e-5
    # set up the tolerance
    tolerance = 0.001

    max_iterations = 100  # Maximum number of iterations
    iterations = 0  # Initialize iteration counter

    while iterations < max_iterations:
        # choose the model
        upper_lim = NEA_AMG510(high*(15/9), gdp_conc*(15/9), gtp_conc)
        lower_lim = NEA_AMG510(low*(15/9), gdp_conc*(15/9), gtp_conc)
        mid= NEA_AMG510((high + low)*(15/9)/2, gdp_conc*(15/9), gtp_conc)
        # find the new range
        if mid > target_value:
            low = (high + low)/2
        else:
            high = (high + low)/2
        # calculate the difference between the target ratio and the mid ratio
        difference = abs(mid - target_value)/target_value
        # check the tolerance
        if difference <= tolerance:
            #print('IC50 found as', (upper_bound + lower_bound)/2)
            #print('The ratio is', mid_ratio)
            conc_IC50 = (high + low)/2
            #ic50 = ic50.append({'Kd': 7.73/on_rate, 'non_co_inact_ic50': conc_IC50}, ignore_index=True)
            print('IC50 found as', conc_IC50)
            AMG510_IC50 = conc_IC50
            break
        iterations += 1
        print('iteration' + str(iterations))
        if iterations == max_iterations:
            #print('excceed iteration')
            #ic50 = ic50.append({'Kd': 7.73/on_rate, 'non_co_inact_ic50': "Excceed Iteration"}, ignore_index=True)
            print('Excceed Iteration')
            AMG510_IC50 = "Excceed Iteration"
            break
    return AMG510_IC50
    


# In[ ]:


import pandas as pd
import numpy as np

# Define the ranges for GDP and GTP concentrations
min_value_gdp = 1e-9  # Minimum GDP concentration
max_value_gdp = 1e-5  # Maximum GDP concentration
min_value_gtp = 1e-8  # Minimum GTP concentration
max_value_gtp = 1e-3  # Maximum GTP concentration
num_points = 20       # Number of points for each concentration

# Generate log-spaced values for GDP and GTP concentrations
testing_gdp_conc = np.logspace(np.log10(min_value_gdp), np.log10(max_value_gdp), num_points)
testing_gtp_conc = np.logspace(np.log10(min_value_gtp), np.log10(max_value_gtp), num_points)

# Create all combinations of GDP and GTP concentrations
gdp_gtp_combinations = [(gdp, gtp) for gdp in testing_gdp_conc for gtp in testing_gtp_conc]

# Create a DataFrame with the combinations
df = pd.DataFrame(gdp_gtp_combinations, columns=["GDP Concentration", "GTP Concentration"])

# Add a third column for simulation results, initialized as NaN
df["Simulation Result"] = np.nan

# run the simulation for AMG 510
df["Simulation Result"] = df.apply(lambda row: search_amg_ic50(row["GDP Concentration"], row["GTP Concentration"]), axis=1)

# Save the updated DataFrame to a CSV file
output_file = "amg510_gdp_gtp_combination.csv"
df.to_csv(output_file, index=False)


