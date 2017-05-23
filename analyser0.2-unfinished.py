#!/usr/bin/env python
#
#  This script uses BoltzTraP output files and creates graphs for different quantites 
#  
#  Usage: ./analyser.py <seedname> <name_of_plot>
#  Additionally a temperature value needs to be defined 'temp_value= <number>' in this scipt
#  where  <number>  must be part of the temperature grid which BoltzTraP simulated.
# 
#  Available <name_of_plot> at the moment: seebeck, dos, sigma, chi, kappa0, all  
#
# version: 0.2
# last edit:  05/04/16 5 pm


import os
import sys
import numpy as np
from scipy.constants import e, m_e
import math
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

def main(argv = None):
    print "=================================================================="
    print "|                BoltzTraP Analyser for CASTEP                   |"
    print "|----------------------------------------------------------------|"
    print "|                           Usage:                               |"
    print "|           ./analyser.py <seedname> <name_of_plots>             |"
    print "|----------------------------------------------------------------|"
    print "|             Available <name_of_plots> arguments:               |"
    print "|            all, seebeck, dos, sigma, kappa0, chi               |"
    print "|----------------------------------------------------------------|"
    print "|                        Version: 0.2                            |"
    print "|                   Last edit: 06/04/2016                        |"
    print "==================================================================" 
    if argv is None:
        argv = sys.argv
    prefix = argv[1]
    temp_value = 300 # Temporary temperature control

    # Read intrans file, from here we will extract the number of temp steps
    intrans_file = prefix + '.intrans'
    intrans_file = open(intrans_file, 'r')
    intrans_data = intrans_file.readlines()
    intrans_file.close()

    # Calculate the number of temperature steps

    for line in intrans_data:
       if 'temperature grid' in line:
           temp_max = float(line.split()[0])
           temp_step = float(line.split()[1])
           num_temp_steps = int((temp_max/temp_step) - 1)
      
         
   

    # Open .outputtrans and gather info about Efermi, doping level and gap size
    
    with open(os.path.join(argv[1] + ".outputtrans"), 'r') as f:
            warning = False
            step = 0
            doping = []
            for line in f:
                if "WARNING" in line:
                    warning = True
                if "Efermi" in line:
                    efermi = (float(line.split()[5]))
                if line.startswith("Doping level number"):
                    doping.append(float(line.split()[6]))
                if line.startswith("Egap:"):
                    gap = float(line.split()[1])



    #---------------------------------------------------------------------------------------------------------------#    
    # Open the .trace file and append data to data_full array                                                       #
    # columns are:     Ef  |  Temp |   N   | DOS(Ef) |  S  |   sigma/t   |  R_H  |   kappa0  |     c     | chi      #
    # units            Ry  |   K   |  e/uc |  e/uc   | V/K | 1/(Ohm*m*s) | m^3/C | W/(m*K*s) | J/(mol*K) | m^3/mol  #
    # column number    0   |   1   |   2   |    3    |  4  |      5      |   6   |     7     |     8     |  9       #
    #---------------------------------------------------------------------------------------------------------------#
    data_full = []
    with open(os.path.join(argv[1] + '.trace'), 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    data_full.append([float(c) for c in line.split()])
    


    # Put each data into its own array for a given temperature    
    Ef_data = []
    N_data = []
    DOS_data = []
    S_data = []
    sigma_t_data = []
    R_H_data = []
    kappa0_data = []
    c_data = []
    chi_data = []
    
    for i in range(len(data_full)):
     if data_full [i][1] == temp_value:                  
       Ef_data.append(float((data_full[i][0])-efermi)*13.605698066) # Convert Ry to eV
       N_data.append(float(data_full[i][2]))
       DOS_data.append(float(data_full[i][3]))
       S_data.append(float(data_full[i][4]))
       sigma_t_data.append(float(data_full[i][5]))
       R_H_data.append(float(data_full[i][6]))
       kappa0_data.append(float(data_full[i][7]))
       c_data.append(float(data_full[i][8]))
       chi_data.append(float(data_full[i][9]))


    # Create new lists where we would store normalised data
    Ef_data_n = [None]*len(Ef_data)    # This is not needed but will be left for consistency 
    N_data_n = [None]*len(N_data)      # This is not needed but will be left for consistency
    DOS_data_n = [None]*len(DOS_data)
    S_data_n = [None]*len(S_data)
    sigma_t_data_n = [None]*len(sigma_t_data)
    R_H_data_n = [None]*len(R_H_data)
    kappa0_data_n = [None]*len(kappa0_data)
    c_data_n = [None]*len(c_data)
    chi_data_n = [None]*len(chi_data)
    
    #=========================================================================================#
    # Create a class which could be used to normalise the data                                #
    # First of all, the magnitude of the maximum value in a given data set will be checked    #
    # If positive -> use this value                                                           #
    # If negative -> decrease by 1                                                            #
    # Secondly, each data set would be divided by 1e+/- magnitude                            #
    # Next, lists defined just above will be populated with normalised values                 #
    # Lastly, a proper y axis label will be generated for each set                            #
    #=========================================================================================# 
    class Yvalue:
     def magnitude_of_y(self,data,data_n):
        self.data = data
        self.mag = int(math.log10(max(data)))
        # Determine if we need 1e+ or 1e-
        if self.mag >= 0:
           self.initial = '1e+'
        elif self.mag < 0:
           self.initial = '1e'
           self.mag = self.mag - 1
        # Populate new lists with normalised values  
        for i in range(len(data)):
            data_n [i] = float(data[i])/float('%s' % self.initial +'%s' % self.mag)
        #-----------------------------   Create y axis labels   ---------------------------# 
        # Density of States
        if self.data == DOS_data:
           if self.mag <> 0: 
              self.yval_string = 'DOS (10\\S%s' % self.mag + '\\0\\Na.u.)'
           else:
              self.yval_string = 'DOS (a.u.)'
        # Seebeck coefficient
        elif self.data == S_data:
           if self.mag <> 0: 
              self.yval_string = 'S (10\\S%s' % self.mag + '\\0\\NV/K)'
           else:
              self.yval_string = 'S (V/K)'
        # Electronic conductivity
        elif self.data == sigma_t_data:
           if self.mag <> 0: 
              self.yval_string = '\\xs/t\\0\\N(10\\S%s' % self.mag + ' \\0\\N/\\xW\\0ms)'
           else:
              self.yval_string = '\\xs/t\\0\\N(1/\\xW\\0ms)'
        # Hall coefficient
        elif self.data == R_H_data:
              pass
        # Electronic thermal conductivity
        elif self.data == kappa0_data:
           if self.mag <> 0: 
              self.yval_string = '\\xk\\s\\0\\el\\N/\\xt \\0\\N(10\\S%s' % self.mag + ' \\0\\NW/mKs)' 
           else:
              self.yval_string = '\\xk\\s\\0\\el\\N/\\xt \\0\\N(\\0\\NW/mKs)'             
        # Specific heat
        elif self.data == c_data:
           if self.mag <> 0: 
              self.yval_string = 'specific heat (10\\S%s' % self.mag + '\\0\\NJ/molK)'
           else:
              self.yval_string = 'specific heat (J/molK)'
        # Pauli magnetic susceptability
        elif self.data == chi_data:
           if self.mag <> 0: 
              self.yval_string = '\\xc\\0\\N(10\\S%s' % self.mag + ' \\0\\Nm\\S3\\N/mol)' 
           else:
              self.yval_string = '\\xc\\0\\N(\\0\\Nm\\S3\\N/mol)'
        return 

    

    # Create a function which will plot data from .trace
    # It needs 7 arguments: <title of the graph>, <legend string>, <xaxis label>, <yaxis label>, <x data>, <y data>, <name of output file>
    def trace_plots(title, legend, xaxis, yaxis, set_used_x, set_used_y, file_name):
           titlename = title
           legend = legend
           xaxis_label = xaxis
           yaxis_label = yaxis
           set_used_y = set_used_y
           set_used_x = set_used_x
           file_name = file_name
           f_text = '@    title "%s"' % titlename +'\n' + '@    s0 line linewidth 2.0' + '\n' + '@    s0 legend  "%s"\n' % legend 
           f_text += '@    s0 line color 2 \n'
           f_text += '@    xaxis  bar linewidth 2.0 \n'
           f_text += '@    xaxis  tick major linewidth 2.0 \n'
           f_text += '@    xaxis  tick minor linewidth 2.0 \n'
           f_text += '@    xaxis  label "%s" \n' % xaxis_label
           f_text += '@    xaxis  ticklabel char size 1.1800000 \n'
           f_text += '@    yaxis  bar linewidth 2.0 \n'
           f_text += '@    yaxis  tick major linewidth 2.0 \n'
           f_text += '@    yaxis  tick minor linewidth 2.0 \n'
           f_text += '@    yaxis  label "%s" \n' % yaxis_label
           f_text += '@    yaxis  ticklabel char size 1.1800000 \n'
           for i in range(len(S_data)):
               f_text += str(set_used_x[i]) + ' ' + str(set_used_y[i]) + '\n'
           foxy = open (file_name + '.agr', 'w')
           foxy.write(f_text)
           foxy.close()
           return f_text
    
    # Use class Yvalue to normalise all data sets
    # Generate strings for y axis labels 
    normaliser = Yvalue()
    normaliser.magnitude_of_y(DOS_data,DOS_data_n)
    yaxis_string_dos = normaliser.yval_string
    normaliser.magnitude_of_y(S_data,S_data_n)
    yaxis_string_S = normaliser.yval_string
    normaliser.magnitude_of_y(sigma_t_data,sigma_t_data_n)
    yaxis_string_sigma_t = normaliser.yval_string
    normaliser.magnitude_of_y(kappa0_data,kappa0_data_n)
    yaxis_string_kappa0 = normaliser.yval_string
    normaliser.magnitude_of_y(c_data,c_data_n)
    yaxis_string_c = normaliser.yval_string
    normaliser.magnitude_of_y(chi_data,chi_data_n)
    yaxis_string_chi = normaliser.yval_string
    
    
    xaxis_all_ef = '\\xe-\\xe\\0\\s0\\0\\N(eV)'
    for plot in argv:
     if plot == 'seebeck':
      trace_plots('Seebeck coefficient','Seebeck',xaxis_all_ef,yaxis_string, Ef_data, S_data, 'Seebeck')
     elif plot == 'dos':
      trace_plots('Density of States','DoS',xaxis_all_ef,'DoS', Ef_data, DOS_data, 'DoS')
     elif plot == 'sigma':
      trace_plots('Electronic conductivity','sigma',xaxis_all_ef,yaxis_string_sigma_t, Ef_data, sigma_t_data, 'sigma')
     elif plot == 'chi':
      trace_plots('Pauli magnetic susceptability','chi',xaxis_all_ef,'chi (m^3/mol)', Ef_data, chi_data, 'magnetic')
     elif plot == 'kappa0':
      trace_plots('Electronic thermal conductivity','kappa0',xaxis_all_ef,yaxis_string_kappa, Ef_data, kappa0_data, 'kappa0')
     elif plot == 'all':
      trace_plots('Seebeck coefficient','Seebeck',xaxis_all_ef,yaxis_string_S, Ef_data, S_data_n, 'Seebeck')
      trace_plots('Density of States','DoS',xaxis_all_ef, yaxis_string_dos, Ef_data, DOS_data_n, 'DoS')
      trace_plots('Electronic conductivity','\\xs',xaxis_all_ef,yaxis_string_sigma_t, Ef_data, sigma_t_data_n, 'elec_cond')
      trace_plots('Pauli magnetic susceptability','chi',xaxis_all_ef,yaxis_string_chi, Ef_data, chi_data_n, 'magnetic')
      trace_plots('Electronic thermal conductivity','\\xk\\s\\0\\el',xaxis_all_ef,yaxis_string_kappa0, Ef_data, kappa0_data_n, 'thermal_cond')
      trace_plots('Electronic specific heat','c',xaxis_all_ef,yaxis_string_c, Ef_data, c_data_n, 'heat_capacity')


    # -------------------------------Hall tensor analyser -------------------------------------#
    # Open halltens file
    with open(os.path.join(argv[1]+".halltens"), 'r') as f:
         data_hall = []
         for line in f:
             if not line.startswith("#"):
                 data_hall.append([float(c) for c in line.split()])
  
    # Create a 3x3x3 array with tensor values
    hall_tens = []
    for d in data_hall:
         if float(d[1]) == temp_value:        
            hall_tens_temp = [[[d[3], d[4], d[5]],
                          [d[6], d[7], d[8]],
                          [d[9], d[10], d[11]]],
                         [[d[12], d[13], d[14]],
                          [d[15], d[16], d[17]],
                          [d[18], d[19], d[20]]],
                         [[d[21], d[22], d[23]],
                          [d[24], d[25], d[26]],
                          [d[27], d[28], d[29]]]]
            hall_tens.append(hall_tens_temp)    



    #---------- Open .condtens file and create arrays for cond, seeback, kappa0 (x,x')  ------------#
    t_steps = set()
    m_steps = set()
    condtens_data = []
    cond_tensor = []
    seebeck_tensor = []
    kappa0_tensor = []
    with open(os.path.join(argv[1]+".condtens"), 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    t_steps.add(int(float(line.split()[1])))
                    m_steps.add(float(line.split()[0]))  
                    condtens_data.append([float(c) for c in line.split()]) 


    for d in condtens_data:
      if float(d[1]) == temp_value:   # gather info for a given temperature
        cond_tensor_temp = [[d[2], d[3], d[4]],
                         [d[5], d[6], d[7]],
                         [d[8], d[9], d[10]]]
        cond_tensor.append(cond_tensor_temp)
        seebeck_tensor_temp = [[d[11], d[12], d[13]],
                            [d[14], d[15], d[16]],
                            [d[17], d[18], d[19]]]
        seebeck_tensor.append(seebeck_tensor_temp)
        kappa0_tensor_temp = [[d[20], d[21], d[22]],
                          [d[23], d[24], d[25]],
                          [d[26], d[27], d[28]]]
        kappa0_tensor.append(kappa0_tensor_temp)  

        
# Calculate power factor - test - tensor dot product
#    results = []
#    print len(cond_tensor), len(m_steps)
#    for i in range(len(m_steps)):
#                full_tensor = np.dot(cond_tensor[i],
#                                     np.dot(seebeck_tensor[i],
#                                            seebeck_tensor[i]))
#                results.append(full_tensor)        
#    print full_tensor
#    test = 1   
#    if test == 1:                
#        seebeck = S_data
#        plt.plot(Ef_data, S_data, linewidth=3.0)
#        plt.legend(['S$_1$', 'S$_2$', 'S$_3$'])
#        plt.ylabel("Seebeck \n coefficient  ($\mu$V/K)", fontsize=30.0)
#        plt.xlabel("E-E$_f$ (eV)", fontsize=30)
#        plt.xticks(fontsize=25)
#        plt.yticks(fontsize=25)
#        plt.show()


#===================================================================================================================#
#============================================      DOPING PART      ================================================#
#===================================================================================================================#

    #---------------------------------------------------------------------------------------------------------------#    
    # Open the trace_fixdoping file and append data to data_full_doping array                                       #
    # columns are:    Temp |   N   | DOS(Ef) |  S  |   sigma/t   |  R_H  |   kappa0  |     c     |   chi   |  Ef    #
    # units            K   |  e/uc |  e/uc   | V/K | 1/(Ohm*m*s) | m^3/C | W/(m*K*s) | J/(mol*K) | m^3/mol |  Ry    #
    # column number    0   |   1   |    2    |  3  |      4      |   5   |     6     |     7     |    8    |  9     #
    #---------------------------------------------------------------------------------------------------------------#
    data_full_doping = []
    with open(os.path.join(argv[1] + '.trace_fixdoping'), 'r') as f:
            for line in f:
                if not line.startswith("#"):
                    data_full_doping.append([float(c) for c in line.split()])
    


   # Put each data into its own array for a given temperature     
    N_data_doping = []
    DOS_data_doping = []
    S_data_doping = []
    sigma_t_data_doping = []
    R_H_data_doping = []
    kappa0_data_doping = []
    c_data_doping = []
    chi_data_doping = []
    Ef_data_doping = []

    
    data_full_doping=filter(None, data_full_doping) #Remove empty entries from the BoltzTraP doping output
   
    for i in range(len(data_full_doping)):
     if data_full_doping [i][0] == temp_value:                   
       Ef_data_doping.append(float((data_full_doping[i][9])-efermi)*13.605698066) # Convert Ry to eV
       N_data_doping.append(float(data_full_doping[i][1]))
       DOS_data_doping.append(float(data_full_doping[i][2]))
       S_data_doping.append(float(data_full_doping[i][3]))
       sigma_t_data_doping.append(float(data_full_doping[i][4]))
       R_H_data_doping.append(float(data_full_doping[i][5]))
       kappa0_data_doping.append(float(data_full_doping[i][6]))
       c_data_doping.append(float(data_full_doping[i][7]))
       chi_data_doping.append(float(data_full_doping[i][8])) #Test

    # This part checks the maximum value in the S_data_doping list and if it is smaller than 1e-02 converts to microV/K  
    
    #IMPORTANT: if there is no data in the doping files this will display an error. Fix this
    if max(S_data_doping) < 1e-02:
       for i in range(len(S_data_doping)):
        S_data_doping [i] = float(S_data_doping[i])*1e+06
       yaxis_string_doping = 'S (\\xm\\0V/K)'
    else: 
       yaxis_string_doping = 'S (V/K)'
       pass
      

    # Create a function which will plot data from .trace_doping
    def trace_plots_dop(title, legend, xaxis, yaxis, set_used, file_name):
           titlename = title
           legend = legend
           xaxis_label = xaxis
           yaxis_label = yaxis
           set_used = set_used
           file_name = file_name
           f_text =  '@    title "%s"' % titlename +'\n' + '@    s0 line linewidth 2.0' + '\n' + '@    s0 legend  "%s"\n' % legend 
           f_text += '@    s0 line color 2 \n'
           f_text += '@    xaxis  bar linewidth 2.0 \n'
           f_text += '@    xaxis  tick major linewidth 2.0 \n'
           f_text += '@    xaxis  tick minor linewidth 2.0 \n'
           f_text += '@    xaxis  label "%s" \n' % xaxis_label
           f_text += '@    xaxis  ticklabel char size 1.1800000 \n'
           f_text += '@    yaxis  bar linewidth 2.0 \n'
           f_text += '@    yaxis  tick major linewidth 2.0 \n'
           f_text += '@    yaxis  tick minor linewidth 2.0 \n'
           f_text += '@    yaxis  label "%s" \n' % yaxis_label
           f_text += '@    yaxis  ticklabel char size 1.1800000 \n'
           for i in range(len(S_data_doping)):
               f_text += str(N_data_doping[i]) + ' ' + str(set_used[i]) + '\n'
           foxy = open (file_name + '.agr', 'w')
           foxy.write(f_text)
           foxy.close()
           return f_text
    trace_plots_dop('Seebeck coefficient doping','Seebeck','n',yaxis_string_doping, S_data_doping, 'Seebeck_doping')
      



if __name__ == "__main__":
    import sys
    sys.exit(main())
