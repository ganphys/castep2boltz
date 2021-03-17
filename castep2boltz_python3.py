#!/usr/bin/env python
#
# CASTEP interface to the BoltzTraP-v1.2.5 program
# Written by Genadi Naydenov <gan503@york.ac.uk> University of York (2016)
#
# USAGE: castep2boltz.py <seedname> <optional arguments>
# available argument: so (for SOC) and down (for spin down calculations)
#
# Input files: <seedname>.castep
#              <seedname>.bands 
#
# <seedname>.castep and <seedname>.bands should be obtained via
# DoS calculation. If a bandstructure calculation is used,
# BoltzTraP would give an error if there is a repetition of any
# of the k-points in the <seedname>.bands file.
#
# Output files: <seedname>.energy or <seedname>.energyso
#               <seedname>.struct
#               <seedname>.intrans
#
# Required packages: spglib and ase
# If ase is missing try "sudo pip install --upgrade ase"
# If spglib is missing try "sudo pip install pyspglib"

import os
import sys
#import numpy as np
try:
 from ase import Atoms
 ase_atoms=True
except:
 ase_atoms=False
try:
 from pyspglib import spglib
 has_spglib=True
except:
 has_spglib=False

symmetry_tol = 0.001  # CASTEP default value for symmetry tolerance is 0.001 ang


def main(argv = None):
    print( '=========================================================')
    print( '||             CASTEP 2 BoltzTraP Interface            ||')
    print( '||                     version 1.3                     ||')
    print( '||                     17 Mar 2021                     ||')
    print( '||-----------------------------------------------------||')
    if not ase_atoms:
       print( '|| ase library not found.                              ||')
       print( '|| TRY installing ase: sudo pip install --upgrade ase  ||')
       print( '|| UNSUCCESSFUL!                                       ||')
       print( '...')
       sys.exit()

    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
       # Avoid ugly errors
       print( '|| Usage: castep2boltz.py <seedname> <opt. arguments>  ||')
       print( '|| optional arguments: "so" (for SOC runs) ...         ||')
       print( '||          and "up/dn" (for spin polarised calc.)     ||')
       print( '||-----------------------------------------------------||')
       print( '|| Needed input:<seedname>.bands and <seedname>.castep ||')
       print( '|| UNSUCCESSFUL! READ Usage above                      ||')
       print( '...')
       sys.exit()
       
    # Define <seedname>
    prefix = argv[1]

    # Help menu, it shows the message and stops the process
    help = ['h', '-h','--h', 'help', '-help', '--help']
    for i in help:
      if i in argv:
            print( '|| Usage: castep2boltz.py <seedname> <opt. arguments>  ||')
            print( '|| optional arguments: "so" (for SOC runs) ...         ||')
            print( '||          and "up/dn" (for spin polarised calc.)     ||')
            print( '|| Needed input:<seedname>.bands and <seedname>.castep ||')
            print( '=========================================================')
            sys.exit()
        

    # Check if an argument for SOC is given
    # so_on is defined here because it is used multiple times 
    if 'so' in argv or '-so' in argv:
      so_on = True
    else:
      so_on = False

      
    # Check if an argument for spin down calculations is given 
    if 'down' in argv or '-down' in argv or 'dn' in argv or '-dn' in argv:
      spin_dn = True
    else:
      spin_dn = False

    # Check if an argument for spin down calculations is given 
    if 'up' in argv or '-up' in argv:
      spin_up = True
    else:
      spin_up = False

    # Set a proper suffix for the .energy file
    if so_on:
      energy_file = prefix + '.energyso' 
    elif spin_up:
      energy_file = prefix + '.energyup'    
    elif spin_dn:
      energy_file = prefix + '.energydn'    
    else:
      energy_file = prefix + '.energy'

    # Names of output files
    def_file = 'BoltzTraP.def'  
    intrans_file = prefix + '.intrans'
    struct_file = prefix + '.struct' 

#========================================================================================#
#                         Begin initial extraction of data from
#                         <seedname>.castep and <seedname>.bands 
#========================================================================================#
  
    # Open the .castep file and read it
    castep_file = prefix + '.castep'
    castep_file = open(castep_file, 'r')
    castep_data = castep_file.readlines()
    castep_file.close()
  
    # Check if there are any symmetry operations. 
    for index, line in enumerate(castep_data):
       if 'Number of symmetry operations' in line:
          n_symm_ops = int(float(line.split()[5]))
          if n_symm_ops != 1:
          	symmetry = True
          else:	
            symmetry = False
          
       elif 'There are no symmetry operations specified' in line:
          symmetry = False

    # Open the .bands file and read it
    bands_file = prefix + '.bands'
    bands_file = open(bands_file, 'r')
    bands_data = bands_file.readlines()
    bands_file.close()

    # Here we will store some of the data
    kpoints_frac_coordinates = []
    eigenenergies = []  # spin 1 (up)
    eigenenergies_spin2 = [] # spin 2 (down)
    unit_cell = [] # Crystal lattice
    
    # Extract number of kpoints, spin components, 
    # electrons, eigenvalues from the <seedname>.bands file
    # Extract values for Fermi energy, kpoints frac coordinates.
    for line in bands_data:
        if 'Number of k-points' in line:
            n_kpoints = float(line.split()[3])
        elif 'Number of spin components' in line:
            spin_components = float(line.split()[4])
        elif 'Number of electrons' in line:
           if spin_components == 1:
              n_electrons = float(line.split()[3])
              n_electrons_down = None
           elif spin_components == 2:
              n_electrons_up = float(line.split()[3])
              n_electrons_down = float(line.split()[4])
              n_electrons = n_electrons_up + n_electrons_down

        # Can you have different number of eigenvalues for
        # spin up and down channels? If yes, this part
        # should be rewritten. n_eigenvalues is used when
        # the .energy file is cooked and might give
        # wrong results if n_eigenvalues != n_eigenvalues_down
        elif 'Number of eigenvalues' in line:
            n_eigenvalues = float(line.split()[3]) 

        # This is present when there is 1 spin component
        elif 'Fermi energy' in line: 
            # castep output is in Hartree, this converts 
            # Fermi energy into Rydberg; Ry=2*Hartree
            efermi = float(line.split()[5])*2 
            # Set Fermi energy for spin down electrons to None
            # This might be used later for a quick check
            efermi_down = None
        # This is present when there are 2 spin components
        elif 'Fermi energies' in line: 
            efermi_up = float(line.split()[5])*2
            efermi_down = float(line.split()[6])*2 
            efermi=efermi_up  
        elif 'K-point' in line:
            kpoints_frac_coordinates.append([float(line.split()[2]),
                                             float(line.split()[3]),
                                             float(line.split()[4])]) 
                  
    # Get eigenenergies and unit cell from .bands file
    for index, line in enumerate(bands_data):
        if 'Spin component 1' in line:         
             eigenenergy_starting_line = index + 1
             for i in range(int(n_eigenvalues)):                        
                eigenenergies.append('{0:3.10f}'.format(
                float(bands_data[eigenenergy_starting_line].split()[0])*2))
                eigenenergy_starting_line += 1  
        elif 'Spin component 2' in line:         
             eigenenergy_starting_line2 = index + 1
             for i in range(int(n_eigenvalues)):                        
                eigenenergies_spin2.append('{0:3.10f}'.format(
                float(bands_data[eigenenergy_starting_line2].split()[0])*2))
                eigenenergy_starting_line2 += 1  
        elif 'Unit cell vectors' in line:  
              start = index + 1
              for j in range (3):
                   unit_cell.append(['{0:3.10f}'.format(float(bands_data[start].split()[0])), 
                                     '{0:3.10f}'.format(float(bands_data[start].split()[1])),
                                     '{0:3.10f}'.format(float(bands_data[start].split()[2]))])
                   start += 1
#========================================================================================#                  
#========================================================================================#



#========================================================================================#
#                                    SPGLIB SECTION                                      #
#----------------------------------------------------------------------------------------#
# spglib will generate symmetry operations even if CASTEP doesn't use symmetry_generate  #
#----------------------------------------------------------------------------------------#
#          If symmetry_generate IS used in CASTEP, spglib SYMMETRY OPERATIONS            #
#                            will be added to the .struct file.                          #
#                                                                                        #
#                If symmetry_generate is NOT used, then the IDENTITY MATRIX              #
#                             will be added to the .struct file.                         #
#----------------------------------------------------------------------------------------#
#  At the end of the section there is an if clause which checks if the number of CASTEP  #
#     symmetry operations is the same as the number of the ones generated by spglib.     #
#                      If they are different, a warning will pop up.                     # 
#========================================================================================#
   
    positions = []
    a_symbols = []
    unit_cell_from_castep = []
        
    # Find the total number of ions
    for line in castep_data:
        if 'Total number of ions' in line:
           num_ions = int(line.split()[7])
        elif 'Cell is a supercell' in line:
          sup_cell = True
          num_sup_cells = int(line.split()[5])
          sup_cell_msg = line
          #print num_sup_cells
    
    # Get unit cell from <seedname>.castep.
    for index, line in enumerate(castep_data):
         # Crystal lattice in (A) units. <seedname>.bands file also contains
         # this information. However, it uses Bohr units and this creates 
         # some problems when symmetry operations are generated with spglib. 
         if 'Real Lattice(A)' in line:  
                           start = index + 1
                           for j in range (3):
                             unit_cell_from_castep.append(
                                              [float(castep_data[start].split()[0]), 
                                               float(castep_data[start].split()[1]),
                                               float(castep_data[start].split()[2])])
                             start += 1  
                           break # avoid double counting  

    # Get atomic positions and symbols from .castep
    # Use the total number of ions and append the position of every atom to positions = []        
    for index, line in enumerate(castep_data):
         if 'Cell Contents' in line:
            for i in range(0, num_ions):
                positions.append([float(castep_data[index+10+i].split()[3]), 
                                  float(castep_data[index+10+i].split()[4]), 
                                  float(castep_data[index+10+i].split()[5])])
                a_symbols.append(str(castep_data[index+10+i].split()[1]))
            break # avoid double counting   
                  # if multiple castep runs are present in one .castep file 
    

    # Create an argument needed by spglib to generate symmetries
    all_atoms = Atoms(symbols = a_symbols,
                cell=unit_cell_from_castep,
                scaled_positions=positions,
                pbc=True)
    # Translational symm ops are not considered by boltztrap
    # Add only unique rotational symm ops
    def compare_elements(element1, element2):
    	absdiff = abs(element1 - element2)
    	value = sum(sum(absdiff))
    	if value > 1.0e-10:
        	return False
    	else:
        	return True

    # Check if spglib is present.
    if has_spglib and symmetry:
                # Use spglib and get rotational operations
                print( '|| spglib found.                                       ||')
		# used for debugging
		#cell_symm = spglib.get_symmetry(all_atoms, symprec=symmetry_tol)
		#translations = spglib.get_symmetry_dataset(all_atoms, symprec=symmetry_tol)['translations']
                rotations = spglib.get_symmetry_dataset(all_atoms, symprec=symmetry_tol)['rotations']
		
                symm_ops_spglib = []

                for i in range(len(rotations)):
                    new_rot = rotations[i]
    
                    newop=True
                    for i in range(len(symm_ops_spglib)):
                        if compare_elements(symm_ops_spglib[i],new_rot):
                            newop=False
                            break
                    if newop:
                        symm_ops_spglib.append(new_rot)
                
                
                def max_length(msg):
                  max_length1 = abs(57-len(msg))
                  return max_length1

                # print number of symmetry operations generated by spglib.
                mess = '|| %s ' %len(symm_ops_spglib) + 'symmetry operations generated by spglib. '
                
                print( '%s' %mess + '{0:>{width}}'.format("||", width=max_length(mess)))

                # Check if the number of symmetry operations generated by spglib
                # is the same as the number of symmetry operations in CASTEP
                if len(symm_ops_spglib) != n_symm_ops:
                  print( '||                                                     ||')
                  print( '||                     WARNING:                        ||')
                  err_mess = '|| %s ' %n_symm_ops + 'symmetry operations found in %s.castep. ' % prefix                  
                  print( '%s' %err_mess + '{0:>{width}}'.format("||", width=max_length(err_mess)))
                  err_mess = '|| %s ' %len(symm_ops_spglib) + 'symmetry operations generated by spglib. '
                  print( '%s' %err_mess + '{0:>{width}}'.format("||", width=max_length(err_mess)))
                  if sup_cell:
                      if len(symm_ops_spglib)*num_sup_cells == n_symm_ops: 
                        sup_cell_msg = '|| Cell is a supercell containing %s ' %num_sup_cells + 'primitive cells '
                        print( '||                                                     ||')
                        print( '%s' %sup_cell_msg + '{0:>{width}}'.format("||", width=max_length(sup_cell_msg)))
                        print( '|| and translational sym. ops. are not needed.         ||')
                        print( '|| Everything should be fine.                          ||')
                      else:
                        sup_cell_msg = '|| Cell is a supercell containing %s ' %num_sup_cells + 'primitive cells. '
                        print( '||                                                     ||')
                        print( '%s' %sup_cell_msg + '{0:>{width}}'.format("||", width=max_length(sup_cell_msg)))
                        print( '|| However, the number of sym. ops. is unexpected.     ||')
                        print( '|| Check symmetry operations tolerance in CASTEP.      ||')
                        print( '|| Proceed with CAUTION.                               ||')

                  
                  else:
                      print( '|| Check symmetry operations tolerance in CASTEP.      ||')
                      print( '|| Proceed with CAUTION.                               ||')
    # If spglib is not found pop a warning            
    else:
       if symmetry:
                print( '||                                                     ||')
                print( '||                     WARNING:                        ||')
                print( '|| Symmetry operations ARE present...                  ||')
                print( '|| ...but spglib was NOT found.                        ||')
                print( '|| Proceed with CAUTION.                               ||')
                print( '|| TRY installing spglib: sudo pip install pyspglib    ||')
       if not symmetry:
                print( '||                                                     ||')
                print( '||                     MESSAGE:                        ||')
                print( '|| Symmetry operations NOT used in CASTEP calculation. ||')
                print( '|| Continue NORMALLY.                                  ||')
#========================================================================================#    
#========================================================================================#
          
    

#========================================================================================#
#
#                          Create the <seedname>.energy file
#
#========================================================================================#         

#---------------------------------- Spin 1/ Spin up PART --------------------------------# 

    # First line is a comment, in this case it is the name of the structure   
    f_energy = prefix + '\n'

    # Number of k-points
    f_energy += str(int(n_kpoints)) + '\n'

    # K-point coordinates and number of eigenvalues
    for ik in range(int(n_kpoints)):
        for j in range(3):
          f_energy += str(kpoints_frac_coordinates[ik][j]) + ' '
        f_energy += str(int(n_eigenvalues)) + '\n'     
        start_pos_for_next_k_point = ik*int(n_eigenvalues)

        # Eigenvalues for a given k-point

        # if you want to exclude bands change range to
        # range(int(n_excluded:n_eigenvalues)) and 
        # add this as offset to start_pos ------> this would only exclude the bottom bands
        for ib in range(int(n_eigenvalues)): 
            f_energy += str(eigenenergies[start_pos_for_next_k_point]) + '\n'
            start_pos_for_next_k_point += 1 
 

#------------------------------- Spin 2/ Spin down PART ---------------------------------#
    # Similar to Spin 1 PART  
    if spin_components == 2:
     f_energy_down = prefix + ' down ' + '\n'
     f_energy_down += str(int(n_kpoints)) + '\n'
     for ik in range(int(n_kpoints)):
        for j in range(3):
          f_energy_down += str(kpoints_frac_coordinates[ik][j]) + ' ' 
        f_energy_down += str(int(n_eigenvalues)) + '\n'     
        start_pos_for_next_k_point = ik*int(n_eigenvalues) 

        # if you want to exclude bands change range to
        # range(int(n_excluded:n_eigenvalues)) and 
        # add this as offset to start_pos ------> this would only exclude the bottom bands
        for ib in range(int(n_eigenvalues)):
            f_energy_down += str(eigenenergies_spin2[start_pos_for_next_k_point]) + '\n'
            start_pos_for_next_k_point += 1 
    else:
     pass

     
    f = open(energy_file, 'w')
    
    # Choose whether spin up or spin down energies
    # will be written to the .energy file
    if 'down' in argv or '-down' in argv or 'dn' in argv or '-dn' in argv:
     # Check if there are spin down eigenvalues...
     if spin_components == 2:
       f.write(f_energy_down)
       f.close()
       print( '||                                                     ||')
       print( '|| Spin down energy file cooked and ready to go.       ||')
     # ...if not, use normal values.
     else:
       print( '||                                                     ||')
       print( '|| Spin polarisation not present, no spin down values. ||')
       print( '|| Default settings used instead.                      ||' )
       f.write(f_energy)
       f.close()
    else:
       f.write(f_energy)
       f.close()
       

# IMPORTANT: If a kpoint is repeated BoltzTraP will give an error, 
# CASTEP BS calculations might have a such point defined in the path
# and it needs to be removed.
#========================================================================================#
#========================================================================================#



#========================================================================================#
#
#                          Create the BoltzTraP.def file
#
#========================================================================================#

    f_def = '5, \'' + prefix + '.intrans\',      \'old\',    \'formatted\',0\n'
    f_def += '6,\'' + prefix + '.outputtrans\',      \'unknown\',    \'formatted\',0\n'
    f_def += '20,\'' + prefix + '.struct\',         \'old\',    \'formatted\',0\n'
    if so_on:
       f_def += '10,\'' + prefix + '.energyso\',         \'old\',    \'formatted\',0\n'
    elif spin_up:
       f_def += '10,\'' + prefix + '.energyup\',         \'old\',    \'formatted\',0\n'
    elif spin_dn:
       f_def += '10,\'' + prefix + '.energydn\',         \'old\',    \'formatted\',0\n'
    else:
       f_def += '10,\'' + prefix + '.energy\',         \'old\',    \'formatted\',0\n'
    f_def += '48,\'' + prefix + '.engre\',         \'unknown\',    \'unformatted\',0\n'
    f_def += '49,\'' + prefix + '.transdos\',        \'unknown\',    \'formatted\',0\n'
    f_def += '50,\'' + prefix + '.sigxx\',        \'unknown\',    \'formatted\',0\n'
    f_def += '51,\'' + prefix + '.sigxxx\',        \'unknown\',    \'formatted\',0\n'
    f_def += '21,\'' + prefix + '.trace\',           \'unknown\',    \'formatted\',0\n'
    f_def += '22,\'' + prefix + '.condtens\',           \'unknown\',    \'formatted\',0\n'
    f_def += '24,\'' + prefix + '.halltens\',           \'unknown\',    \'formatted\',0\n'
    f_def += '30,\'' + prefix + '_BZ.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '31,\'' + prefix + '_fermi.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '32,\'' + prefix + '_sigxx.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '33,\'' + prefix + '_sigyy.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '34,\'' + prefix + '_sigzz.dx\',           \'unknown\',    \'formatted\',0\n'
    f_def += '35,\'' + prefix + '_band.dat\',           \'unknown\',    \'formatted\',0\n'
    f_def += '36,\'' + prefix + '_band.gpl\',           \'unknown\',    \'formatted\',0\n'
    f_def += '37,\'' + prefix + '_deriv.dat\',           \'unknown\',    \'formatted\',0\n'
    f_def += '38,\'' + prefix + '_mass.dat\',           \'unknown\',    \'formatted\',0\n'

    f = open(def_file, 'w')
    f.write(f_def)
    f.close()
#========================================================================================#
#========================================================================================#



#========================================================================================#
#
#                          Create the <seedname>.intrans file
#
#========================================================================================#

    deltae = 0.0005
    ecut = 0.4
    lpfac = 5
    efcut = 0.15
    tmax = 800.0
    deltat = 50.0
    ecut2 = -1.0

    # efermi is only for spin up at the moment, add an if statement if spin down is present
    f_intrans = 'GENE                      # Format of DOS\n'
    f_intrans += '0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap\n'
    


    if 'down' in argv or '-down' in argv:
        if spin_components == 2:
            f_intrans += str(efermi_down) + ' ' + str(deltae) + ' ' + str(ecut) + ' ' + str(n_electrons) + '    # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n'
        else:
            # This could be regarded as an error or a typo. 
            # This part is also present in the .energy file section and should pop a warning if triggered.
            f_intrans += str(efermi) + ' ' + str(deltae) + ' ' + str(ecut) + ' ' + str(n_electrons) + '    # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n'

    else:
    	f_intrans += str(efermi) + ' ' + str(deltae) + ' ' + str(ecut) + ' ' + str(n_electrons) + '    # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n'
    f_intrans += 'CALC                      # CALC (calculate expansion coeff), NOCALC read from file\n'
    f_intrans += str(lpfac) + '                         # lpfac, number of latt-points per k-point\n'
    f_intrans += 'BOLTZ                     # run mode (only BOLTZ is supported)\n'
    f_intrans += str(efcut) + '                      # (efcut) energy range of chemical potential\n'
    f_intrans += str(tmax) + ' ' + str(deltat) + '                # Tmax, temperature grid\n'
    f_intrans += str(ecut2) + '                      # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)\n'
    f_intrans += 'HISTO                     # Scheme to obtain DOS. HISTO/TETRA: histogram/thetrahedron sampling\n'
    f_intrans += '#0 0 0 0 0                # tau-model. Not documented\n'
    f_intrans += '#14                       # number of fixed dopings\n'
    f_intrans += '#1E20 5E20 1E21 2E21 3E21 5E21 1E22 5E22 -1E20 -5E20 -1E21 -2.5E21 -5E21 -1E22 # example of fixed doping levels in cm^3\n'

    f = open(intrans_file, 'w')
    f.write(f_intrans)
    f.close()
#========================================================================================#
#========================================================================================#
 


#========================================================================================#
#
#                          Create the <seedname>.struct file
#
#========================================================================================#

    # First line is a comment, name of the structure in this case.
    f_struct = prefix + '\n'
    
    # Crystal lattice from <seedname>.bands file
    for i in range(3):
     for j in range(3):
        f_struct += str(unit_cell[i][j]) + ' ' 
     f_struct +='\n'

    # Check if symmetry operations and spglib are present
    if symmetry is True and has_spglib:
         f_struct += str(len(symm_ops_spglib)) + '\n'
         for i in range(int(len(symm_ops_spglib))):
          for j in range(3):
           for k in range(3):
            f_struct += str(symm_ops_spglib[i][j][k]) + ' '
          f_struct += '\n'
        
    
    # If there are no symmetry operations, append an identity matrix
    # to the .struct file. BoltzTraP doesn't run otherwise.
    elif symmetry is False or not has_spglib:
            f_struct += '1' + '\n'
            f_struct += '1 0 0 0 1 0 0 0 1'
    
      
    f = open(struct_file, 'w')
    f.write(f_struct)
    f.close() 
#========================================================================================#
#========================================================================================#

    if so_on:
       print( '||                                                     ||')
       print( '|| SOC .energyso file cooked and ready to go.          ||')

    print( '||                                                     ||')
    print( '|| Done.                                               ||')
    print( '=========================================================')


if __name__ == "__main__":
    import sys
    sys.exit(main())
