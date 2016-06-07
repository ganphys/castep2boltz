#!/usr/bin/env python
#
# Fledgling CASTEP interface to the BOLTZTRAP program
#
# usage:
# castep2boltz <seedname> <optional arguments>
# available argument: so (for SOC) and down (for spin down calculations)
#
# Input files: <seedname>.castep
#              <seedname>.bands 
#
# Output files: <seedname>.energy or <seedname>.energyso
#               <seedname>.struct
#               <seedname>.intrans
#
# Required packages: spglib and ase
# If ase is missing try "sudo pip install --upgrade ase"
# If spglib is missing try "sudo pip install pyspglib"
#
# 
# edit 15.04.2016: Bug fixed in create .struct file section...
# Now the output writes the symm op matrices and not their transpose


import os
import sys
import numpy as np
from ase import Atoms
#from atoms import Atoms
try:
 from pyspglib import spglib
 has_spglib=True
except:
 has_spglib=False




def main(argv = None):
    print '========================================================='
    print '||             CASTEP 2 BoltzTraP Interface            ||'
    print '||                     version 1.0                     ||'
    print '||                    14 April 2016                    ||'
    print '||-----------------------------------------------------||'
    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
       # Avoid ugly errors
       print '|| Usage: castep2boltz <seedname> <optional arguments> ||'
       print '|| optional arguments: "so" (for SOC runs) ...         ||'
       print '||          and "down" (for spin down calculations)    ||'
       print '||-----------------------------------------------------||'
       print '|| UNSUCCESSFUL! READ Usage above                      ||'
       print '...'
       sys.exit()
       
    # Define <seedname>
    prefix = argv[1]

    # Help menu, it shows the message and stops the process
    help = ['h', '-h','--h', 'help', '-help', '--help']
    for i in help:
      if i in argv:
            print '|| Usage: castep2boltz <seedname> <optional arguments> ||'
            print '|| optional arguments: "so" (for SOC runs) ...         ||'
            print '||          and "down" (for spin down calculations)    ||'
            print '========================================================='
            sys.exit()
        

    # Check if an argument for SOC is given
    # so_on is defined here because it is used multiple times 
    if 'so' in argv or '-so' in argv:
     so_on = 'so'
    else  :
     so_on = None

    # Set a proper suffix for the .energy file
    if so_on == 'so':
      energy_file = prefix + '.energyso'    
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
          symmetry = True
          
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
              n_electrons = float(line.split()[3])
              n_electrons_down = float(line.split()[4])

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
                eigenenergies.append('{:3.10f}'.format(
                float(bands_data[eigenenergy_starting_line].split()[0])*2))
                eigenenergy_starting_line += 1  
        elif 'Spin component 2' in line:         
             eigenenergy_starting_line2 = index + 1
             for i in range(int(n_eigenvalues)):                        
                eigenenergies_spin2.append('{:3.10f}'.format(
                float(bands_data[eigenenergy_starting_line2].split()[0])*2))
                eigenenergy_starting_line2 += 1  
        elif 'Unit cell vectors' in line:  
              start = index + 1
              for j in range (3):
                   unit_cell.append(['{:3.10f}'.format(float(bands_data[start].split()[0])), 
                                     '{:3.10f}'.format(float(bands_data[start].split()[1])),
                                     '{:3.10f}'.format(float(bands_data[start].split()[2]))])
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
    

   # For debugging 
   # print 'positions',positions, 'atoms_symbols',a_symbols,'total num of ions', num_ions
   # print 'Unit cell', unit_cell
   # print 'unit_cell_from_castep' ,unit_cell_from_castep

  
    # Create an argument needed by spglib to generate symmetries
    all_atoms = Atoms(symbols = a_symbols,
                cell=unit_cell_from_castep,
                scaled_positions=positions,
                pbc=True)
        
    # Taken from the VASP interface, it checks if two symm ops are identical.
    # If they are, then don't add them to the list of operations.     
    def cmp_mat(mat1, mat2):
         absdiff = abs(mat1 - mat2)
         value = sum(sum(absdiff))
         if value > 1.0e-10:
            return False
         else:
            return True	
    
############################################################################################
# Debugging part
#    dataset = spglib.get_symmetry_dataset( all_atoms )
#    for i, (rot,trans) in enumerate( zip( dataset['rotations'], dataset['translations'] ) ):
#     print "  --------------- %4d ---------------" % (i+1)
#     print "  rotation:"
#     for x in rot:
#      print "     [%2d %2d %2d]" % (x[0], x[1], x[2])
#     print "  translation:"
#     print "     (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2])
#############################################################################################

    # Time for spglib
    # Check if spglib is present.
    if has_spglib and symmetry:
                # Use spglib and get rotational operations
                print '|| spglib found.                                       ||'
                rotations = spglib.get_symmetry_dataset(all_atoms, symprec=1e-5)['rotations']
                symm_ops_spglib = []
                
                # Append every value from rotations to symm_ops_spglib ...
                for i in range(len(rotations)):
                    new_rot = rotations[i]
                    newop = True
                    # ... but first use the VASP function to see 
                    # if the symmetry operation is unique ... 
                    # [I don't think this is needed but it doesn't hurt to have it]
                    for k in range(len(symm_ops_spglib)):
                        if cmp_mat(symm_ops_spglib[k],new_rot):
                            # If the operation is not unique don't add it
                            newop = False
                            break
                    # ... eventually when a new operation is found 
                    # append it symm_ops_spglib
                    if newop:
                        symm_ops_spglib.append(new_rot)
                if len(symm_ops_spglib) < 10:
                   print '||',len(symm_ops_spglib), 'symmetry operations generated by spglib.          ||'
                elif len(symm_ops_spglib) >= 10:
                   print '||',len(symm_ops_spglib), 'symmetry operations generated by spglib.         ||'

                # Check if the number of symmetry operations generated by spglib
                # is the same as the number of symmetry operations in CASTEP
                if len(symm_ops_spglib) != n_symm_ops:
                  print '||                                                     ||'
                  print '||                     WARNING:                        ||'
                  if n_symm_ops < 10:
                     print '||', n_symm_ops, 'symmetry operations found in %s.castep.     ||' % prefix
                  elif n_symm_ops >= 10: 
                     print '||', n_symm_ops, 'symmetry operations found in %s.castep.    ||' % prefix
                  print '||',len(symm_ops_spglib), 'symmetry operations generated by spglib.          ||'
                  print '|| CASTEP tolerance might be too loose.                ||'
                  print '|| Proceed with CAUTION.                               ||'
    # If spglib is not found pop a warning            
    else:
       if symmetry:
                print '||                                                     ||'
                print '||                     WARNING:                        ||'
                print '|| Symmetry operations ARE present...                  ||'
                print '|| ...but spglib was NOT found.                        ||'
                print '|| Proceed with CAUTION.                               ||'
                print '|| TRY installing spglib: sudo pip install pyspglib    ||'
       if not symmetry:
                print '||                                                     ||'
                print '||                     MESSAGE:                        ||'
                print '|| Symmetry operations NOT used in CASTEP calculation. ||'
                print '|| Continue NORMALLY.                                  ||'
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
    if 'down' in argv or '-down' in argv:
     # Check if there are spin down eigenvalues...
     if spin_components == 2:
       f.write(f_energy_down)
       f.close()
       print '||                                                     ||'
       print '|| Spin down energy file cooked and ready to go.       ||'
     # ...if not, use normal values.
     else:
       print '||                                                     ||'
       print '|| Spin polarisation not present, no spin down values. ||'
       print '|| Default settings used instead.                      ||' 
       f.write(f_energy)
       f.close()
    else:
       f.write(f_energy)
       f.close()
       

# IMPORTANT: If a kpoint is repeated BoltzTraP will give an error, 
# CASTEP BS calculations might have such point defined in the path
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
    if so_on == 'so':
       f_def += '10,\'' + prefix + '.energyso\',         \'old\',    \'formatted\',0\n'
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
            f_intrans += str(efermi_down) + ' ' + str(deltae) + ' ' + str(ecut) + ' ' + str(n_electrons_down) + '    # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons\n'
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
    f_intrans += 'HISTO                      #Scheme to obtain DOS. HISTO/TETRA: histogram/thetrahedron sampling\n'

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

    if so_on == 'so':
       print '||                                                     ||'
       print '|| SOC .energyso file cooked and ready to go.          ||'

    print '||                                                     ||'
    print '|| Done.                                               ||'
    print '========================================================='


if __name__ == "__main__":
    import sys
    sys.exit(main())

