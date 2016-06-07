# castep2boltz.py
# CASTEP interface to the BOLTZTRAP program
#
# usage:
# castep2boltz [seedname] [optional arguments]
# available argument: so (for SOC) and down (for spin down calculations)
#
# Input files: [seedname].castep
#              [seedname].bands 
#
# Output files: [seedname].energy or [seedname].energyso
#               [seedname].struct
#               [seedname].intrans
#
# Required packages: spglib and ase
# If ase is missing try "sudo pip install --upgrade ase"
# If spglib is missing try "sudo pip install pyspglib"
#
# 
# edit 15.04.2016: Bug fixed in create .struct file section...
# Now the output writes correctly the symm op matrices rather than their transpose
