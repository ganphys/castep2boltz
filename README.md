# CASTEP2BoltzTraP interface

This directory contains CASTEP interface files to the BoltzTraP package for the solution of the semi-classical Boltzmann Transport Equation for electrons. The interface allows the user to convert the output files generated from a CASTEP run into a format which can be used by BoltzTraP.
BoltzTraP is NOT part of CASTEP and the CASTEP authors are not responsible for it.
You can find further information on BoltzTraP, including features and downloads, at http://www.imc.tuwien.ac.at/division_theoretical_chemistry/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap/EN/

## Getting Started

You can obtained the interface either as part of the CASTEP distribution, or download the most recent version from https://github.com/ganphys/castep2boltz 

The current version of the interface (castep2boltz.py v.1.1) is compatible with BoltzTraP v.1.2.5

### Prerequisites

The interface is a python script which uses python2.7 and requires the following libraries: numpy, ase and spglib  

Before running the script the user needs to obtain <seedname>.castep and <seedname>.bands files from a CASTEP DoS calculation.

### Installing

There is no need to install the script. If any of the prerequisite python libraries are missing, please try to install them through the terminal.

numpy

```
sudo apt-get install python-numpy
```
ase

```
sudo pip install --upgrade ase
```

spglib

```
sudo pip install pyspglib
```

## Running castep2boltz.py

You can run the interface by simply calling the script:

```
castep2boltz.py <seedname> <optional arguments>
```

The optional arguments can be: so (for SOC calculations) or down (for spin down calculations)
Output files for spin-up calculations are prepared without using any optional arguments. 

If executed successfully, the script will create 4 output files which can be used later on in BoltzTraP:

```
<seedname>.energy or <seedname>.energyso
<seedname>.struct
<seedname>.intrans
BoltzTraP.def
```

### BolztTraP notes

Once BoltzTraP is installed (please refer to BoltzTraP website for instructions), it needs to be called using the x_trans script:

```
x_trans BolztTraP -f <seedname>
```

The <seedname.intrans> file contains the parameters for the BolzTraP run and can be modified by the user. Temperature range and step size can be changed from there. If the user wants to add doping levels, please uncomment and edit accordingly the last 3 lines of <seedname.intrans>. The tau model line needs to be uncommented but can be left unchanged.
 
## Author

* **Genadi Naydenov** - (https://github.com/ganphys/castep2boltz)
