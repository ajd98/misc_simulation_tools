# misc_simulation_tools
#------------------------------------------------------------------------------
# coloraccess.py
#
# Utility for checking the compatibility of color schemes with color vision 
# deficiencies.
# 
# How to use:
#
# This utility simulates a variety of color vision 
# deficiencies in order to aid the design of accessible 
# images and figures.  To use, specify an image file and  
# an output file (an image or .pdf). This utility will 
# tile the image, with different color vision deficiencies 
# simulated in each tile.  Note that this utility requires 
# matplotlib, NumPy, and SciPy.
#
# Run this script as "python coloraccess.py -i <my_image> -o <outputpath>"
# where <outputpath> has an image extension, such as .png, .jpg, or .pdf.
# coloraccess.py will deduce the output format from this extension.
# For more information, run `python coloraccess.py --help`.
#
#------------------------------------------------------------------------------
# cm.py
# 
# Calculate residue-residue contact frequency maps.
# See comments in cm.py for more information
#
#------------------------------------------------------------------------------
# fraction_native_contacts.py
#
# Calculate the fraction of native contacts for a protein, given a
# coordinate file in the fort.23 format, as well as a *.go.parameters
# file.
# Run "python fraction_native_contacts.py --help" for more information.
#
#------------------------------------------------------------------------------
# rmsdfinder.py 
# 
# Calculate the root-mean-square deviation of a set of coordinates relative to 
# a reference structure.
# Run "python rmsdfinder.py --help" for information on arguments and usage.
#
#------------------------------------------------------------------------------ 
# visualizer.py
#
# Visualize trajectories from coarse-grained protein models. Supports
# coordinate input from "west.h5" WESTPA data files, in addition to the fort.23
# format.  Requires the BBQ java program.  Run "python visualizer.py --help"
# for information on usage. 
