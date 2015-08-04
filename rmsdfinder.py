#!/usr/bin/env python
#
# Usage: <myself>.py --reference folded.pdb --coords fort.23 
# See --help for more options.
# By Alex DeGrave: July 14, 2015.
# Updated by Alex DeGrave, August 4, 2015
#
# This utility calculates the RMSD of each timestep from a coordinate file
# (in the fort.23 format) specified by --coords, using the structure specified
# by --reference as the reference structure.  

# Import necessary utilities
import numpy
import sys
import argparse
import math


################################ Main Class ####################################

class RMSDFinder:
    def __init__(self):
        '''
        Initialize class. Parse arguments, open files, prepare the reference 
        structure, and create a mask for optionally skipping some atoms in the
        RMSD calculation.
        '''
        # Parse command line arguments.
        self._parse_args()

        # Open reference file.
        ref_path = self.args.reference_path
        reffile = open(ref_path,'r')

        # Prepare native structure, to which other states are compared.
        self.ref_structure = self.pdb_to_numpy_array(reffile)
        reffile.close()

        # Load the list of atom indices to skip.
        self._load_skip_list()

        self.ref_structure = self.ref_structure[self.skip_mask]
        self.ref_structure = self.ref_structure.reshape(self.ref_structure.shape[0]/3, 3)

    def _parse_args(self):
        '''
        Parse command line arguments. See ``--help`` for options.
        '''
        parser = argparse.ArgumentParser()
        parser.add_argument('--reference', dest='reference_path',
                            default='folded.pdb', type=str,
                            help='Use REFERENCE_PATH as the structure against which'
                                 'to calculate RMSD.  The file should be pdb '
                                 'format, with the same number and ordering of '
                                 'atoms as appears in the specified fort.23 file.')
        parser.add_argument('--coords', dest='coord_path',
                            default='fort.23', type=str,
                            help='Extract coordinates for calculating RMSD from '
                                 'the file specified by COORD_PATH.  The file '
                                 'should be in the fort.23 format (that is, a text'
                                 ' file where each line represents a stucture, and'
                                 ' entries are buffered by whitespace.  The first '
                                 'entry on each line should be the time, and each '
                                 'successive group of three entries represents the'
                                 ' x y z coordinates of each atom.  In other words'
                                 ', each line should have 3*N+1 whitespace-'
                                 'buffered entries, for a protein with N atoms. '
                                 'The atom ordering should correspond to that in '
                                 'the reference file.')
        parser.add_argument('--skip', dest='skip_list',
                            default=None, type=str,
                            help='Ignore atom indices specified by SKIP_LIST during'
                                 ' the RMSD calculation (both during least-squares '
                                 'alignment and root-mean-square distance '
                                 'calculation). SKIP_LIST should be a str that can '
                                 'be parsed as a python list or list comprehension.'
                                 ' The numpy and math libraries are made available '
                                 'during string evaluation. Atom indices should be ' 
                                 'zero-indexed')
        parser.add_argument('--low-memory', dest='low_memory_mode',
                            default=False, type=bool,
                            help='Use low memory mode if the user specifies '
                                 '``True`` for this option.  If enabled, low '
                                 'memory mode reads the coordinate file line-'
                                 'by-line, rather than all at once. For small '
                                 'files, this may be slower, as it no longer '
                                 'uses the well optimized ``loadtxt`` method '
                                 'of the numpy library.  For very large files,'
                                 ' this option may be preferable, as this '
                                 'utility will use far less memory.')
        self.args = parser.parse_args()

    def _load_skip_list(self):
        '''
        Load the list of coordinates to skip, and save a mask for coordinates 
        to ``self.skip_mask``.  Atoms with indices corresponding to elements 
        of this array set to one will be used in the RMSD calculation.
        '''
        temp_namespace = {'numpy': numpy,
                          'math' : math  }
        self.skip_mask = numpy.ones(self.ref_structure.shape, dtype=bool)
        if self.args.skip_list is not None:
            skip_list = eval(self.args.skip_list, temp_namespace)
            skip_list = numpy.array(skip_list)
            self.skip_mask[skip_list,:] = 0

    def pdb_to_numpy_array(self, pdbfile):
        '''
        Takes a pdb file (a python file object) and returns a numpy array.  
        Axis one is the atom number (in the case of UIOWA-BD's folded.pdb files,
        this is also the bead number).  Axis two is (x,y,z) coordinates.
        '''
        list_of_coordinates = []
        
        for line in pdbfile:
            if line[0:4] == 'ATOM':
                list_of_coordinates.append( line.split()[5:8] )
        return numpy.array(list_of_coordinates, dtype = numpy.dtype(float))
    
    def distance_squared(self, p1, p2):
        '''
        Returns the distance between two arrays (or tuples, or lists) of (x,y,z)
        coordinates, squared.
        '''
        return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
    
    def centroid(self, coordinate_array):
        '''
        Given a numpy array of x,y,z coordinates, this function returns the 
        centroid of these coordinates (points are weighted equally)
        '''
        x_sum = 0
        y_sum = 0
        z_sum = 0
        for point in coordinate_array:
            x_sum += point[0]
            y_sum += point[1]
            z_sum += point[2]
        # Calculates and returns (x,y,z) tuple for centroid.  
        return numpy.array((float(x_sum)/float(len(coordinate_array)), 
                         float(y_sum)/float(len(coordinate_array)), 
                         float(z_sum)/float(len(coordinate_array)) ))
    
    def RMSD(self, native_structure, new_structure): 
        '''
        Calculates the RMSD between two structures.  Structures are passed to this 
        function as numpy arrays.  Dimension one should be the atom or bead number 
        (these must correspond in each structure), and dimension two is x,y,z 
        coordinates.  Note that both structures must have the same dimensions.
        '''
    
        # Center both native_structure and new_structure on centroid.
        c = self.centroid(native_structure) 
        native_structure -= c
        c = self.centroid(new_structure)
        new_structure -= c
    
        # Use Kabsch algorithm to calculate optimal rotation matrix.
        # Calculate covariance matrix.
        covariance_matrix = numpy.dot(numpy.transpose(native_structure), 
                                      new_structure)
    
        # Singular Value Decomposition.
        V, S, Wt = numpy.linalg.svd(covariance_matrix)
        d = numpy.sign(numpy.linalg.det(numpy.dot(numpy.transpose(Wt),
                                                  numpy.transpose(V)
                                                  )
                                        )
                       )
    
        U = numpy.dot(numpy.transpose(Wt), 
                   numpy.dot(numpy.array(((1,0,0),
                                          (0,1,0),
                                          (0,0,d))), 
                             numpy.transpose(V)
                             )
                      )
    
        # Multiplying new_structure (n*3 matrix) by 3*3 optimal rotation matrix
        # ``U`` gives least_squares alignment.
        l_aligned = new_structure.dot(U)
    
        # Sum distances squared over each particle, and take the square root to 
        # return RMSD.
        square_sum = 0
        for i in range(len(l_aligned)):
            square_sum += self.distance_squared(l_aligned[i],native_structure[i])
        av = square_sum/len(l_aligned)
        rmsd = numpy.sqrt(av)
        return rmsd
    
    def run(self):
        '''
        Calculate root-mean-square deviation of structures in the coordinate
        file after alignment on the reference structure.  
        '''

        coord_path = self.args.coord_path
        # If low memory mode is enabled:
        if self.args.low_memory_mode is True:
            coordfile = open(coord_path,'r')
            for line in coordfile:

                # Create n*3 array for the final datapoint in the coordfile.
                line_array = numpy.array(line.split()[1:], 
                                         dtype=numpy.dtype(float)
                                         ).reshape(len(line.split()[1:])/3,3)
 
                line_array = line_array[self.skip_mask]
                line_array = line_array.reshape(line_array.shape[0]/3, 3)
                # Calculate rmsd between the reference structure and the  
                # stucture specified by the current line of the coordfile. 
                rmsd = self.RMSD(self.ref_structure, line_array)

                # Print the rmsd value to the screen.
                print(rmsd)
            coordfile.close()

        # If low memory mode is not enabled (default)
        else:
            coords = numpy.loadtxt(coord_path)
            coords = coords[:,1:]
            coords = coords.reshape(coords.shape[0], coords.shape[1]/3, 3)
            coords = coords[:,self.skip_mask]
            coords = coords.reshape(coords.shape[0], coords.shape[1]/3, 3)
            for timestep in coords:
                rmsd = self.RMSD(self.ref_structure, timestep)
                print(rmsd)

################################## Main Code ###################################

def main():    
    rmsdfinder = RMSDFinder()
    rmsdfinder.run()

if __name__ == '__main__':
    main()
