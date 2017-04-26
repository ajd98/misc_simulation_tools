#!/usr/bin/env python
#
# Usage: <myself>.py --reference <reference.pdb> --coords <file1.pdb> [<file2.pdb> ... ]  
# See --help for more options.
# By Alex DeGrave: June 8, 2016.
#
# This utility calculates the RMSD of each timestep from coordinate file(s)
# (in the PDB format) specified by --coords, using the structure specified
# by --reference as the reference structure.  
#
# Make sure the atoms are in the same order in the reference and coordinate
# files!

# Import necessary utilities
import numpy
import argparse

################################ Main Class ####################################

class RMSDTool(object):
    def __init__(self):
        '''
        Initialize class. Parse arguments, open files, prepare the reference 
        structure.
        '''
        # Atom names for backbone residues; Python sets are much faster for "in" 
        #self.backboneres = {'N', 'C', "CA", 'O', 'OT1'} 
        self.backboneres = {'N', 'C', "CA"} 

        # Parse command line arguments.
        self._parse_args()

        # Prepare native structure, to which other states are compared.
        ref_path = self.args.reference_path
        self.ref_structure = self.pdb_to_numpy_array(ref_path)


    def _parse_args(self):
        '''
        Parse command line arguments. See ``--help`` for options.
        '''
        parser = argparse.ArgumentParser()
        parser.add_argument('--reference', dest='reference_path',
                            default='reference.pdb', type=str,
                            help='Use REFERENCE_PATH as the structure against '
                                 'which to calculate RMSD.  The file should be '
                                 'pdb format, with the same number and ordering'
                                 ' of atoms as appears in the specified '
                                 'coordinate files.')
        parser.add_argument('--coords', dest='coord_paths',
                            default=None, type=str, nargs='+',
                            help='Extract coordinates for calculating RMSD '
                                 'from the file(s) specified by COORD_PATHS. '
                                 'The file(s) should be in the PDB format.')
        parser.add_argument('--backbone', dest='backbone',
                            action='store_true', 
                            help='Include only backbone atoms (N, Ca, C, and O'
                                 ' in the RMSD calculation')
        self.args = parser.parse_args()

        if self.args.coord_paths is None:
            print("Please specify a coordinate file for finding RMSD.")
            exit(1)
        return


    def pdb_to_numpy_array(self, pdbpath):
        '''
        Takes a pdb file (a python file object) and returns a numpy array.  
        Axis one is the atom number. Axis two is (x,y,z) coordinates.
        '''
        list_of_coordinates = []
        
        pdbfile = open(pdbpath, 'r') 
        for line in pdbfile:
            if line[0:4] == 'ATOM':
                split = line.split()
                if self.args.backbone:
                    atomname = split[2] 
                    if atomname in self.backboneres:
                        list_of_coordinates.append(split[5:8])
                else:
                   list_of_coordinates.append(split[5:8])
        pdbfile.close()
        return numpy.array(list_of_coordinates, dtype=numpy.dtype(float))
    
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
        return coordinate_array.mean(axis=0)
    
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

        for pdbpath in self.args.coord_paths:
            coords = self.pdb_to_numpy_array(pdbpath) 
            rmsd = self.RMSD(self.ref_structure, coords)
            print(rmsd)

################################## Main Code ###################################
def main():    
    rmsdtool = RMSDTool()
    rmsdtool.run()

if __name__ == '__main__':
    main()
