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

class RMSDCalc:
    def __init__(self):
        pass

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
    

