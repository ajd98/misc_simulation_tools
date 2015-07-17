#!/usr/bin/env python
# By Alex DeGrave, July 16 2015
'''
This script is designed to calculate the fraction of native contacts for a 
single frame of calbindin-AFF. Arguments should be supplied as follows:

python <myself.py> <fort.23> <calbindin.go.parameters> <output file name>

The script outputs to a text file.  Each line corresponds to a time point in 
the given fort.23 file.The value on the line is the fraction of native contacts. 
'''

import numpy 
import sys
from scipy.spatial.distance import cdist

def build_contact_list(contacts_file):
    '''
    Build a list of native contacts from a <protein>.go.parameters file.
    Return a list of contacts, where is contact is a tuple of the form:
    <resid 1> <resid 2> <sigma_native> 
    The variables ``resid 1`` and ``resid 2`` are one-index residue indices.
    '''
    # The first line of the contacts file is junk.
    contacts_file.readline()

    # The rest of the lines contain data, in the form:
    # <chain id for residue 1>   <residue 1 id>  <chain id for residue 2>  <residue 2 id> <native distance>  <well depth for go potential> 
    # We need on the two residue indices and the native distance.
    # Start by initializing the contact list
    contact_list = []
    for line in contacts_file:
        # Each field is buffered by 8 spaces. Split on these spaces to 
        # isolate each field
        split = line.split()
        # We care about items 2,4, and 5 (as explained above)
        contact = ( int(split[1]), int(split[3]) , float(split[4]) )
        
        # We do not want duplicates in the contact list.  For example, the go 
        # contact list supplied to this script will have a line both for the 
        # contact between residues 1 and 2, AND residues 2 and 1.  These are 
        # indeed the same contact, and it saves computational time if we 
        # exclude them.
        if not ( contact[1], contact[0], contact[2]) in contact_list:
            contact_list.append(contact)

    return contact_list

def line_to_coord_array(line):
    '''
    Converts a line from a fort.23 file to an N by 3 array of coordinates. 
    This method is useful for working on a coordinate file line by line, in 
    case loading the whole thing into memory at once would be problematic.
    Returns an array of coordinates.
    '''
    # Coordinates are buffered by spaces; we split on spaces.
    split = line.split()
    
    # The first number is the time point.  We do not care about it.
    data = split[1:]
   
    # Find the number of coordinates in the array
    number_of_coordinates = len(data)/3

    # Build and reshape the array.  The first dimension is the residue index 
    # (0 indexed!) and the second dimension is x,y,z coordinates
    array = numpy.array(data).reshape(number_of_coordinates, 3)
    
    return array

def calculate_fraction_of_native_contacts(coord_array, contact_list):
    '''
    Given a zero-indexed array of coordinates ``coord_array`` and a one-indexed
    array of coordinates ``contact_list``, return the fraction of native 
    contacts.  This method is specific to simulations with a Go-type potential,
    and most specifically the UIOWA-BD dynamics engine.  By definition, a 
    contact forms when the distance between pseudoatoms is less than or equal
    to 1.2 times the native sigma.
    Returns a float describing the fraction of native contacts.
    '''
    
    # First, find the distance matrix for the coordinate array.
    dm = cdist(coord_array, coord_array)
    
    # Initialize a counting variable.  This method will increment it every time
    # a native contact is found.
    native_contact_count = 0

    # Iterate over every native contact, and check if it exists.
    for contact in contact_list:
       
        # Convert residue indices to zero indexing scheme, for use with the 
        # coordinate array.
        zero_indexed_resid_1 = contact[0] - 1
        zero_indexed_resid_2 = contact[1] - 1

        # Extract the native_contact distance from the "contact" variable.
        native_distance = contact[2]

        # Check if the residues are in contact at the time point in question.  
        # By definition, residues are in contact if they are native contacts, 
        # and they are within 1.2 angstroms of their native distance.
        
        if dm[zero_indexed_resid_1, zero_indexed_resid_2] <= 1.2 * native_distance:
            native_contact_count += 1

    # Find the total number of contacts possible
    total_contacts = len(contact_list)
    
    # Now find the fraction of native contacts.
    Q = float(native_contact_count)/total_contacts

    return Q

def main():
    fort23_name        = sys.argv[1]
    contacts_file_name = sys.argv[2]
    output_file_name   = sys.argv[3]

    fort23 = open(fort23_name)
    contacts_file = open(contacts_file_name)
    output_file = open(output_file_name, 'w+')    

    # Extract the list of native contacts, with their associated distances.  
    # Each contact is of the form (resid 1, resid 2, distance, epsilon) .
    contact_list = build_contact_list(contacts_file)

    # Iterate over each line and calculate the fraction of native contacts.
    for line in fort23:
        
        # Extract the coordinates from the line
        coord_array = line_to_coord_array(line)
        Q = calculate_fraction_of_native_contacts(coord_array, contact_list)
        output_file.write('%f\n' % Q)

    output_file.close()
    fort23.close()
    contacts_file.close()  

if __name__=="__main__":
    main()
