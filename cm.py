#!/usr/bin/env python
import numpy
import h5py
from scipy.spatial.distance import cdist

'''
This script calculates the average residue-residue contact map for the 
apo N': holo N calbindin AFF construct, using data from five simulations.
In this version of the script, state definitions are based on an assignments
file (for the N and N' states), and the isocommittor bins output by
w_comprobs.py (for the transition state). 
'''

def load_transition_state_bins():
    '''
    Loads the bins for the transition state into a properly formatted array.
    The bin boundaries are loaded from the file specified in this function.
    '''
    comprobH5 = h5py.File('../comprob.h5','r') 
    bin_idx_arr = comprobH5['isocommittor_bins'] 
    bin_idx_list = list(bin_idx_arr)
    return bin_idx_list

def load_contact_list(go_parameters_file_path):
    '''
    Loads a list of contacts from a go.parameters file. Format the list as: 

            [[particle_idx_1, particle_idx_2, distance_in_native_structure], 
             [particle_idx_1, particle_idx_2, distance_in_native_structure], 
             ...
             [particle_idx_1, particle_idx_2, distance_in_native_structure]] 
    '''
    contacts_file = open(go_parameters_file_path,'r')
    # The first line of the contacts file is junk.
    contacts_file.readline()
    # The rest of the lines contain data, in the form:
    # <chain id for residue 1>   <residue 1 id>  <chain id for residue 2>  <residue 2 id> <native distance>  <well depth for go potential> 
    # We need on the two residue indices and the native distance.
    # Start by initializing the contact list
    contact_list = []
    for line in contacts_file:
        # Each field is buffered by 8 spaces. We split on these spaces to 
        # isolate each field
        split = line.split()
        # We care about items 2,4, and 5 (as explained above)
        contact = ( int(split[1]), int(split[3]) , float(split[4]) )
        
        # We do not want duplicates in the contact list.  For example, the go 
        # contact list supplied to this script will have a line both for the 
        # contact between residues 1 and 2, AND residues 2 and 1.  These are 
        # indeed the same contact, and it saves computational time if we exclude 
        # them.
        if not ( contact[1], contact[0], contact[2]) in contact_list:
            contact_list.append(contact)
    return contact_list

def load_distance_reference_array_and_mask(contact_list):
    '''
    Creates a numpy array representing a square distance reference matrix from
    a list of contacts.  This functino accepts a contact list formatted as:

            [[particle_idx_1, particle_idx_2, distance_in_native_structure], 
             [particle_idx_1, particle_idx_2, distance_in_native_structure], 
             ...
             [particle_idx_1, particle_idx_2, distance_in_native_structure]] 

    Return a square numpy array of dimension 113*113 (the number of particles 
    in the calbindin-AFF system, which is currently hardcoded).

    Elements C_i,j that are nonzero in the corresponding upper triangular
    matrix represent the distance between pseudoatoms i and j multiplied by
    a margin of 1.2, if they are involved in a native contact. 

    This function sets elements C_i,j that are nonzero in the corresponding 
    lower triangular matrix to 6.5, the global cutoff for nonnative contacts,
    if they are not involved in a native contact (and therefore may be involved
    in a nonnative contact.

    Finally, this function sets elements corresponding to pairs of psuedoatoms
    that may be involved in local interactions to zero (ie, i->i+1, i->i+2, 
    i->i+3, and i->i).
    '''
    reference_array = numpy.zeros((113,113),dtype=numpy.float)
    mask = numpy.zeros((113,113),dtype=numpy.bool_)
    # Iterate over the contacts in the variable contact_list
    for contact in contact_list:
        # The part of the array with the lower residue index first is reserved
        # for native contacts. For example, reference_array[0,1] must refer to a
        # native contact, while reference_array[1,0] must refer to a nonnative 
        # contact. 
        indices = sorted([contact[0],contact[1]])
        # Convert from 1 indexing in the go.parameters file to zero indexing
        # for numpy arrays. Also, apply the 1.2 error margin.
        reference_array[indices[0]-1,indices[1]-1] = 1.2*contact[2]
        # For the index i,j describing a native contact, set the mask[j,i] to 1.
        # Later, we can call contact_map[mask] = 0 in order to make sure that
        # native contacts are not included in the nonnative portion of the map. 
        mask[indices[1]-1,indices[0]-1] = 1

    # Certain elements of the mask need to be set to one in order to exclude 
    # interactions accounted for by bonds (i,i+1) and angles (i,i+2). Dihedrals
    # (i,i+3) are included, because even though the are taken care of by the 
    # dihedral terms (a local interaction), it still may be useful to draw 
    # parallels with things such as alpha-helicity.
    for i in range(mask.shape[0]):
        if i+1 < mask.shape[0]:
            mask[i,i+1] = 1
            mask[i+1,i] = 1
        if i+2 < mask.shape[0]:
            mask[i,i+2] = 1
            mask[i+2,i] = 1

    # Set the cutoff distance for nonnative contacts to 6.5 angstroms.
    upper_triangular_ones = numpy.tril(numpy.ones(reference_array.shape,dtype=numpy.bool_))
    upper_triangular_ones_with_zero_diagonal = upper_triangular_ones - numpy.identity(reference_array.shape[0],dtype=numpy.bool_)
    reference_array[upper_triangular_ones_with_zero_diagonal] = 6.5
    reference_array[mask] = 0
    return reference_array

def gen_contact_map(coordinate_array):
    '''Generates a contact map, a binary numpy array representing a square
    matrix, where elements C_i,j that are nonzero in the corresponding upper 
    triangular matrix indicate a native contact currently formed between 
    residues i and j if ``1`` and no current contact if ``0``, and elements
    C_i,j that are nonzero in the corresponding lower triangular matrix
    indicate a nonnative contact currently formed between residues i and j if 
    ``1`` and no current nonnative contact if ``0``. '''
    distance_array = cdist(coordinate_array,coordinate_array)
    contact_map_array = distance_array < distance_reference_array 
    return contact_map_array

def calculate_average_contact_map(westh5, assignh5, binlist=[], 
                                  statelist=None, fi=1, li=None):
    '''
    Calculate the average contact map for structures in one of the specified 
    bins.  

    Arguments:
      westh5: a WESTPA main HDF5 file containing coordinates stored in
        "auxdata/coords"
      assignh5: a WESTPA assignments file corresponding to the given west.h5
        file.  
      binlist: a Python list containing the indices (integers) of bins to be
        included in the contact map calculation. 
      statelist: a Python list containing the indices (integers) of states to be
        included in the contact map calculation. This list will be converted to 
        bin indices and unioned with `binlist`. 
      fi: (integer) include weighted ensemble iterations beginning with `fi`.
        (Default: 1)
      li: (integer) include weighted ensemble iterations up to and including 
          `li` (Default: last iteration)

    Returns:
      (1) for a system with N particles, return an N*N numpy array, where 
        element i,j corresponds to the frequency of observing a contact between
        residues i and j, normalized by the total probability observed in the 
        selected bins or states. 
    '''
    # Convert states to bins.
    binset = set(binlist)
    if statelist is not None:
        state_map = numpy.array(assignh5['state_map'][...])
        for state_idx in statelist:
            # Check the index
            newbins = numpy.where(state_map==state_idx)[0]
            newbins = set(newbins)
            binset = binset.union(newbins)
    binlist = list(binset)

    total_weight = 0
    contact_map = None
    for iter_idx in xrange(fi,li+1):
        iter_key = "iter_{:08d}".format(iter_idx)
        print("   calculating map for iteration %s"%(iter_idx))
        iter_group = westh5['iterations/iter_{:08d}/'.format(iter_idx)]
        iteration_coords_arr = numpy.array(iter_group['auxdata/coords'])
        # Check the indexing here!
        iter_assignments = numpy.array(assignh5['assignments'][iter_idx-1])
        weights = numpy.array(iter_group['seg_index']['weight'])

        n_particles = iter_assignments.shape[0]
        n_tp = iter_assignments.shape[1]
        # Build a filter array for which walkers fall into one of the selected 
        # bins. A "1" in this array means that the particle of corresponding 
        # index is IN one of the selected bins. 
        state_filter_arr = numpy.zeros((n_particles, n_tp))
        for bin_idx in binlist:
            state_filter_arr = numpy.logical_or(state_filter_arr,
                                                iter_assignments==bin_idx)

        # Also check the indexing in these next few lines.  It's probably
        # not right.                                               V
        desired_structure_idx_list = numpy.where(state_filter_arr)
        if len(desired_structure_idx_list[0]) > 0: 
            struct_idx = (desired_structure_idx_list[0][0],
                          desired_structure_idx_list[1][0])
            weight = weights[struct_idx[0]]
            contact_map = weight*gen_contact_map(iteration_coords_arr[struct_idx].squeeze()) 
            total_weight += weight
            if len(desired_structure_idx_list[0]) > 1: 
                for i in xrange(1, len(desired_structure_idx_list[0])):
                    struct_idx = (desired_structure_idx_list[0][i],
                                  desired_structure_idx_list[1][i])
                    weight = weights[struct_idx[0]]
                    contact_map += weight*\
                            gen_contact_map(iteration_coords_arr[struct_idx].squeeze()) 
                    total_weight += weight
            
    # Normalize by the total weight
    if contact_map is not None:
        contact_map /= total_weight
    else:
        print('No structures found in selected state!')
    return contact_map
    
    

def main():
   
    # Set some initial variables. 
    global distance_reference_array
    global westh5
    global simulation_prefix

    contact_list = load_contact_list('./afCAL.go.parameters')
    Tbins = load_transition_state_bins()
    distance_reference_array = load_distance_reference_array_and_mask(contact_list)
    
    #simulations = [ '69_60_N2NP_1',
    #                '69_60_NP2N_1',
    #                '69_60_N2NP_2',
    #                '69_60_NP2N_2',
    #                '69_60_N2NP_3',
    #                '69_60_NP2N_3',
    #                '69_60_N2NP_4',
    #                '69_60_NP2N_4',
    #                '69_60_N2NP_5',
    #                '69_60_NP2N_5'  ]
    simulations = [ '69_60_N2NP_5',
                    '69_60_NP2N_5'  ]
    
    for simulation_prefix in simulations: 
        westH5 = h5py.File("/mnt/NAS2/SS_SWITCH/RADIAL_SIMS/{:s}_west.h5".format(simulation_prefix),
                           'r')
        assignH5 = h5py.File('../assignments/{:s}_assign.h5'.format(simulation_prefix),
                             'r')


        print("Working on {:s}".format(simulation_prefix))
        T_cm = calculate_average_contact_map(westH5,
                                             assignH5,
                                             binlist=Tbins,
                                             fi=1,
                                             li=2000         )
        numpy.save("contact_maps/{:s}_T_cm".format(simulation_prefix), T_cm) 


        # Check this!
        N_cm = calculate_average_contact_map(westH5,
                                             assignH5,
                                             statelist=[0],
                                             fi=1,
                                             li=1999         )
        numpy.save("contact_maps/{:s}_N_cm".format(simulation_prefix), N_cm) 


        Np_cm = calculate_average_contact_map(westH5,
                                              assignH5,
                                              statelist=[1],
                                              fi=1,
                                              li=1999         )
        numpy.save("contact_maps/{:s}_Np_cm".format(simulation_prefix), Np_cm) 
       
        westH5.close()
        assignH5.close()
    return
    
          
if __name__ == '__main__':
    main()
