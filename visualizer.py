#!/usr/bin/env python
# Written 7/2015 by Alex DeGrave


################################# FUNCTIONS ###################################
from __future__ import print_function
import argparse
import MDAnalysis
import numpy
import h5py
import re
import time
import os
import subprocess
import sys


class Visualizer:
    def __init__(self):
        '''
        Initialize the Visualizer class.
        '''
        self.argparser()
        self.BBQ_script_path = "./BBQ_script.sh" 
        self.pymol_script_path = "./pymol_script.py" 
        if self.args.westh5_name is not None:
            self.coordinate_file_type = 'westpa'
        if self.args.fort23_name is not None:
            self.coordinate_file_type = 'fort23'
        if self.coordinate_file_type == 'westpa':
            self.westh5_file = h5py.File(self.args.westh5_name,'r')

    def argparser(self):
        '''
        Parses command-line arguments.
        '''    
        parser = argparse.ArgumentParser()
        # Coordinate file types
        parser.add_argument('--westh5', dest='westh5_name', default=None,
                            help='Inpuyt coordinates from the specified '
                                 'west.h5 file. Ex: west.h5')
        parser.add_argument('--fort23', dest='fort23_name', default=None,
                            help='Input coordinates from the specified '
                                 'fort.23 file. Ex: fort.23')
        # WESTPA specific flags 
        parser.add_argument('--iter', dest='root_iteration', default=None, type=int,
                            help='The iteration from which to start tracing. Ex: 500')
        parser.add_argument('--seg-id', dest='root_seg_id', default=None, type=int,
                            help='The segment id of the segment to trace. Ex: 100')
        parser.add_argument('--tp-length', dest='tp_length', default=3, type=int,
                            help='The number of timepoints per iteration. Ex: 3')
        # Topology information
        parser.add_argument('--topo', dest='topology_file', default=None, 
                            help='The topology file. Ex: 3icb.pdb')
        parser.add_argument('--num-atoms', dest='num_atoms', default=113, type=int,
                            help='The number of atoms in the system. Ex: 113')
        # Flag for images
        parser.add_argument('--timepoint-id', dest='tp_id', default=None, type=int,
                            help='The timepoint id of the structure for which'
                                 ' to generate an image. Ex: 2')
        # Output options
        parser.add_argument('--output-mode', dest='output_mode', type=str,
                            help='The output mode.\n    Use ``movie`` for a ' +
                                 'movie of the trajectory up to the iteration' +
                                 ' specified by --iter. \n' +
                                 '    Use ``image`` for a still image of the' +
                                 ' structure specified by ``--iter`` and ' +
                                 '``--seg-id``')
        parser.add_argument('--raw-outdir', dest='raw_output_dir', type=str,
                            help='The output directory for raw pdb files. Ex: ./RAW_PDBS/')
        parser.add_argument('--bbq-outdir', dest='BBQ_output_dir', type=str,
                            help='The output directory for the pdb files with'+
                                 ' added backbones. Ex: ./BBQ_PDBS/')
        self.args = parser.parse_args()
    
    def check_args(self):
        '''
        Check that the supplied arguments are sufficient and that they
        unambiguously inform the Visualizer class of what to do. 
        '''
        if self.args.westh5_name and self.args.fort23_name:
            print('Error.  Only one coordinate file may be specified. \n'
                  'Both ``--westh5 %s`` and ``--fort23 %s were specified.'
                  % (self.args.westh5_name, self.args.fort23_name)) 
            exit(0)

    def centroid(self,coordinate_array):
        '''
        Given a numpy array of x,y,z coordinates, this function returns the 
        centroid of these coordinates (points are weighted equally)
        '''
        sums = numpy.sum(coordinate_array, axis=1)
        return numpy.array((
                            float(sums[0])/float(len(coordinate_array)),
                            float(sums[1])/float(len(coordinate_array)),
                            float(sums[2])/float(len(coordinate_array)) 
                            ))

    def align(self,native_structure,new_structure): 
        '''
        Performs a least-squares superimposition of ``new_structure`` on
        ``native_structure`` and returns an aligned ``new_structure`` array,
        with its centroid at the origin.
        ``native_stucture`` and ``new_structure`` should be numpy arrays with 
        dimension zero as the particle number, and dimension 1 as x,y,z 
        coordinates.
        '''
    
        # Center both native_structure and new_structure on centroid.
        c = self.centroid(native_structure) 
        for i in range(len(native_structure)):
            native_structure[i] = numpy.subtract(native_structure[i], c)
        c = self.centroid(new_structure)
        for i in range(len(new_structure)):
            new_structure[i] = numpy.subtract(new_structure[i], c)
    
        ### Use Kabsch algorithm to calculate optimal rotation matrix. ###
        # Calculate covariance matrix.
        covariance_matrix = numpy.dot(numpy.transpose(native_structure),
                                      new_structure)
    
        # Singular Value Decomposition.
        V, S, Wt = numpy.linalg.svd(covariance_matrix)

        # Correct for the ``handedness`` of the rotation.
        d = numpy.sign(numpy.linalg.det(numpy.dot(numpy.transpose(Wt),
                                                  numpy.transpose(V)))
                       )
        # Calculate the optimal rotation matrix
        U = numpy.dot(numpy.transpose(Wt),
                   numpy.dot(numpy.array(((1,0,0),(0,1,0),(0,0,d))),
                   numpy.transpose(V)))
    
        # Multiplying new_structure (n*3 matrix) by 3*3 optimal rotation matrix
        # ``U`` gives least_squares alignment.
        l_aligned = new_structure.dot(U)
        return l_aligned

    def get_coordinates(self, iteration_id=None, segment_id=None):
        '''
        Return the coordinates from coordinate file.
            For westpa coordinate file: Return 3-dimensional numpy array 
                (dim 0 = timepoint_id, dim 1 = atom, dim 2 = x/y/z) 
                for the set of structures specified by segment_id, iteration_id.
                Data is returned auxdata/coords.
            For fort23 coordinate file: Return 3-dimensional numpy array
                (dim 0 = timepoint_id, dim 1 = atom, dim 2 = x/y/z) 
                for the entire set of structures specified in the fort.23 file.
        '''
        if self.coordinate_file_type == 'westpa':
            coords = self.westh5_file['iterations/iter_%08d/auxdata/coords'\
                                       % iteration_id][segment_id]
            coords = numpy.array(coords)
        if self.coordinate_file_type == 'fort23':
            coords = numpy.loadtxt(self.args.fort23_name)[:,1:] 
            coords = coords.reshape((coords.shape[0],coords.shape[1]/3,3)) 
        return coords  
 
    def trace(self, iteration_id, segment_id):
        '''
        Trace a trajectory up the history tree, starting from segment
        ``segment_id`` of iteration ``iteration_id.``  Return a list of tuples,
        formatted as (iteration_id, segment_id, coordinates). The first element
        corresponds to the highest iteration id, while the last element
        corresponds to iteration id 1.
        '''
        trace = []
        coordinates = self.get_coordinates(iteration_id, segment_id)
        trace.append( (iteration_id, segment_id, coordinates) )
        while iteration_id > 1:
            parent_segment_id = self.westh5_file[
                'iterations/iter_%08d/seg_index' % iteration_id][segment_id][1]
            iteration_id -= 1
            coordinates = self.get_coordinates(iteration_id, parent_segment_id)
            trace.append( (iteration_id, parent_segment_id, coordinates) )
            segment_id = parent_segment_id
            sys.stdout.flush()
            print('\rTracing segment {seg}: working on iteration {i}          '.format(
                  seg=self.args.root_seg_id,i=iteration_id), end='')
        return trace

    def get_traj_coordinates(self):
        '''
        Return the coordinates for the trajectory traced upwards from 
        segment ``segment_id`` of iteration ``iteration_id.``
        '''
        # Trace the trajectory and pull the coordinates.
        if self.coordinate_file_type == 'westpa':
            trajectory = self.trace(self.args.root_iteration,
                                    self.args.root_seg_id)
            # Reverse the trace so that the first element corresponds to iteration 1 
            trajectory.reverse()
            # Make an array of all the coordinates
            coordinates = numpy.vstack([point[2][0:self.args.tp_length-1]\
                                        for point in trajectory]) 

        # Add code for fort.23 files and xtcs here
        elif self.coordinate_file_type == 'fort23':
             coordinates = self.get_coordinates()
        elif self.coordinate_file_type == 'xtc':
            print('Methods for xtc files are not yet implemented.') 
            sys.exit(1) 
         
        return coordinates  

    def write_BBQ_script(self, max_pdb_index):
        '''
        Writes a bash script for running the BBQ program on every structure.
        Input structures are taken from ``self.args.raw_output_dir``, and output
        structures are placed in ``self.args.BBQ_output_dir``
        '''
        BBQ_script = open(self.BBQ_script_path,'w+')
        for i in range(max_pdb_index):
            BBQ_script.write('java -classpath /home/ajd98/apps/BBQ/:'\
                             '/home/ajd98/apps/BBQ/jbcl.jar BBQ '\
                             '-d=/home/ajd98/apps/BBQ/q_50_xyz.dat '\
                             '-r=%s/%06d.pdb 1> %s/%06d.pdb 2>> BBQ.err\n' \
                             % (self.args.raw_output_dir, 
                                i, 
                                self.args.BBQ_output_dir,i))
        BBQ_script.close()
        return

    def run_BBQ_script(self):
        '''
        Runs the BBQ_script written by the function ``write_BBQ_script``
        using subprocessing. 
        Something odd with the method can cause the BBQ output file not to 
        finish.  Running the BBQ script from the command line, however, does
        not produce this error.
        '''
        subprocess.Popen(["bash",self.BBQ_script_path]).wait()
        
        return

    def run_BBQ(self,raw_pdb_name,BBQ_pdb_name=None):
        '''
        Uses subprocess.Popen to run the BBQ utility on the file 
        ``raw_pdb_name`` located in the directory ``self.args.raw_output_dir``.
        A new PDB file is output to 
        ``<self.args.BBQ_output_dir>/<raw_pdb_name>``.
        '''
        if BBQ_pdb_name:
            stdout_file = open("%s/%s"%(self.args.BBQ_output_dir,BBQ_pdb_name),'w')
        else:
            stdout_file = open("%s/%s"%(self.args.BBQ_output_dir,raw_pdb_name),'w')
        stderr_file = open("BBQ.err",'w+') 
        subprocess.Popen(["java","-classpath",
                          "/home/ajd98/apps/BBQ:/home/ajd98/apps/BBQ/jbcl.jar",
                          "BBQ",
                          "-d=/home/ajd98/apps/BBQ/q_50_xyz.dat",
                          "-r=%s/%s"%(self.args.raw_output_dir,raw_pdb_name)],
                          stdout=stdout_file, stderr=stderr_file).wait() 
        return

    class PDBFormatter:
        def __init__(self, input_file_path):
            '''
            Initialize the PDBFormatter class. Read the input file into memory
            and initialize a list ``self.output_lines`` to hold lines for 
            output. 
            '''
            input_file = open(input_file_path,'r')
            self.input_lines = input_file.readlines()
            self.output_lines = self.input_lines 
            return

        def add_header_information(self,coordinate_file_path,current_iteration=None,current_segment=None,current_timepoint=None,root_iteration=None,root_segment=None):
            '''
            Add header lines describing the time of creation and source file for coordinates. 
            '''
            if current_timepoint != None:
                self.output_lines.insert(0, 'HEADER Current timepoint: %d \n' % current_timepoint)
            if current_segment   != None:
                self.output_lines.insert(0, 'HEADER Current segment  : %d \n' % current_segment)
            if current_iteration != None:
                self.output_lines.insert(0, 'HEADER Current iteration: %d \n' % current_iteration)

            self.output_lines.insert(0, 'HEADER **PDB FILE CREATED BY VISUALIZER.PY, WRITTEN BY ALEX DEGRAVE** \n')
            self.output_lines.insert(1, 'HEADER created on %s\n' % time.strftime('%c'))
            self.output_lines.insert(2, 'HEADER from coordinate file %s\n' % os.path.abspath(coordinate_file_path)) 
            if root_iteration and root_segment:
                self.output_lines.insert(3, 'HEADER This structure is part of a trace from iteration %d, segment %d \n' \
                                         % (root_iteration, root_segment)) 
            return

        def add_ss_information(self):
            '''
            Add header lines describing secondary structure for Calbindin-AFF.
            Helical residues are based on PDB 3ICB.
            '''
            self.output_lines.insert(0, 'HELIX    4 III THR      3  ASP     11  1                                  10   \n') 
            self.output_lines.insert(1, 'HELIX    5  IV SER     19  GLN     32  1IRREGULAR, 3/10 FOR 66-70         14   \n') 
            self.output_lines.insert(2, 'HELIX    1   I SER     40  LYS     54  1BECOMES 3/10 RES 13-16            15   \n') 
            self.output_lines.insert(3, 'HELIX    2  II SER     62  PHE     74  1                                  13   \n') 
            self.output_lines.insert(4, 'HELIX    3   L PHE     74  LYS     79  1IRREGULAR                          6   \n') 
            self.output_lines.insert(5, 'HELIX    4 III THR     83  ASP     92  1                                  10   \n') 
            self.output_lines.insert(6, 'HELIX    5  IV SER    100  GLN    113  1IRREGULAR, 3/10 FOR 66-70         14   \n') 
            return

        def make_BBQ_residue_ordering_consistent(self):
            '''
            Rearrange the first four ATOM entries in the PDB file so that 
            ordering is consistent (N, CA, C, O, N, CA, C, O, ...)
            This method assumes that the file contains no header information,
            and should therefore be called before 
            ``self.add_header_information``.
            '''
            self.output_lines[0] = self.input_lines[3]
            self.output_lines[1] = self.input_lines[1]
            self.output_lines[2] = self.input_lines[0]
            self.output_lines[3] = self.input_lines[2]
            return

        def renumber(self):
            '''
            Renumbers the residue and atom indices to be consecutive
            one-indexed integers. Only ``ATOM`` and ``HETATM`` entries are 
            renumbered.
            '''
            atom_index = 1
            residue_index = 1
            for iline, line in enumerate(self.output_lines):
                if (re.match('ATOM',line) or re.match('HETATM',line)):
                    try:
                        if previous_residue_number != int(line[22:26]):
                            residue_index += 1
                    except NameError:
                        pass
                    atom_number_field    = str(atom_index).rjust(7)
                    residue_number_field = str(residue_index).rjust(4)
                    reformatted_line     = line[0:4] + atom_number_field +\
                                           line[11:22] + residue_number_field +\
                                           line[26:]
                    self.output_lines[iline] = reformatted_line
                    atom_index += 1
                    previous_residue_number = int(line[22:26])
            return
        
        def add_TER_record(self):
            '''
            Adds a TER record after the last ``ATOM`` or ``HETATOM`` entry. 
            '''
            TER_index = None
            for iline, line in enumerate(self.output_lines):
                if line.startswith('ATOM') or line.startswith('HETATOM'):
                    TER_index = iline
            if TER_index != None: self.output_lines.insert((TER_index + 1),'TER \n')    
            return

        def add_CONECT_flags(self):
            '''
            Adds CONECT flags to the end of the PDB file.  This method assumes
            that ``ATOM`` entries exist only for the backbone, and are ordered 
            by atom type as:
                     N, CA, C, O, N, CA, C, O, ..., N, CA, C, O
            '''
            num_atoms = len(self.input_lines)
            for i in range(0,(num_atoms)/4):
                 self.output_lines.append('CONECT' + str(4*i+1).rjust(5)+str(4*i+2).rjust(5)+'\n')
                 self.output_lines.append('CONECT' + str(4*i+2).rjust(5)+str(4*i+3).rjust(5)+'\n')
                 self.output_lines.append('CONECT' + str(4*i+3).rjust(5)+str(4*i+4).rjust(5)+'\n')
                 if i != ((num_atoms)/4 -1):
                     self.output_lines.append('CONECT' + str(4*i+3).rjust(5)+str(4*i+5).rjust(5)+'\n')
            return

        def write_to_file(self,output_file_path):
            '''
            Writes ``self.output_lines`` to the file specified by 
            ``output_file_path``.
            '''
            output_file = open(output_file_path,'w')
            for line in self.output_lines:
                output_file.write(line)
            output_file.close()    
            return

        def run(self):
            '''
            Run the formatter; reorder residues, add header information, 
            renumber, and add CONECT flags.  ``self.write_to_file`` must be 
            called separately to write to an output file. 
            '''
            self.make_BBQ_residue_ordering_consistent()
            self.add_ss_information()
            self.renumber()
            self.add_CONECT_flags()
            self.add_TER_record()
            return
            
    def write_and_run_pymol_script(self,max_pdb_index):
        '''
        Writes a pymol script that loads every frame, adds secondary structure
        information for calbindin-AFF (based on PDB 3ICB), highlights
        residues at which fluorophores are attached, and colors the individual
        EF-hands. This function then uses subprocessing to run the script.
        '''
        pymol_script = open(self.pymol_script_path,'w+')
        pymol_script.write("from pymol import cmd\n")
        for i in range(max_pdb_index):
            pymol_script.write("cmd.load('%s/%06d-renumbered.pdb','mov')\n" \
                               % (self.args.BBQ_output_dir,i))
        pymol_script.write("cmd.hide('everything')\n")
        pymol_script.write('''cmd.alter("resi 3-10", "ss='H'")\n''')
        pymol_script.write('''cmd.alter("resi 20-31", "ss='H'")\n''')
        pymol_script.write('''cmd.alter("resi 41-53", "ss='H'")\n''')
        pymol_script.write('''cmd.alter("resi 63-72", "ss='H'")\n''')
        pymol_script.write('''cmd.alter("resi 84-91", "ss='H'")\n''')
        pymol_script.write('''cmd.alter("resi 101-112", "ss='H'")\n''')
        pymol_script.write("cmd.show('spheres', 'id 2')\n")
        pymol_script.write("cmd.show('spheres', 'id 322')\n")
        pymol_script.write("cmd.show('cartoon')\n")
        pymol_script.write("cmd.color('green', 'resi 1-32')\n")
        pymol_script.write("cmd.color('gray', 'resi 33-38')\n")
        pymol_script.write("cmd.color('marine', 'resi 39-81')\n")
        pymol_script.write("cmd.color('firebrick', 'resi 82-113')\n")
        pymol_script.write("cmd.color('yellow', 'id 2')\n")
        pymol_script.write("cmd.color('yellow', 'id 322')\n")
        pymol_script.close()
        subprocess.call(['pymol','-r',self.pymol_script_path])
        return
    


     
    def make_image(self):
        # Get the coordinates of the specified structure. Returns 
        # multiple timesteps, but only the first one is actually used.  This
        # is due to a bug were BBQ does not complete correctly for the last
        # structure provided.
        if self.coordinate_file_type == 'westpa':
            coord_array = self.get_coordinates(self.args.root_iteration, 
                                               self.args.root_seg_id) 
        if self.coordinate_file_type == 'fort23':
            coord_array = self.get_coordinates()
      
        # Prepare the MDAnalysis universe for loading coordinates, by first loading
        # a topology file.
        universe = MDAnalysis.Universe(self.args.topology_file)
        
        # A horribly hacky fix for the weird BBQ errors. 
        print('\nOutputting the selected frame PDB file to ``{0}``'.format(\
              self.args.BBQ_output_dir + '/000000-renumbered.pdb')) 
        temp_coord_array = numpy.empty([2]+[i for i in coord_array.shape[1:]])
        temp_coord_array[0] = coord_array[self.args.tp_id]
        temp_coord_array[1] = coord_array[self.args.tp_id]
        for icoords, coords in enumerate(temp_coord_array): 
	    # Calling ``align`` with the same structure twice simply centers
	    # ``coord_array`` on the origin. Otherwise, structures are actually
            # aligned.
	    aligned_array = self.align(coord_array[0],coords)
	
	    # Write the new pdb
	    writer = universe.trajectory.Writer(
		self.args.raw_output_dir + '/%06d.pdb' % icoords,
		numatoms=self.args.num_atoms)
	    universe.atoms.set_positions(aligned_array)
	    writer.write(universe)
        
        # Write a script for BBQ to run on every structure and fill in the 
        # backbone coordinates, from the C-alpha input coordinates.
        self.write_BBQ_script(coord_array.shape[0])
        
        # Run the BBQ_script
        self.run_BBQ_script()
        #self.run_BBQ('000000.pdb')
    
    
        # Reformat the PDB file
        pdb_formatter = self.PDBFormatter(self.args.BBQ_output_dir+'/%06d.pdb'%0) 
        pdb_formatter.run()
        if self.coordinate_file_type == 'westpa':
            pdb_formatter.add_header_information(self.args.westh5_name,
                              root_iteration=self.args.root_iteration,
                              root_segment=self.args.root_seg_id)
        if self.coordinate_file_type == 'fort23':
            pdb_formatter.add_header_information(self.args.fort23_name)
        pdb_formatter.write_to_file(
            self.args.BBQ_output_dir+'/%06d-renumbered.pdb'%0)
        
        # Run pymol to visualize the trajectory. 
        self.write_and_run_pymol_script(1) 
        return

    def make_movie(self):
        # Trace the trajectory and pull the coordinates.
        coordinates = self.get_traj_coordinates()
      
        # Prepare the MDAnalysis universe for loading coordinates, by first loading
        # a topology file.
        universe = MDAnalysis.Universe(self.args.topology_file)
    
        for i_coord, coord_array in enumerate(coordinates):
    
            # Align the structures
            aligned_array = self.align(coordinates[0],coord_array)
    
            # Write the new pdb
            writer = universe.trajectory.Writer(
                self.args.raw_output_dir + '/%06d.pdb' % i_coord,
                numatoms=self.args.num_atoms)
            universe.atoms.set_positions(aligned_array)
            writer.write(universe)
        
        # Write a script for BBQ to run on every structure and fill in the 
        # backbone coordinates, from the C-alpha input coordinates.
        self.write_BBQ_script(coordinates.shape[0])
        
        # Run the BBQ_script
        #self.run_BBQ_script()
        print('\n')
        for i in range(coordinates.shape[0]):
            self.run_BBQ('%06d.pdb'%i)
            sys.stdout.flush()
            print('\rRunning BBQ: {0}%'.format(int(float(i)/coordinates.shape[0]*100)), end='') 
    
        # Reformat the PDB files
        for i in range(coordinates.shape[0]):
            pdb_formatter = self.PDBFormatter(self.args.BBQ_output_dir+'/%06d.pdb'%i) 
            pdb_formatter.run()
            if self.coordinate_file_type == 'westpa':
                pdb_formatter.add_header_information(self.args.westh5_name,
                                  root_iteration=self.args.root_iteration,
                                  root_segment=self.args.root_seg_id)
            elif self.coordinate_file_type == 'fort23':
                pdb_formatter.add_header_information(self.args.fort23_name)
            pdb_formatter.write_to_file(
                self.args.BBQ_output_dir+'/%06d-renumbered.pdb'%i)
        
        # Run pymol to visualize the trajectory. 
        self.write_and_run_pymol_script(coordinates.shape[0]) 
        return

    def run(self):
        if self.args.output_mode == 'movie':
            self.make_movie()
        elif self.args.output_mode == 'image':
            self.make_image()
        else:
            print("\nUnknown output mode: %s"%self.args.output_mode) 
            print("Please choose either ``movie`` or ``image``")
            exit(1)
        return
            
if __name__ == '__main__':
    visualizer = Visualizer()
    visualizer.run() 
