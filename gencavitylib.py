#!/usr/bin/env python
'''
Generate a cavity.lib file for use with Voidoo*, based on an Amber parm7 
topology file.

*http://www.msg.ucsf.edu/local/programs/ono/manuals/voidoo_man.html

-----------
How to use:
  Import the CavityLibGen object into a Python script, then initialize the 
  object as CavityLibGen(parmpath), where parmpath is a (str) file path to
  the parm7 file.

  The cavity.lib file will be written to "cavity.lib.autogen".
'''
import numpy

class CavityLibGen(object):
    def __init__(self, parmpath):
        '''
        parmpath: path to Amber parm file
        '''
        self.parmpath = parmpath
        self._load()
        self.build_map()
        self.print_cavity_lib()

    def _load_atom_name(self):
        arr = []
        found_next_section = False

        # skip over the format line
        self.lineidx += 2
        line = self.lines[self.lineidx]

        linecount = 0
        while not found_next_section:
            linecount += 1
            split = [line[i:i+4] for i in range(0, len(line)-1, 4)]
            arr += split
            self.lineidx += 1
            line = self.lines[self.lineidx]
            if line.startswith('%'):
                found_next_section = True 

        return numpy.array(arr)


    def _load_until_next_section(self, dtype=None):
        arr = []
        found_next_section = False

        # skip over the format line
        self.lineidx += 2
        line = self.lines[self.lineidx]
        while not found_next_section:
            split = line.split()
            arr += split
            self.lineidx += 1
            line = self.lines[self.lineidx]
            if line.startswith('%'):
                found_next_section = True 

        if dtype is not None:
            return numpy.array(arr, dtype=dtype)
        else:
            return numpy.array(arr)
            
    def _load(self):
        parmfile = open(self.parmpath, 'r')
        self.lines = parmfile.readlines()
        parmfile.close()
        self.lineidx = 0
        while self.lineidx < len(self.lines):
            line = self.lines[self.lineidx]
            if line.startswith('%FLAG ATOM_TYPE_INDEX'):
                self.atom_type_index = self._load_until_next_section(dtype=int)
            elif line.startswith('%FLAG NONBONDED_PARM_INDEX'):
                self.nonbonded_parm_index = self._load_until_next_section(dtype=int)
            elif line.startswith('%FLAG LENNARD_JONES_ACOEF'):
                self.acoef = self._load_until_next_section(dtype=float)
            elif line.startswith('%FLAG LENNARD_JONES_BCOEF'):
                self.bcoef = self._load_until_next_section(dtype=float)
            elif line.startswith('%FLAG RESIDUE_LABEL'):
                self.residue_label = self._load_until_next_section()
            elif line.startswith('%FLAG ATOM_NAME'):
                self.atom_name = self._load_atom_name()
            elif line.startswith('%FLAG RESIDUE_POINTER'):
                self.residue_pointer = self._load_until_next_section(dtype=int)
            else: self.lineidx += 1
        self.ntypes = int(numpy.sqrt(self.nonbonded_parm_index.shape[0]))



    def get_residue_label(self,iatom):
        '''
        iatom: (int) the index of atom
        '''
        residx = 0

        if iatom < self.residue_pointer[-1]:
            for residx, respointer in enumerate(self.residue_pointer):
                if iatom < self.residue_pointer[residx+1]:
                    break
        else:
            residx = self.residue_pointer.shape[0]-1

        return self.residue_label[residx]

    def get_atom_name(self, iatom):
        '''
        iatom: (int) the index of the atom
        '''
        return self.atom_name[iatom]

    def get_nonbonded_parm_index(self, atomtypeindex1, atomtypeindex2):
        '''
        '''
        return self.nonbonded_parm_index[self.ntypes*(atomtypeindex1-1)+atomtypeindex2-1]


    def get_van_der_waals_radius(self, atomtypeindex):
        '''
         
        '''
        nonbonded_parm_index = self.get_nonbonded_parm_index(atomtypeindex,
                                                             atomtypeindex)
        if nonbonded_parm_index > 0:
            acoef = self.acoef[nonbonded_parm_index-1]
            bcoef = self.bcoef[nonbonded_parm_index-1]
            r = (2*acoef/bcoef)**(1./6)/2
        else:
            r = 0
        if numpy.isnan(r):
            print("Radius is zero for atom type index: {:d}".format(atomtypeindex) )
            r = 0
        return r

    def build_map(self):
        '''
        Build a map between a the pair (residue_label, atom_name) and the van 
        der Waals radius.
        '''
        print("length of ATOM_NAME: {:d}".format(self.atom_name.shape[0]))
        print("length of ATOM_TYPE_INDEX: {:d}".format(self.atom_type_index.shape[0]))
        self.map = dict()
        # iterate over every atom
        for iatom, atomtypeindex in enumerate(self.atom_type_index):
            reslabel = self.get_residue_label(iatom)
            atomname = self.get_atom_name(iatom)
            vdw_r = self.get_van_der_waals_radius(atomtypeindex)
            self.map[(reslabel, atomname)] = vdw_r

    def print_cavity_lib(self):
        outfile = open('cavity.lib.autogen','w+')
        for (reslabel, atomname), vdw_r in self.map.iteritems():
            atomname = atomname.strip()
            if len(atomname) == 1:
                atomname = " {:s}  ".format(atomname)
            elif len(atomname) == 2:
                atomname = " {:s} ".format(atomname)
            elif len(atomname) == 3:
                atomname = " {:s}".format(atomname)
            outfile.write("SPAT '{:s}*{:s}' {:.03f}\n".format(reslabel, atomname, vdw_r)) 
        
        outfile.close()

if __name__ == "__main__":
   CavityLibGen('VILLIN.parm7')
