#!/usr/bin/env python
#
# symmetrycorrect.py
#
# From an input PDB file, create new PDB files with reordered atoms to represent
# every possible labeling of symmetric groups.
#
# Currently, corrects for:
#   PHE symmetry
#   O/OXT symmetry
#
# Files are saved to names such as:
#     0000.pdb, 0001.pdb, 0010.pdb, 0011.pdb, ..., 1111.pdb
# where each 1 represents that a residue's atom labeling was flipped relative to
# the input PDB.  Flips are ordered by (i) residue index, followed by (ii) side
# chain before backbone (to resolve any conflicts).
import argparse
import copy

class Atom(object):
    def __init__(self, pdbstr):
        self.index = int(pdbstr[6:11]) 
        self.atomname = pdbstr[12:16] 
        self.resname = pdbstr[17:20]
        self.chainid = pdbstr[21]  
        self.resid = int(pdbstr[23:26])  
        self.coords = (float(pdbstr[30:38]),
                       float(pdbstr[38:46]),
                       float(pdbstr[46:54]))
        self.other = pdbstr[54:]

    def pdb_entry(self):
        s = "ATOM  "+\
            "{:d}".format(self.index).rjust(5)+" "+\
            "{:s}".format(self.atomname).rjust(4)+" "+\
            "{:s}".format(self.resname).rjust(3)+" "+\
            "{:s}".format(self.chainid)+\
            "{:d}".format(self.resid).rjust(4)+\
            "{:.03f}".format(self.coords[0]).rjust(12)+\
            "{:.03f}".format(self.coords[1]).rjust(8)+\
            "{:.03f}".format(self.coords[2]).rjust(8)+\
            self.other
        return s

class Residue(object):
    def __init__(self):
        self.atoms = [] 
    def add_atom(self, atom):
        if len(self.atoms) == 0:
            self.atoms.append(atom)
        else:
            if atom.resid == self.get_resid() \
            and atom.resname == self.get_resname(): 
                self.atoms.append(atom)

    def get_resname(self):
        return self.atoms[0].resname.strip()

    def get_resid(self):
        return self.atoms[0].resid

    def get_atomnames(self):
        return [atom.atomname.strip() for atom in self.atoms]

    def get_atom_by_name(self, atomname):
        for atom in self.atoms:
            if atom.atomname.strip() == atomname:
                return atom
        raise ValueError("No atom found with name {:s}".format(atomname))
        

#class StructureMeta(type):
#    def __iter__(self):
#        return self.residues.__iter__() 

class Structure(object):
   # __metaclass__ = StructureMeta
    def __init__(self):
        self.residues = []

    def get_resids(self):
        return [res.get_resid() for res in self.residues] 
        
    def add_atom(self, atom):
        # Does a corresponding residue already exist?
        if atom.resid in self.get_resids():
            # That residue already exists, so we need to find it and add the atom
            for ires, res in enumerate(self.residues):
                if res.get_resid() == atom.resid:
                    self.residues[ires].add_atom(atom)
                    return
            raise Exception("Resid in get_resids(), but no such resid found!")
        else:
            self.residues.append(Residue())
            self.residues[-1].add_atom(atom)

    def get_c_terminal_res(self):
        cres = self.residues[0]
        for residue in self.residues:
            if residue.get_resid() > cres.get_resid():
                cres = residue
        return cres

    def _get_internal_resid(self, resid):
        for ires, residue in enumerate(self.residues):
            if residue.get_resid() == resid:
                return ires
        raise ValueError("Residue with index {:s} not found!".format(repr(resid)))

    def flip_labeling(self, resid, fliptype):
        internal_resid = self._get_internal_resid(resid)
        if fliptype == "PHE":
            # Possible issues with mutability
            CD1_coords = self.residues[internal_resid].get_atom_by_name('CD1').coords
            CD2_coords = self.residues[internal_resid].get_atom_by_name('CD2').coords
            CE1_coords = self.residues[internal_resid].get_atom_by_name('CE1').coords
            CE2_coords = self.residues[internal_resid].get_atom_by_name('CE2').coords

            self.residues[internal_resid].get_atom_by_name('CD1').coords = CD2_coords
            self.residues[internal_resid].get_atom_by_name('CD2').coords = CD1_coords
            self.residues[internal_resid].get_atom_by_name('CE1').coords = CE2_coords
            self.residues[internal_resid].get_atom_by_name('CE2').coords = CE1_coords
            return

        elif fliptype == "OXT":
            OXT_coords = self.residues[internal_resid].get_atom_by_name('OXT').coords
            O_coords = self.residues[internal_resid].get_atom_by_name('O').coords
            self.residues[internal_resid].get_atom_by_name('OXT').coords = O_coords
            self.residues[internal_resid].get_atom_by_name('O').coords = OXT_coords
            return

        else:
            raise ValueError("flip type {:s} is not implemented!".format(fliptype))
            
    def write_to_pdb(self, outpath='out.pdb', headerlines=None):
        # Sort by resid
        self.residues.sort(key=lambda res: res.get_resid())

        outfile = open(outpath, 'w+')
        if headerlines is not None:
            for line in headerlines:
                outfile.write(line)
        for residue in self.residues:
            residue.atoms.sort(key=lambda atom: atom.index)
            for atom in residue.atoms:
                outfile.write(atom.pdb_entry())
    
    def __len__(self):
        return len(self.residues)
  

class SymmetryCorrect(object):
    def __init__(self):
        self.output_structures = []
        self._parse_args()
        self._load_pdb()
        self._find_symmetric_groups()
        self._generate_relabelings()

    def _parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('--pdb', dest="input_pdb", required=True)
        self.args = parser.parse_args()

    def _load_pdb(self):
        self.inputstructure = Structure()
        f = open(self.args.input_pdb,'r')
        # Go to first ATOM entry in the pdb
        for line in f:
            if line.startswith('ATOM'):
                atom = Atom(line)
                self.inputstructure.add_atom(atom)

        print('Done parsing input pdb file {:s}. Found {:d} residues.'.format(
              self.args.input_pdb, len(self.inputstructure)))

    def _find_symmetric_groups(self):
        count = 0
        # A list of tuples:  (resid, type) where type is 'PHE' or 'OXT'
        self.sym_groups = []

        # PHE symmetry
        for residue in self.inputstructure.residues:
            if residue.get_resname() == "PHE":
                count += 1
                self.sym_groups.append((residue.get_resid(), 'PHE'))

        # OXT symmetry
        cres = self.inputstructure.get_c_terminal_res()
        if 'OXT' in cres.get_atomnames() and 'O' in cres.get_atomnames():
            count += 1
            self.sym_groups.append((cres.get_resid(), 'OXT'))
        print("Found {:d} symmetric groups for which to correct:".format(len(self.sym_groups)))
        for sym_group in self.sym_groups:
            print("    Resid"+"{:d}".format(sym_group[0]).rjust(4)+\
                  " with symmetry type '"+sym_group[1] +"'.")

    def _generate_relabelings(self):
        n_sym = len(self.sym_groups) 
        fmt_str = "0{:d}b".format(n_sym) # for formatting binary strings like 0001, 0010, 0011, etc.

        # Start doing the flips! We need 2^n_sym in total.
        for i_flip in range(2**n_sym): 
            structure = copy.deepcopy(self.inputstructure)
            flip_str = format(i_flip, fmt_str)
            for i_bit, bit in enumerate(flip_str):
                if bit == '1':
                    structure.flip_labeling(self.sym_groups[i_bit][0], self.sym_groups[i_bit][1]) 

            headerlines = ['HEADER PDB for symmetry corrections\n',
                           'HEADER\n',
                           'HEADER Generated by symmetrycorrection.py.\n'
                           'HEADER Written 2017.07.10 by Alex DeGrave\n'
                           'HEADER\n',
                          ]
            if '1' in flip_str:
                headerlines.append('HEADER The following residues have had atom coordinates switched to account for symmetry:\n')
                for i, sym_group in enumerate(self.sym_groups):
                    if flip_str[i] == '1':
                        headerlines.append('HEADER     resid {:d}; symmetry type {:s}\n'
                                           .format(sym_group[0], sym_group[1]))
            else:
                headerlines.append('HEADER    This file represents the original coordinates.\n')
                headerlines.append('HEADER    no residues have had atomic coordinates altered\n')
                            
            structure.write_to_pdb(flip_str+'.pdb', headerlines=headerlines) 
        
            
             
if __name__ == "__main__":
    SymmetryCorrect()
