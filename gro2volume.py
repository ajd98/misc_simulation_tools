#!/usr/bin/env python
#
# gro2volume.py
#
# Calculate volume of triclinic box, based on vectors specified in a .gro file
# Syntax: python gro2volume.py -i <myfile>.gro
#
# Written 17.09.07 by Alex DeGrave
import numpy
import argparse

class VolumeCalc(object):
    def __init__(self):
        self._parse_args()
        self.read_box_vectors()
        self.calculate_volume()

    def _parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', dest='inputpath',
                            help='The .gro file containing the box vectors on '
                                 'which the volume calculation is based.'
                            )
        self.args = parser.parse_args()

    def read_box_vectors(self):
        f = open(self.args.inputpath, 'r')
        for line in f:
            pass

        # Box vectors in nm
        # xx, yy, zz, xy, xz, yx, yz, zx, zy
        self.box_vectors = numpy.array([float(v) for v in line.split()])
        f.close()

    def calculate_volume(self):
        M = self.box_vectors

        a = M[0]
        b = M[1]
        c = M[2]
        xy = M[3]
        xz = M[4]
        yx = M[5]
        yz = M[6]
        zx = M[7]
        zy = M[8]

        # See hoomd-blue.readthedocs.io/en/stable/box.html
        # and 
        # https://www.fxsolver.com/browse/formulas/Triclinic+crystal+system+(Unit+cell's+volume)
        # for references on the calculation

        # Cosine of alpha
        cos_alpha = (xy*xz+yz)/numpy.sqrt((1+xy**2)*(1+xz**2+yz**2))
        # Cosine of beta
        cos_beta = xz/numpy.sqrt(1+xz**2+yz**2)
        # Cosine of theta
        cos_gamma = xy/numpy.sqrt(1+xy**2)

        V = a*b*c*numpy.sqrt(1-cos_alpha**2-cos_beta**2-cos_gamma**2+2*cos_alpha*cos_beta*cos_gamma)
        print("Volume of box is {:.05f} nm^3, corresponding to {:.03f} Angstroms^3".format(V,V*1000))


if __name__ == "__main__":
    VolumeCalc()
