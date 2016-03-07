#!/usr/bin/env python
#
# coloraccess.py
#
# Utility for checking the compatibility of color schemes with color vision 
# deficiencies.
# 
# Written 3-6-16 by Alex DeGrave
#
# Color matrices based on http://www.colorjack.com (open-source; link now dead)
#
# 

import argparse
import numpy
import time 
import matplotlib.pyplot as plt
from scipy import misc

class ColorMatrices:
    '''Define matrices for RGB color conversions'''
    normal = numpy.array([[100, 0, 0], 
                          [0, 100, 0],
                          [0,   0, 100]])

    protanopia = numpy.array([[56.667, 43.333, 0], 
                              [55.833, 44.167, 0], 
                              [0, 24.167, 75.833]])

    protanomaly = numpy.array([[81.667, 18.333, 0], 
                               [33.333, 66.667, 0], 
                               [0, 12.5, 87.5]])

    deuteranopia = numpy.array([[62.5, 37.5, 0], 
                                [70, 30, 0], 
                                [0, 30, 70]])

    deuteranomaly = numpy.array([[80, 20, 0], 
                                 [25.833, 74.167, 0], 
                                 [0, 14.167, 85.833]])

    tritanopia = numpy.array([[95, 5, 0], 
                              [0, 43.333, 56.667], 
                              [0, 47.5, 52.5]])

    tritanomaly = numpy.array([[96.667, 3.333, 0], 
                               [0, 73.333, 26.667], 
                               [0, 18.333, 81.667]])

    achromatopsia = numpy.array([[29.9, 58.7, 11.4], 
                                 [29.9, 58.7, 11.4], 
                                 [29.9, 58.7, 11.4]])

    achromatomaly = numpy.array([[61.8, 32, 6.2], 
                                 [16.3, 77.5, 6.2], 
                                 [16.3, 32.0, 51.6]]) 

class ColorAccessTool:
    convmat = {'normal' : numpy.array([[100, 0, 0], 
                                       [0, 100, 0],
                                       [0,   0, 100]]),
               'protanopia' : numpy.array([[56.667, 43.333, 0], 
                                           [55.833, 44.167, 0], 
                                           [0, 24.167, 75.833]]),
               'protanomaly' : numpy.array([[81.667, 18.333, 0], 
                                            [33.333, 66.667, 0], 
                                            [0, 12.5, 87.5]]),
               'deuteranopia' : numpy.array([[62.5, 37.5, 0], 
                                             [70, 30, 0], 
                                             [0, 30, 70]]),
               'deuteranomaly' : numpy.array([[80, 20, 0], 
                                              [25.833, 74.167, 0], 
                                              [0, 14.167, 85.833]]),
               'tritanopia' : numpy.array([[95, 5, 0], 
                                           [0, 43.333, 56.667], 
                                           [0, 47.5, 52.5]]),
               'tritanomaly' : numpy.array([[96.667, 3.333, 0], 
                                            [0, 73.333, 26.667], 
                                            [0, 18.333, 81.667]]),
               'achromatopsia' : numpy.array([[29.9, 58.7, 11.4], 
                                              [29.9, 58.7, 11.4], 
                                              [29.9, 58.7, 11.4]]),
               'achromatomaly' : numpy.array([[61.8, 32, 6.2], 
                                              [16.3, 77.5, 6.2], 
                                              [16.3, 32.0, 51.6]])
               } 
    for label in convmat.keys():
        convmat[label] = convmat[label] * .01
    def __init__(self):
        self.parse_args()
        self.create_image()
    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', type=str, required=True,
                            dest='imagepath',
                            help='Input file.  Read the image located at '
                            'IMAGEPATH.'
                            )
        parser.add_argument('-o', '--output', type=str, required=True,
                            dest='outpath',
                            help='Save tiled image to OUTPATH.'
                            )

        self.args = parser.parse_args()

    def create_image(self):
        orig_arr = misc.imread(self.args.imagepath)
        fig, axes = plt.subplots(3,3)
        labels = ['normal', 'deuteranomaly', 'deuteranopia', 'protanomaly', 
                  'protanopia', 'tritanomaly', 'tritanopia', 
                  'achromatomaly', 'achromatopsia']
        #labels = ['normal', 'protanopia', 'protanomaly', 'deuteranopia',
        #          'deuteranomaly', 'tritanopia', 'tritanomaly', 
        #          'achromatopsia', 'achromatomaly']
        for idx, label in enumerate(labels):
            iplot = idx//3
            jplot = idx%3
            cm = self.convmat[label]
            new_arr = self.process(orig_arr, cm)

            ax = axes[iplot, jplot]
            ax.imshow(new_arr)
            ax.tick_params(
                axis='both',
                which='both',
                bottom='off',
                top='off',
                left='off',
                right='off',
                labelbottom='off',
                labelleft='off'
                           )
            ax.set_title(label)
            
        fig = plt.gcf()
        fig.set_size_inches(10,8)
        plt.savefig(self.args.outpath, dpi=300)

    def process(self, imarr, rgb_T):
        '''Apply the RGB conversion matrix to imarr and return a new matrix'''
        newarr = numpy.zeros(imarr.shape, dtype=numpy.uint8)
        if newarr.shape[2] == 4:
            rgb_T = numpy.array([[rgb_T[0,0], rgb_T[0,1], rgb_T[0,2], 0],
                                 [rgb_T[1,0], rgb_T[1,1], rgb_T[1,2], 0],
                                 [rgb_T[2,0], rgb_T[2,1], rgb_T[2,2], 0],
                                 [0         , 0         , 0         , 1]],
                                dtype=numpy.float64) 
        st = time.time()
        newarr = numpy.einsum('ijk,lk', imarr, rgb_T)
        newarr = numpy.array(newarr, dtype=numpy.uint8)
        return newarr
        
if __name__ == '__main__':
    coloraccesstool = ColorAccessTool() 
