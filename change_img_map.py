#!/usr/bin/env python

import pylab
import numpy

class ColorMapChanger:
    def __init__(self,input_file_name, input_cmap, output_cmap):
        self.input_file_name = input_file_name 
        self.input_cmap = input_cmap
        self.output_cmap = output_cmap
        self.imdata = pylab.imread(self.input_file_name)
        self.outdata = numpy.zeros(self.imdata.shape, dtype=float)
        self.N = 256
        self.threshold = .01

    def _nearest(self, _lut, cval):
        diff = abs(_lut-cval)
        sqrdiff = numpy.multiply(diff,diff)
        sums = sqrdiff.sum(axis=1)
        idx = sums.argmin()
        return idx, sums

    def _change_cmap(self):
        # Get colormaps with self.N discrete values.
        icmap = pylab.cm.get_cmap(self.input_cmap, self.N)
        # make self._lut available.
        icmap._init()
        ocmap = pylab.cm.get_cmap(self.output_cmap, self.N)
        ocmap._init()

        for i in range(self.imdata.shape[0]):
            for j in range(self.imdata.shape[1]):
                # Check if pixel is grayscale
                if numpy.all(self.imdata[i,j,:] == self.imdata[i,j,0]):
                    new_cvals = self.imdata[i,j]

                # If not grayscale
                else:
                    idx, sums = self._nearest(icmap._lut, self.imdata[i,j])
                    if sums[idx] < self.threshold:
                        new_cvals = ocmap._lut[idx]
                    else:
                        new_cvals = self.imdata[i,j]
                self.outdata[i,j] = new_cvals 

    def _save(self):
        pylab.imsave('test_out.png',self.outdata)

    def run(self):
        self._change_cmap()
        self._save()

def main():
    cmapchanger = ColorMapChanger('/home/ajd98/brandon.png', 'jet', 'hot')
    cmapchanger.run()

if __name__ == '__main__':
    main()

        
