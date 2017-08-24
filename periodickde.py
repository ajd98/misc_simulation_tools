#!/usr/bin/env python
from __future__ import print_function
import numpy

class PeriodicKDE(object):
    '''
    Estimate the probability density over a two-dimensional, periodic space
    (with 360-degree periodicity along both dimensions), using kernel density
    estimation. This class explicitly accounts for the periodicity.

    The kernel is a truncated and shifted cosine function which is zero outside 
    of the central waveform.  

    ---------
    Arguments
    ---------

    phi: A vector (numpy array or list) of values along the first dimension. 
        Values should fall in [-180,180]

    psi: A vector (numpy array or list) of values along the second dimension. 
        Values should fall in [-180,180]

    gridspacing: The spacing between edges of the grid corresponding to points
        at which to evaluate the kernel density estimate.

    wavelength: The wavelength of the cosine function used as the kernel.
    '''
    def __init__(self, phi, psi, gridspacing=2, wavelength=10):
        self.data = numpy.array((phi,psi))
        self.data = numpy.swapaxes(self.data, 0, 1)
        self.gridspacing = gridspacing
        self.wavelength = wavelength

    def _1d_kernel(self, x):
        '''
        Evalute the 1-dimensional kernel at x.
        The kernel is f: [-180,180] -> [0,inf)
                    
            f(x) = cos(x/wavelength*2*2*pi)+1 ; x < wavelength/2
                   0                          ; x >=wavelength/2
                                  .....
                               ...     ...
                             ..           ..
                            .               .
                           .                 .
                         ..                   ..
        ______________...                       ..._______________
                  
                      ^--------wavelength---------^

        Integral = wavelength
        '''
        result = numpy.zeros(x.shape)
        half_w = self.wavelength/2
        nonzero_idxs = numpy.absolute(x)<half_w
        result[nonzero_idxs] = numpy.cos(x[nonzero_idxs]/half_w*2*numpy.pi)+1

        # normalize
        return result/self.wavelength

    def _kernel(self, x):
        '''
        Evaluate the kernel at x; result is product of two one-dimension
        kernels.
        '''
        xs = self._1d_kernel(x[:,0])
        ys = self._1d_kernel(x[:,1])
        return numpy.multiply(xs,ys).sum()

    def evaluate(self):
        '''
        Evaluate the kernel density estimate at each point in a grid.
        '''
        # Image width
        iw = 360/self.gridspacing
        phigrid = numpy.linspace(-540, 540, 3*iw+1) # fencepost condition for endpoint
        psigrid = phigrid
        points = numpy.dstack(numpy.meshgrid(phigrid, psigrid))
        estimate = numpy.empty(points.shape[:2], dtype=float)
        print(points.shape)
        print(self.data.shape)
        total_steps = points.shape[0]*points.shape[1]
        for irow, row in enumerate(points):
            for icol, p in enumerate(row):
                progress = float(irow*len(row)+icol)/total_steps
                print("\r{:.02f}".format(progress), end='')
                difference = self.data - p
                estimate[irow, icol] = self._kernel(difference)
        estimate /=  self.data.shape[0]

        phigrid = numpy.linspace(-180, 180, num=360./self.gridspacing, endpoint=False)
        psigrid = phigrid
        points = numpy.array(numpy.meshgrid(phigrid, psigrid))
        # fold along 0th axis
        print("Shape of estimate: " + repr(estimate.shape))
        print("iw: " + repr(iw))
        estimate[iw:2*iw,:] = estimate[:iw,:] + estimate[iw:2*iw,:] + estimate[2*iw:3*iw,:]

        # fold along 1st axis
        estimate[iw:2*iw,iw:2*iw] = estimate[iw:2*iw,:iw] + estimate[iw:2*iw,iw:2*iw] + estimate[iw:2*iw,2*iw:3*iw]
        return points, estimate[iw:2*iw,iw:2*iw]

