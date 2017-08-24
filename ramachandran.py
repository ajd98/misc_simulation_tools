#!/usr/bin/env python
import numpy
from matplotlib import pyplot
from periodickde import PeriodicKDE

class RamachandranPlot(object):
    '''
    Create a Ramachandran plot on a given matplotlib axis object, based on 
    vectors of phi and psi angles. An estimate of the probability density is
    performed using a kernel density that accounts for periodicity of the  
    phi/psi space. The kernel is a truncated cosine function, which is zero
    outside the first wave.

    Note that after instantiating the class, the ContourSet (needed for creation
    of a colorbar) may be accessed via the ``p`` attribute.

    ---------
    Arguments
    ---------

    ax: A matplotlib axis object on which to create the plot.

    phi: A length-n numpy array or list of phi angles. Values should be in
        degrees, and fall between -180 and +180 degrees

    psi: A length-n numpy array or list of psi angles. Values should be in 
        degrees, and fall between -180 and +180 degrees

    cmap: (default='magma_r') A string denoting a matplotlib colormap.

    energyscale: (default=True) If true, plot "- ln P", where P is the estimate
        of the probability density.  Otherwise plot the probability denisty.

    kernelwidth: (default=10) The half-width of the kernel used for the kernel 
        density estimate, in degrees. The kernel density estimate will only be 
        correct if kernelwidth<360.

    gridspacing: (default=5) The spacing in degrees between edges of the grid 
        over which the kernel density estimate is performed. Lower values are
        more computationally intensive (dividing the gridspacing by a factor of
        two increase the computational cost by a factor of four), but provide
        higher-resolution estimates.
    '''
    def __init__(self, ax, phi, psi, cmap='magma_r', energyscale=True, 
                 kernelwidth=10, gridspacing=5):
        kde = PeriodicKDE(phi, psi, wavelength=kernelwidth*2) 
        self.points, self.estimates = kde.evaluate()
        self.ax = ax 
        self.energyscale = energyscale
        self.cmap = cmap
        self.plot_contourf_l()
        return 

    def plot_contourf_l(self):
        '''
        Plot contour levels with black lines between then. 
        '''
        ax = self.ax
        X = self.points[0,:]
        Y = self.points[1,:]
        if self.energyscale:
            Z = -1*numpy.log(self.estimates)
            nanmin = Z[numpy.isfinite(Z)].min()
            nanmax = Z[numpy.isfinite(Z)].max()
            Z -= nanmin 
            Z[numpy.logical_not(numpy.isfinite(Z))] = nanmax*2
            zedges = numpy.arange(0, nanmax+1, 1, dtype=int)
            zedges = numpy.arange(0,11,1)
        else:
            Z = self.estimates
            zedges = numpy.arange(0, Z.max()*1.1, Z.max()/10, dtype=float)
        print(self.estimates.shape)

        # Plot contours.
        p = ax.contourf(X, Y, Z, zedges, 
                        cmap=pyplot.cm.get_cmap(self.cmap, len(zedges)-1)
                        )

        # Plot contour lines.
        ax.contour(X, Y, Z, zedges, 
                   linewidths=0.3,
                   colors='k') 
        self.p = p

