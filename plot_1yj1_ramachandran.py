#!/usr/bin/env python
import matplotlib
import numpy
from matplotlib import pyplot
from ramachandran import RamachandranPlot

class MultiRamachandranPlot(object):
    '''
    Example of how to use the RamachandranPlot object.  Here, the backbone
    dihedral angles for several variants of a protein are analyzed.
    '''
    def __init__(self):
        matplotlib.rcParams['font.size'] = 6
        self.fig, self.axes = pyplot.subplots(3, 3)
        self.fig.subplots_adjust(left=0.2, bottom=0.2, top=0.8, right=0.8, wspace=0.0, hspace=0.0)
        self.fig.set_size_inches(8.8/2.54, 8.8/2.54)
        
        self.datawt = numpy.loadtxt('/mnt/NAS2/UBIQUITIN/analysis/ff15ipq.1ubq/dihedrals.dat',
                                    skiprows=1,
                                    usecols=(1,2,3,4,5,6)
                                    )[::10]
        self.data1yiw = numpy.loadtxt('/mnt/NAS2/UBIQUITIN/analysis/ff15ipq.1yiw/dihedrals.dat',
                                      skiprows=1,
                                      usecols=(1,2,3,4,5,6)
                                      )[::1]
        self.data1yj1 = numpy.loadtxt('/mnt/NAS2/UBIQUITIN/analysis/ff15ipq.1yj1/dihedrals.dat',
                                      skiprows=1,
                                      usecols=(1,2,3,4,5,6)
                                      )[::1]

        RamachandranPlot(self.axes[0,0], self.datawt[:,0], self.datawt[:,1])
        RamachandranPlot(self.axes[0,1], self.datawt[:,2], self.datawt[:,3])
        RamachandranPlot(self.axes[0,2], self.datawt[:,4], self.datawt[:,5])

        RamachandranPlot(self.axes[1,0], self.data1yiw[:,0], self.data1yiw[:,1])
        RamachandranPlot(self.axes[1,1], self.data1yiw[:,2], self.data1yiw[:,3])
        RamachandranPlot(self.axes[1,2], self.data1yiw[:,4], self.data1yiw[:,5])

        RamachandranPlot(self.axes[2,0], self.data1yj1[:,0], self.data1yj1[:,1])
        RamachandranPlot(self.axes[2,1], self.data1yj1[:,2], self.data1yj1[:,3])
        p = RamachandranPlot(self.axes[2,2], self.data1yj1[:,4], self.data1yj1[:,5]).p

        cmapax = self.fig.add_axes([0.82, 0.300, 0.018,0.40])
        self.cbar = self.fig.colorbar(p, cax=cmapax, orientation='vertical')
        self._format()
        pyplot.savefig('ramachandran.pdf')

    def _format(self):
        self.cbar.set_label(u'\u2212RT ln P')
        self.linewidth = 0.8
        for i in range(3):
            for j in range(3):
                ax = self.axes[i,j]
                ax.set_xlim(-180,180)
                ax.set_ylim(-180,180)
                ax.tick_params(direction='out', width=self.linewidth)
                ax.set_xticks((-180,0,180))
                ax.set_yticks((-180,0,180))

                if i != 2:
                    ax.set_xticklabels(['' for k in range(5)]) 
                elif j != 0:
                    ax.set_xticklabels(['', '0', '180'])

                if j != 0:
                    ax.set_yticklabels(['' for k in range(5)]) 
                elif i != 0:
                    ax.set_yticklabels(['-180','0','']) 
        self.axes[2,1].set_xlabel(u'\u03c6')
        self.axes[1,0].set_ylabel(u'\u03c8')

        self.axes[0,0].set_title('Glu34')
        self.axes[0,1].set_title('Gly35/(D)-Gln35')
        self.axes[0,2].set_title('Ile36')

        self.axes[0,0].text(-0.8,0.5,'WT',transform=self.axes[0,0].transAxes, rotation=90, va='center', ha='center')
        self.axes[1,0].text(-0.8,0.5,'Met1Leu',transform=self.axes[1,0].transAxes, rotation=90, va='center', ha='center')
        self.axes[2,0].text(-0.8,0.5,'Met1Leu/Gly35(D)-Gln',transform=self.axes[2,0].transAxes, rotation=90, va='center', ha='center')


if __name__ == "__main__":
    MultiRamachandranPlot()
