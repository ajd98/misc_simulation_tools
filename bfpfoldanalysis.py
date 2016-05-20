#!/usr/bin/env python
import math
import matplotlib.pyplot as pyplot
import numpy
import os
import scipy.optimize
import scipy.stats

class PFoldAnalysis:
    '''
    Evaluate similarity of brute-force and WE committor probabilites using 
    Clopper-Pearson CI.
    '''
    def __init__(self, dirlist, alpha=0.05):
        self.dirlist = dirlist
        self.success_status = 0
        self.alpha = alpha

    def clopper_pearson(self, k, n, alpha):
        lb = scipy.stats.beta.ppf(alpha/2, k, n-k+1)
        ub = scipy.stats.beta.ppf(1-alpha/2, k+1, n-k)
        if k==0:
            lb = 0
        if k==n:
            ub = 1
        return lb, ub

    def find_successes(self, directory):
        '''
        Return the number of successes (endpoint status = success_status) and
        number of trials given a path (d) to a directory of pfold simulations.
        '''
        dirs = [os.path.join(directory, d) for d in os.listdir(directory)\
                if os.path.isdir(os.path.join(directory, d))]
        successes = 0
        trials = 0
        for d in dirs:
            status = numpy.loadtxt(os.path.join(d, 'endpoint_status')) 
            if status == self.success_status:
                successes += 1
            trials += 1
        return successes, trials

    def load_WE_pfold(self, directory):
        pfold = numpy.loadtxt(os.path.join(directory, 'pfoldWE.txt'))
        return pfold 

    def plot(self, xs, ys, errs):
        fig, ax = pyplot.subplots(figsize=(3.25, 3.25))
        #ax.errorbar(xs, ys, yerr=errs, fmt='o', color='black') 
        ax.errorbar(xs[:6], ys[:6], yerr=errs[:,:6], fmt='o', color='black') 
        ax.errorbar(xs[6:], ys[6:], yerr=errs[:,6:], fmt='o', color='#aaaaaa') 
        ax.plot([0,1], [0,1], color='black', ls='--')
        ax.set_xlabel(r'$p_{N}$ estimate from WE')
        ax.set_ylabel(r'$p_{N}$ estimate from brute force')
        ax.set_xticks(numpy.arange(0,1.2, 0.2))
        ax.set_yticks(numpy.arange(0,1.2, 0.2))
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.tick_params(direction='out', top='off', right='off', width=1.5)
        for kw in ['left', 'right', 'top', 'bottom']:
            ax.spines[kw].set_linewidth(1.5)
        fig.subplots_adjust(left=0.2, bottom=0.2)
        pyplot.savefig('pfold_comparison.pdf')
        return
 
    def run(self):
        xs = numpy.zeros(len(self.dirlist)) 
        ys = numpy.zeros(len(self.dirlist)) 
        errs = numpy.zeros((2, len(self.dirlist))) 
        for i, d in enumerate(self.dirlist): 
            k, n = self.find_successes(d)
            lb, ub = self.clopper_pearson(k, n, self.alpha)
            we_estimate = self.load_WE_pfold(d) 
 
            xs[i] = we_estimate
            mean = float(k)/n
            ys[i] = mean 
            errs[0,i] = mean-lb
            errs[1,i] = ub-mean

            print("mean: {:.05f};    le: {:.05f};     ue: {:.05f}".format(mean, mean-lb, ub-mean))
        self.plot(xs, ys, errs)
        return
            

def main():
    dirs = ["00_001472_325_0",
            "02_000970_449_2",
            "03_001426_386_1",
            "02_000344_344_0",
            "02_001442_300_1",
            "09_000980_263_1",
            "10_000723_346_2",
            "13_001215_315_1",
            "19_000585_402_2",
            "10_000727_356_1",
            "18_000625_317_1",
            "19_000633_470_2"]

    pfoldanalysis = PFoldAnalysis(dirs)
    pfoldanalysis.run()

if __name__ == '__main__':
    main()
