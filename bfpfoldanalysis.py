#!/usr/bin/env python
import math
import matplotlib.pyplot
import numpy
import os
import scipy.optimize

class PFoldAnalysis:
    '''
    Evaluate similarity of brute-force and WE committor probabilites using 
    Clopper-Pearson CI.
    '''
    def __init__(self, dirlist, alpha=0.05):
        self.dirlist = dirlist
        self.success_status = 0
        self.alpha = alpha
       
    def nCr(self, n, r):
        '''n choose r'''       
        if r > n:
            raise ValueError("nCr cannot be evaluated for r>n!")
        return math.factorial(n)/(math.factorial(r)*math.factorial(n-r))

    def binomial_pm(self, p, k, n):
        '''Return probability mass at k for binomial distribution with n
        samples and success probability p'''
        return nCr(n,k)*(p**k)*(1-p)**(n-k) 

    def error_in_p(self, p, k, n, alpha):  
        '''Returns a function representing the error in the probability 
        estimate'''
        return abs(self.binomial_pm(p, k, n) - alpha/2.0)

    def clopper_pearson(self, k, n, alpha):
        '''Find the Clopper-Pearson CI of width (1-alpha) based on an
        experiment with k successes out of n trials.'''
        x_bar = float(k)/n
        lb = scipy.optimize.fminbound(self.error_in_p, 0, x_bar, 
                                      args=(k, n, self.alpha))
        ub = scipy.optimize.fminbound(self.error_in_p, x_bar, 1, 
                                      args=(k, n, self.alpha))
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
            if status[0] == self.success_status:
                successes += 1
            trials += 1
        return successes, trials

    def load_WE_pfold(self, directory):
        pfold = numpy.loadtxt(os.path.join(directory, 'pfoldWE.txt'))
        return pfold 

    def plot(self, xs, ys, errs):
        ax, fig = pyplot.subplots(figsize=(3.25, 3.25))
        ax.errorbar(xs, ys, yerr=errs, fmt='0', color='black') 
        ax.plot([0,1], [0,1], color='black', ls='--')
        ax.set_xlabel(r'$p_{fold}$ estimate from WE')
        ax.set_ylabel(r'$p_{fold}$ estimate from brute force')
        ax.set_xticks(numpy.arange(0,1.2, 0.2))
        ax.set_yticks(numpy.arange(0,1.2, 0.2))
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
            ys[i] = float(k)/n
            errs[0,i] = lb
            errs[1,i] = ub
        plot(xs, ys, errs)
        return
            

def main():
    dirs = ["00_001472_325_0",
            "02_000970_449_2",
            "03_001426_386_1",
            "10_000723_346_2",
            "13_001215_315_1",
            "19_000585_402_2",
            "02_000344_344_0",
            "02_001442_300_1",
            "09_000980_263_1",
            "10_000727_356_1",
            "18_000625_317_1",
            "19_000633_470_2"]

    pfoldanalysis = PFoldAnalysis(dirs)
    pfoldanalysis.run()

if __name__ == '__main__':
    main()
