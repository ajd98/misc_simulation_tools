#!/usr/bin/env python
import numpy
import shutil
import subprocess

class PFoldSim:
    def __init__(self):
        self.bdpath = "/share/work/ajd98/pfold/tools/uiowa_bd_openmp.32-07-12.BMM.exe"
        self.rmsdfinder_N_path = "/share/work/ajd98/pfold/tools/rmsd_non_perm.py"
        self.rmsdfinder_Np_path = "/share/work/ajd98/pfold/tools/rmsd_perm.py"
        self.reference = '/share/work/ajd98/pfold/reference/folded.pdb'

        self.status_N = 0
        self.status_Np = 1

        self.seedfilepath = './random_seed'
        seedfile = open(self.seedfilepath, 'w+')
        seedfile.close()

    def gen_randomseed(self):
        self.randomseed = numpy.random.randint(2**16)
        self.randomseedfile = open(self.seedfilepath, 'rw+')
        self.randomseedfile.write(str(self.randomseed) + '\n')
        self.randomseedfile.close()
        return

    def _run_bd(self):
        self.gen_randomseed()
        p = subprocess.Popen(['bash','module_load.sh'])
        p.wait()
        siminp = open('sim.inp', 'r')
        p = subprocess.Popen([self.bdpath, 
                              str(self.randomseed),
                              '1', '1'],
                             stdin=siminp
                             )
        p.wait()
        siminp.close()
        return
 
    def _check_status(self, n_rmsd_arr, np_rmsd_arr):
        for i in xrange(n_rmsd_arr.shape[0]):
            if n_rmsd_arr[i] < 2.0 and np_rmsd_arr[i] > 2.0:
                return self.status_N
            if n_rmsd_arr[i] > 2.0 and np_rmsd_arr[i] < 2.0:
                return self.status_Np
        return 'continue'

    def _calc_N_rmsd(self):
         out = subprocess.check_output([self.rmsdfinder_N_path,
                                       self.reference, 
                                       'fort.23']
                                      )
         s = out.split()
         vals = numpy.array([float(val) for val in s])
         return vals

    def _calc_Np_rmsd(self):
         out = subprocess.check_output([self.rmsdfinder_Np_path,
                                       self.reference, 
                                       'fort.23']
                                      )
         s = out.split()
         vals = numpy.array([float(val) for val in s])
         return vals

    def write_rmsds(self):
        n_rmsd_arr = numpy.array(self.n_rmsds)
        np_rmsd_arr = numpy.array(self.np_rmsds)
        numpy.save('n_rmsd', n_rmsd_arr)
        numpy.save('np_rmsd', np_rmsd_arr)
        return

    def run(self):
        not_in_state = True
        iter_idx = 0
        self.n_rmsds = []
        self.np_rmsds = []
        while not_in_state:
            self._run_bd()
            print('finished running BD')
            shutil.copy('fort.23', 'coords/fort.23.{:06d}'.format(iter_idx))
            n_rmsd_arr =  self._calc_N_rmsd()
            np_rmsd_arr =  self._calc_Np_rmsd()
            self.n_rmsds.append(n_rmsd_arr)
            self.np_rmsds.append(np_rmsd_arr)
            status = self._check_status(n_rmsd_arr, np_rmsd_arr)
            self.write_rmsds()
            if status == 0 or status == 1:
                not_in_state = False 
            iter_idx += 1
        endpointfile = open('endpoint_status', 'w+') 
        endpointfile.write(str(status))
        endpointfile.close()
          
def main():
    pfoldsim = PFoldSim()
    pfoldsim.run()

if __name__ == "__main__":
    main()
