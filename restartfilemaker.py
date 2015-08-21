#!/usr/bin/env python
#
# Written Aug 7 2015, by Alex DeGrave

'''
This utility generates a UIOWA-BD restart file (a.k.a., 'restart.file') from a
user-specified line of the fort.23 coordinate file.  See ``--help`` for 
options.
'''

import argparse
import numpy
import sys

class RestartFileMaker:
    def __init__(self):
        '''
        Initialize the RestartFileMaker class.  Parse command line arguments.
        '''
        self._parse_args()
        return

    def _parse_args(self):
        '''
        Parse command line arguments.  See <myself>.py --help for options.
        '''
        parser = argparse.ArgumentParser()
        parser.add_argument('--fort23', 
                            type=str, 
                            default='fort.23', 
                            dest='fort23',
                            help='Use FORT23 as the coordinate input file.'
                            )

        parser.add_argument('--restart-ref',
                            type=str,
                            default='restart.file.initial',
                            dest='restart_ref',
                            help='Use RESTART_REF as a reference restart file. '
                                 'It should be from the same system as the one '
                                 'for which you are creating a new restart '
                                 'file. For example, you may the file '
                                 '``restart.file.initial``.'
                            )

        parser.add_argument('--output',
                            type=str,
                            default='restart.file',
                            dest='outpath',
                            help='Write output file to OUTPATH'
                            )

        parser.add_argument('--line',
                            type=int,
                            default=None,
                            dest='linenum',
                            help='Use coordinates from line number LINENUM '
                                 'of the coordinate input file.  Specify '
                                 'this option or use the ``--time`` option, '
                                 'but not both. By default, this utility will '
                                 'use the final line of the specified '
                                 'coordinate input file.'
                            )

        parser.add_argument('--time',
                            type=float,
                            default=None,
                            dest='time',
                            help='Use coordinates from time ``TIME`` of the '
                                 'coordinate output file.  Specify this option '
                                 'or use the ``--line`` option, but not both.'
                            )
        parser.add_argument('--zero-time',
                            default=False,
                            type=bool,
                            dest='zero_time',
                            help='If True, zero the time line of the restart '
                                 'file.  By default, use the time from the '
                                 'coordinate file that the structure came from.'
                            )
        self.args = parser.parse_args()

        # Check for invalid input combinations.
        if self.args.linenum is not None and self.args.time is not None:
            print('Please specify only one of the options ``--line`` and '
                  '``--time``.')
            exit(1)

    def _open_files(self):
        '''Open input and output files.'''
        self.inputfile = open(self.args.fort23, 'r')
        self.outputfile = open(self.args.outpath, 'w')
        self.reffile = open(self.args.restart_ref, 'r')

    def _close_files(self):
        '''Close input and output files.'''
        self.inputfile.close()
        self.outputfile.close()
        self.reffile.close()

    def _get_coordinates(self):
        '''
        Get the coordinates from the user-specified line or timepoint and
        save them as a numpy array.
        '''
        # By default, go to the last line of the file.
        if self.args.linenum is None and self.args.time is None:
            for i, line in enumerate(self.inputfile):
                coordline = line

        # If the user specified a line number, go to that line.
        elif self.args.linenum is not None:
            existencecheck = False
            for i, line in enumerate(self.inputfile):
                # Use one-indexed lines (for compatibility with Vim, etc)
                if i+1 == self.args.linenum:
                    coordline = line
                    existencecheck = True
            if not existencecheck:
                print('Line number ``%d`` does not exist in the specified '
                      'coordinate input file, ``%s``.'%(self.args.linenum,
                                                        self.args.fort23)
                      )
                exit(1)

        # If the user specified a certain time, go to that time.
        else:
            existencecheck = False
            for i, line in enumerate(self.inputfile):
                time = float(line.split()[0])
                if time == self.args.time:
                    coordline = line
                    existencecheck = True
            if not existencecheck:
                print('Time ``%f`` does not exist in the specified '
                      'coordinate input file, ``%s``.'%(self.args.time,
                                                        self.args.fort23)
                      )
                exit(1)

        # Convert data to a numpy array.
        coordarr = numpy.fromstring(coordline, dtype=float, sep=' ')
        self.time = coordarr[0] # The first value is the time
        coordarr = coordarr[1:] # Remove the time entry
        self.coordarr = coordarr.reshape(coordarr.shape[0]/3, 3)

    def _write_output(self): 
        '''Write coordinates to a new restart file.'''
        self.reffile.readline()
        # The second line of the restart file is somewhat mysterious.  
        # It never seems to change with time, and refers to a nonexistent bead zero.
        # As it is always the same as the cryptic line of the restart.file.initial, 
        # we may simply copy it, without actually having any idea what it means.
        mysterious_line = self.reffile.readline()

        # Write the header line
        if self.args.zero_time:
            self.outputfile.write('time ={:>21.5f} ps\n'.format(0))
        else:
            self.outputfile.write('time ={:>21.5f} ps\n'.format(self.time))

        # Print the cryptic line to the new restart file.
        self.outputfile.write(mysterious_line)
        
        # Now print each bead's position to the file.
        bead_index = 1
        for x, y, z in self.coordarr:
            xf = '{:>20.5f}'.format(x) 
            yf = '{:>20.5f}'.format(y) 
            zf = '{:>20.5f}'.format(z) 
            chain_str = '       1'
            bead_number_str = '{:>8}'.format(bead_index) 
            self.outputfile.write(chain_str + 
                                  bead_number_str + 
                                  xf + yf + zf + '\n')
            bead_index += 1

    def run(self):
        '''
        Parse the input files and write a new restart file.  This is the 
        main public class of the RestartFileMaker class.
        '''
        self._open_files()
        self._get_coordinates()
        self._write_output()
        self._close_files()
        return
    
def main():
    restartfilemaker = RestartFileMaker()
    restartfilemaker.run()

if __name__ == '__main__':
    main()
