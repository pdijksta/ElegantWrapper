import os
import glob

import numpy as np

from GenesisWrapper.simulation import InputParser
from .watcher import Watcher, FileViewer

class ElegantSimulation:

    comment_chars = ('!',)

    def __init__(self, input_file, _file_=None):
        if not input_file.endswith('.ele'):
            raise ValueError('File is not an .ele')
        if _file_ is not None:
            input_file = os.path.join(os.path.dirname(_file_), input_file)
        self.input = InputParser(input_file, self.comment_chars, {})
        self.directory = os.path.abspath(os.path.dirname(input_file))

        self.rootname = self.input['rootname']
        self.filename = input_file

        self.sig = self._get_FileViewer('%s.sig' % self.rootname)
        self.twiss = self._get_FileViewer('%s.twi' % self.rootname)
        try:
            self.out = self._get_Watcher('%s.out' % self.rootname, s=self.sig['columns/s'][-1])
        except:
            print('No out file!')

        self.cen = self._get_FileViewer('%s.cen' % self.rootname)
        self.mag = self._get_FileViewer('%s.mag' % self.rootname)

        watch_files = glob.glob(self.directory+'/%s*.w1' % self.rootname)
        self.watch = [self._get_Watcher(f) for f in watch_files]

    def _get_FileViewer(self, filename):
        processed_file = self._convert(filename)
        return FileViewer(processed_file, self)

    def _get_Watcher(self, filename, *args, **kwargs):
        processed_file = self._convert(filename)
        return Watcher(processed_file, self, *args, **kwargs)

    def _convert(self, filename):
        raw_file = os.path.join(self.directory, filename)
        processed_file = raw_file + '.h5'

        if not os.path.isfile(processed_file) or \
                os.path.getmtime(processed_file) < os.path.getmtime(raw_file):
            status = os.system('sdds2hdf %s %s 2>&1 >/dev/null' % (raw_file, processed_file))
            if status != 0:
                #raise SystemError(status)
                print(status)
            print('Generated %s' % processed_file)
        return processed_file

    def get_element_position(self, elementName):
        s = self.mag['columns/s']
        indices = np.argwhere(self.mag['columns/ElementName'] == bytes(elementName, 'utf-8')).squeeze()
        return s[indices.min()], s[indices.max()]

