import os
import glob

import numpy as np

from GenesisWrapper.simulation import InputParser
from .watcher import Watcher, FileViewer

class ElegantWrapperError(Exception):
    pass

class ElegantSimulation:

    comment_chars = ('!',)

    def __init__(self, input_file, _file_=None, rootname=None):
        if not input_file.endswith('.ele'):
            raise ElegantWrapperError('File is not an .ele')
        if _file_ is not None:
            input_file = os.path.join(os.path.dirname(_file_), input_file)
        self.input = InputParser(input_file, self.comment_chars, {})
        self.directory = os.path.abspath(os.path.dirname(input_file))

        if rootname is None:
            self.rootname = self.input['rootname']
        else:
            self.rootname = rootname
        self.filename = input_file

        self.sig = self._get_FileViewer('%s.sig' % self.rootname)
        if self.sig == 'no_file':
            import pdb; pdb.set_trace()
        self.twiss = self.twi = self._get_FileViewer('%s.twi' % self.rootname)
        try:
            self.out = self._get_Watcher('%s.out' % self.rootname, s=self.sig['columns/s'][-1])
        except:
            self.out = None
            print('No out file!')
        try:
            self.bun = self._get_Watcher('%s.bun' % self.rootname, s=0.)
        except:
            self.bun = None
            print('No bun file!')

        matfile = os.path.join(self.directory, '%s.mat' % self.rootname)
        if os.path.isfile(matfile):
            self.mat = self._get_FileViewer('%s.mat' % self.rootname)
        else:
            self.mat = None

        self.cen = self._get_FileViewer('%s.cen' % self.rootname)
        self.mag = self._get_FileViewer('%s.mag' % self.rootname)

        watch_files = glob.glob(self.directory+'/%s-*.w1' % self.rootname)
        self.watch = [self._get_Watcher(f) for f in watch_files]

    def _get_FileViewer(self, filename):

        try:
            processed_file = self._convert(filename)
            return FileViewer(processed_file, self)
        except ElegantWrapperError:
            return 'no_file'

    def _get_Watcher(self, filename, *args, **kwargs):
        processed_file = self._convert(filename)
        return Watcher(processed_file, self, *args, **kwargs)

    def _convert(self, filename):
        raw_file = os.path.join(self.directory, filename)
        processed_file = raw_file + '.h5'

        if not os.path.isfile(raw_file):
            raise ElegantWrapperError('File %s does not exist!' % raw_file)

        if not os.path.isfile(processed_file) or \
                os.path.getmtime(processed_file) < os.path.getmtime(raw_file):
            cmd = 'sdds2hdf %s %s 2>&1 >/dev/null' % (raw_file, processed_file)
            status = os.system(cmd)
            if status != 0:
                raise SystemError('Status %i for file %s\nCommand: %s' % (status, filename, cmd))
            print('Generated %s' % processed_file)
        return processed_file

    def get_element_position(self, elementName, mean=False):
        s = self.mag['columns/s']
        indices = np.argwhere(self.mag['columns/ElementName'] == elementName).squeeze()
        if mean:
            return (s[indices.max()]+s[indices.min()])/2.
        else:
            return s[indices.min()], s[indices.max()]

    def get_watcher_list(self, key):
        output = []
        for w in self.watch:
            output.append(w[key])
        return np.array(output)

    def get_watcher_func(self, funcname, *args, **kwargs):
        output = []
        for w in self.watch:
            func = getattr(w, funcname)
            output.append(func(*args, **kwargs))
        return np.array(output)

    def get_geometric_emittance(self, dimension):
        assert dimension in ('x', 'y')
        get_entry = lambda x: self.sig['columns/'+x]
        if dimension == 'x':
            x, xp, xxp = get_entry('s1')**2, get_entry('s2')**2, get_entry('s12')
        elif dimension == 'y':
            x, xp, xxp = get_entry('s3')**2, get_entry('s4')**2, get_entry('s34')

        return np.sqrt(x*xp - xxp**2)

