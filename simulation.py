import os
import glob

import numpy as np

from GenesisWrapper.simulation import InputParser
from .watcher import Watcher, FileViewer, sdds2hdf

class ElegantWrapperError(Exception):
    pass

class ElegantSimulation:

    comment_chars = ('!',)

    def __init__(self, input_file, _file_=None, rootname=None, add_watch=None):
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

        try:
            self.sig = self._get_FileViewer('%s.sig' % self.rootname)
        except:
            #print('Weird sig error, continue!')
            self.sig = None
        self.twiss = self.twi = self._get_FileViewer('%s.twi' % self.rootname)
        try:
            self.out = self._get_Watcher('%s.out' % self.rootname, s=self.sig['columns/s'][-1])
        except:
            self.out = None
            #print('No out file!')
        try:
            self.bun = self._get_Watcher('%s.bun' % self.rootname, s=0.)
        except:
            self.bun = None
            #print('No bun file!')

        try:
            self.par = self._get_FileViewer('%s.par' % self.rootname)
        except:
            self.par = None
            #print('No par file!')

        matfile = os.path.join(self.directory, '%s.mat' % self.rootname)
        if os.path.isfile(matfile):
            self.mat = self._get_FileViewer('%s.mat' % self.rootname)
        else:
            self.mat = None

        try:
            self.cen = self._get_FileViewer('%s.cen' % self.rootname)
        except:
            self.cen = None
        self.mag = self._get_FileViewer('%s.mag' % self.rootname)

        watch_files = list(glob.glob(self.directory+'/%s[-_]*.w1' % self.rootname))
        if add_watch:
            watch_files.extend(add_watch)

        watch = [self._get_Watcher(f) for f in watch_files]
        full_path_list = [w.filename for w in watch]
        self.watch = watch = [watch[nn] for nn, path in enumerate(full_path_list) if path not in full_path_list[:nn]]
        if watch:
            # Only unique entries in watch list
            s_list = np.array(self.get_watcher_list('s'))
            try:
                watch = list(list(zip(*sorted(zip(s_list, watch))))[1])
            except Exception as e:
                print(e)
                import pdb; pdb.set_trace()
        self.watch = watch

    def __repr__(self):
        return os.path.basename(self.filename)
    __str__ = __repr__

    def _get_FileViewer(self, filename):
        try:
            processed_file = self._convert(filename)
            return FileViewer(processed_file, self)
        except ElegantWrapperError:
            return 'no_file'

    def _get_Watcher(self, filename, *args, **kwargs):
        processed_file = self._convert(filename)
        w = Watcher(processed_file, self, *args, **kwargs)
        return w

    def _convert(self, filename):
        raw_file = os.path.join(self.directory, filename)
        processed_file = raw_file + '.h5'

        sdds2hdf(raw_file, processed_file)
        return processed_file

    def get_element_position(self, elementName, mean=False):
        s = self.mag['columns/s']
        indices = np.argwhere(self.mag['columns/ElementName'] == elementName).squeeze()
        if len(indices) == 0:
            raise KeyError('%s not found' % elementName)
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

