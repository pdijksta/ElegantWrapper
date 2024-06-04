import shutil
import os
import glob

import numpy as np

from GenesisWrapper.simulation import InputParser
from .watcher import Watcher, FileViewer, sdds2hdf

class ElegantWrapperError(Exception):
    pass

class ElegantSimulation:

    comment_chars = ('!',)

    def __init__(self, input_file, _file_=None, rootname=None, add_watch=None, del_sim=False):
        self.del_sim = del_sim
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
        try:
            self.bun = self._get_Watcher('%s.bun' % self.rootname, s=0.)
        except:
            self.bun = None

        try:
            self.par = self._get_FileViewer('%s.par' % self.rootname)
        except:
            self.par = None

        try:
            self.opt = self._get_FileViewer('%s.opt' % self.rootname)
        except:
            self.opt = None

        matfile = os.path.join(self.directory, '%s.mat' % self.rootname)
        if os.path.isfile(matfile):
            self.mat = self._get_FileViewer('%s.mat' % self.rootname)
        else:
            self.mat = None

        self.cen = self._get_FileViewer('%s.cen' % self.rootname)
        self.mag = self._get_FileViewer('%s.mag' % self.rootname)

        watch_files = list(glob.glob(self.directory+'/%s.*.w1' % self.rootname))
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
        self.watch = watch

    def __del__(self):
        if self.del_sim:
            directory = os.path.dirname(self.filename)
            if os.path.isdir(directory):
                shutil.rmtree(directory)
                print('Deleted %s' % directory)

    def __repr__(self):
        return os.path.basename(self.filename)
    __str__ = __repr__

    def _get_FileViewer(self, filename):
        try:
            processed_file = self._convert(filename)
            return FileViewer(processed_file, self)
        except ElegantWrapperError:
            #import pdb; pdb.set_trace()
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
        if indices.shape == ():
            indices = np.array([indices])
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

    def get_k1ls(self, combine_q1q2=True, remove_mqup=True):
        quad_indices = self.par['ElementType'] == 'QUAD'
        k1_indices = self.par['ElementParameter'] == 'K1'
        length_indices = self.par['ElementParameter'] == 'L'
        quad_names = np.unique(self.par['ElementName'][quad_indices])
        k1 = np.zeros_like(quad_names, dtype=float)
        length = k1.copy()
        for ctr, quad_name in enumerate(quad_names):
            k1_index = np.logical_and(self.par['ElementName'] == quad_name, k1_indices)
            assert k1_index.sum() == 1
            k1[ctr] = self.par['ParameterValue'][k1_index]

            length_index = np.logical_and(self.par['ElementName'] == quad_name, length_indices)
            assert length_index.sum() == 1
            length[ctr] = self.par['ParameterValue'][length_index]

        if remove_mqup:
            indices = np.array(['MQUP' in x for x in quad_names], bool)
            quad_names = quad_names[~indices]
            k1 = k1[~indices]
            length = length[~indices]

        if combine_q1q2:
            quad_names0 = quad_names
            k10 = k1
            length0 = length
            quad_names, k1, length = [], [], []
            for quad_name, this_k1, this_length in zip(quad_names0, k10, length0):
                if quad_name.endswith('.Q1'):
                    quad_names.append(quad_name[:-3])
                    k1.append(this_k1)
                    length.append(this_length*2)
                elif quad_name.endswith('.Q2'):
                    pass
                else:
                    quad_names.append(quad_name)
                    k1.append(this_k1)
                    length.append(this_length)
            quad_names = np.array(quad_names)
            k1 = np.array(k1)
            length = np.array(length)

        return {
                'quad_names': quad_names,
                'k1': k1,
                'length': length,
                'k1l': k1*length,
                }

