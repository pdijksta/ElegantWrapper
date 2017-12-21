import os
import glob
import h5py
import numpy as np

from GenesisWrapper.simulation import InputParser

class ElegantSimulation:
    def __init__(self, input_file):
        if not input_file.endswith('.ele'):
            raise ValueError('File is not an .ele')
        self.input = InputParser(input_file)
        self.directory = os.path.abspath(os.path.dirname(input_file))

        self.rootname = self.input['rootname']
        self.filename = input_file

        self.sig = self._get_FileViewer('%s.sig' % self.rootname)
        self.twiss = self._get_FileViewer('%s.twi' % self.rootname)
        self.out = self._get_FileViewer('%s.out' % self.rootname)
        self.cen = self._get_FileViewer('%s.cen' % self.rootname)

        watch_files = glob.glob(self.directory+'/%s*.w1' % self.rootname)
        self.watch = [self._get_FileViewer(f) for f in watch_files]

    def _get_FileViewer(self, filename):
        raw_file = os.path.join(self.directory, filename)
        processed_file = raw_file + '.h5'

        if not os.path.isfile(processed_file) or \
                os.path.getmtime(processed_file) < os.path.getmtime(raw_file):
            os.system('sdds2hdf %s %s 2>&1 >/dev/null' % (raw_file, processed_file))
            print('Generated %s' % processed_file)

        return FileViewer(processed_file)


class FileViewer:
    def __init__(self, filename):
        if not filename.endswith('.h5'):
            raise ValueError('File is not an .h5')
        self.filename = filename
        self._dict = {}

    def __getitem__(self, key):
        if key not in self._dict:
            with h5py.File(self.filename, 'r') as ff:
                self._dict[key] = np.array(ff['page1/'+key])
        return self._dict[key]

    def print_tree(self):
        def name_and_size(key, ff):
            try:
                print((key, ff[key].shape, ff[key].dtype))
            except:
                print(key)

        with h5py.File(self.filename, 'r') as ff:
            ff.visit(lambda x: name_and_size(x, ff))

