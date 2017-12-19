import os
import h5py
import numpy as np

class ElegantSimulation:
    def __init__(self, input_file):
        if not input_file.endswith('.ele'):
            raise ValueError('File is not an .ele')
        self.basename = input_file[:-4]
        self.filename = input_file

        self.sig = self._get_FileViewer('.sig')
        self.twiss = self._get_FileViewer('.twi')
        self.out = self._get_FileViewer('.out')
        self.cen = self._get_FileViewer('.cen')

    def _get_FileViewer(self, ending):
        raw_file = self.basename + ending
        processed_file = raw_file + '.h5'

        if not os.path.isfile(processed_file):
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

