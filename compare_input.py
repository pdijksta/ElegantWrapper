import random
import os
import shutil
from . import watcher

def compare(files):
    converted_files = []
    base_dir = '/tmp/compare_input_%i' % random.randint(1000,10000)
    os.makedirs(base_dir)
    for f in files:
        if f.endswith('.h5'):
            converted_files.append(f)
        else:
            new_file = os.path.join(base_dir, os.path.basename(f)+'.twiss')
            cmd = 'sddsanalyzebeam %s %s' % (f, new_file)
            print(cmd)
            os.system(cmd)
            new_h5_file = new_file + '.h5'
            cmd = 'sdds2hdf %s %s' % (new_file, new_h5_file)
            print(cmd)
            os.system(cmd)
            converted_files.append(new_h5_file)

    w_list = [watcher.FileViewer(f) for f in converted_files]

    w0 = w_list[0]

    for col in w0.columns:
        str_list = [col]
        print_ = True
        for w in w_list:
            if col not in w.columns:
                print_ = False
                print('%s not found in %s' % (col, w.filename))
                break
            str_list.append(str(w[col].squeeze()))
        if print_:
            print(' '.join(str_list))

    shutil.rmtree(base_dir)



