import os
import pdb

# Generic
def copy_and_replace_file(dirname, filename, replace, call_pdb=False):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            full_input = ''.join(lines)
    except FileNotFoundError:
        print(filename, 'not found')
        if call_pdb:
            pdb.set_trace()

    for keyword in replace:
        if keyword not in full_input:
            print('%s not part of %s!' % (keyword, filename))
            if call_pdb:
                pdb.set_trace()
            raise ValueError

    if replace:
        new_lines = []
        for line in lines:
            for keyword, replacement in replace.items():
                line = line.replace(keyword, replacement)
            new_lines.append(line)
    else:
        new_lines = lines

        for line in new_lines:
            if '__' in line:
                print('Not all converted:\n', line)
                if call_pdb:
                    pdb.set_trace()
                raise ValueError

    try:
        filename_out = os.path.join(dirname, filename)
        with open(filename_out, 'w') as f:
            f.writelines(new_lines)
    except FileNotFoundError:
        print(filename_out, 'not found')
        if call_pdb:
            pdb.set_trace()
        raise ValueError

def split_runfile(file_template, bash_lines, n_files):
    run_of_runs = []
    outputs = []
    for n_file in range(n_files):
        outputs.append([])
    for ctr, line in enumerate(bash_lines):
        n_file = ctr % n_files
        outputs[n_file].append(line)
    for n_file in range(n_files):
        filename = file_template % (n_file+1)
        with open(filename, 'w') as f:
            f.writelines(outputs[n_file])
            print('Wrote %s' % filename)
        run_of_runs.append('bash %s &' % filename)

    filename = file_template % 0
    with open(filename, 'w') as f:
        f.writelines(run_of_runs)
        print('Wrote %s' % filename)

