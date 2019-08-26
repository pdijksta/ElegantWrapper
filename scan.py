import os
# Generic
def copy_and_replace_file(dirname, filename, replace):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            full_input = ''.join(lines)
    except FileNotFoundError as e:
        print(filename, 'not found')
        import pdb; pdb.set_trace()

    for keyword in replace:
        if keyword not in full_input:
            raise ValueError('%s not part of %s!' % (keyword, filename))

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
                raise ValueError('Not all converted:\n', line)

    try:
        filename_out = os.path.join(dirname, filename)
        with open(filename_out, 'w') as f:
            f.writelines(new_lines)
    except FileNotFoundError as e:
        print(filename_out, 'not found')
        import pdb; pdb.set_trace()

