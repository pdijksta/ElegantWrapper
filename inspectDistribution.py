import sys
import matplotlib.pyplot as plt
import myplotstyle as ms

import GenesisWrapper.match_particle_dist as match_particle_dist

try:
    import watcher
except ImportError:
    from . import watcher

def inspect(watcher, bins=(100,100), show=True, title=None, charge=200e-12):
    dimensions = 'x', 'xp', 'y', 'yp', 't', 'p'

    if title is None:
        title = watcher
    fig = ms.figure(title)
    plt.subplots_adjust(hspace=0.35, wspace=0.25, bottom=0.1)
    subplot = ms.subplot_factory(4,4)
    sp_ctr = 1

    for dim_ctr1, dim1 in enumerate(dimensions):
        for dim_ctr2, dim2 in enumerate(dimensions[dim_ctr1+1:], dim_ctr1+1):
            sp = subplot(sp_ctr, title='%s-%s' % (dim1, dim2), xlabel=dim1, ylabel=dim2, scix=True, sciy=True, grid=False)
            sp_ctr += 1
            x_arr = watcher[dim1]
            y_arr = watcher[dim2]
            if dim1 == 't':
                x_arr = x_arr - x_arr.mean()
            if dim2 == 't':
                y_arr = y_arr - y_arr.mean()
            sp.hist2d(x_arr, y_arr, bins=bins)

    sp = subplot(sp_ctr, title='Beam current', xlabel='t', ylabel='I [A]', scix=True, sciy=True, grid=False)
    xx, yy = watcher.get_current('t', bins=bins[0], charge=charge)
    sp.step(xx, yy)

    if show:
        plt.show()

    return fig


if __name__ == '__main__':
    in_ = sys.argv[1]
    if not in_.endswith('.h5'):
        new_in = in_ + '.h5'
        watcher.sdds2hdf(in_, new_in)
    else:
        new_in = in_
    try:
        watcher = watcher.Watcher(new_in)
    except:
        dist = match_particle_dist.h5_in_genesis(new_in)
        watcher = watcher.Watcher2({}, dist)
    inspect(watcher)

