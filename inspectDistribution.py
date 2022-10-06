import sys
import numpy as np
import matplotlib.pyplot as plt

from GenesisWrapper import match_particle_dist
try:
    import watcher
    import myplotstyle as ms
except ImportError:
    from . import watcher
    from . import myplotstyle as ms

def inspect(watch, bins=(100,100), show=True, title=None, charge=200e-12, center_time=True):
    dimensions = 'x', 'xp', 'y', 'yp', 't', 'p'
    units = '$\mu$m', '$\mu$rad', '$\mu$m', '$\mu$rad', 'fs', 'GeV'
    factors = 1e6, 1e6, 1e6, 1e6, 1e15, 511e3/1e9

    if title is None:
        title = str(watch)
    fig = ms.figure(title)
    plt.subplots_adjust(hspace=0.35, wspace=0.25, bottom=0.1)
    subplot = ms.subplot_factory(4,5,False)
    sp_ctr = 1

    for dim_ctr1, dim1 in enumerate(dimensions):
        for dim_ctr2, dim2 in enumerate(dimensions[dim_ctr1+1:], dim_ctr1+1):
            xlabel = '%s (%s)' % (dim1, units[dim_ctr1])
            ylabel = '%s (%s)' % (dim2, units[dim_ctr2])
            xfactor = factors[dim_ctr1]
            yfactor = factors[dim_ctr2]

            sp = subplot(sp_ctr, title='%s-%s' % (dim1, dim2), xlabel=xlabel, ylabel=ylabel)
            sp_ctr += 1
            x_arr = watch[dim1]
            y_arr = watch[dim2]
            if center_time:
                if dim1 == 't':
                    x_arr = x_arr - x_arr.mean()
                if dim2 == 't':
                    y_arr = y_arr - y_arr.mean()
            sp.hist2d(x_arr*xfactor, y_arr*yfactor, bins=bins)

    sp = subplot(sp_ctr, title='Beam current', xlabel='t (fs)', ylabel='I (kA)')
    sp_ctr += 1
    xx, yy = watch.get_current('t', bins=bins[0], charge=charge, center_time=center_time)
    sp.step(xx*1e15, yy/1e3)

    sp = subplot(sp_ctr, title='Norm. emittance', xlabel='t (fs)', ylabel='$\epsilon_n$ (nm)')
    sp_ctr += 1

    slices = watcher.SliceCollection(watch.slice_beam(bins[0]), watch)
    tt = np.array([s['t'].mean() for s in slices.slices])
    for dim in ('x', 'y'):
        em = slices.get_slice_func('get_emittance_from_beam', dim, True)
        sp.plot(tt*1e15, em*1e9, label=dim)
    sp.legend()

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
        watch = watcher.Watcher(new_in)
    except:
        dist = match_particle_dist.h5_in_genesis(new_in)
        watch = watcher.Watcher2({}, dist)
    inspect(watch)

