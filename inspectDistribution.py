import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c, m_e, e

m_e_eV = m_e*c**2/e

from GenesisWrapper import match_particle_dist
try:
    import watcher
    import myplotstyle as ms
except ImportError:
    from . import watcher
    from . import myplotstyle as ms

def inspect(watch, bins=(100,100), show=True, title=None, charge=200e-12, center_time=True, fig_kwargs={}, hspace=0.35, wspace=0.25, bottom=0.1):
    dimensions = 'x', 'xp', 'y', 'yp', 't', 'p'
    units = '$\mu$m', '$\mu$rad', '$\mu$m', '$\mu$rad', 'fs', 'GeV/c'
    factors = 1e6, 1e6, 1e6, 1e6, 1e15, m_e_eV/1e9

    if title is None:
        title = str(watch)
    fig = ms.figure(title, **fig_kwargs)
    fig.set_tight_layout(True)
    plt.subplots_adjust(hspace=hspace, wspace=wspace, bottom=bottom)
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
    curr_time, curr = watch.get_current('t', bins=bins[0], charge=charge, center_time=center_time)
    sp.step(curr_time*1e15, curr/1e3)

    sp = subplot(sp_ctr, title='Norm. emittance', xlabel='t (fs)', ylabel='$\epsilon_n$ (nm)')
    sp_ctr += 1
    sp_curr = sp.twinx()
    #sp_curr.set_ylabel('I (kA)')
    sp_curr.step(curr_time*1e15, curr/1e3, color='black')
    sp_curr.set_yticks([])

    sp_beta = subplot(sp_ctr, title='Slice optics', xlabel='t (fs)', ylabel=r'$\beta$ (m)')
    sp_ctr += 1
    sp_alpha = sp_beta.twinx()
    sp_alpha.set_ylabel(r'$\alpha$')

    slices = watcher.SliceCollection(watch.slice_beam(bins[0], method='const_size'), watch)
    tt = np.array([s['t'].mean() for s in slices.slices])
    for dim in ('x', 'y'):
        em = slices.get_slice_func('get_emittance_from_beam', dim, True)
        beta = slices.get_slice_func('get_beta_from_beam', dim)
        alpha = slices.get_slice_func('get_alpha_from_beam', dim)
        sp.plot(tt*1e15, em*1e9, label=dim)
        color = sp_beta.plot(tt*1e15, beta, label=dim)[0].get_color()
        sp_alpha.plot(tt*1e15, alpha, ls='--', color=color)

        beta = watch.get_beta_from_beam(dim)
        alpha = watch.get_alpha_from_beam(dim)
        color = {'x': 'black', 'y': 'gray'}[dim]
        sp_beta.axhline(beta, color=color, ls='solid')
        sp_alpha.axhline(alpha, color=color, ls='dotted')
    sp.legend()
    sp_beta.legend()

    sp = subplot(sp_ctr, title='Energy spread', xlabel='t (fs)', ylabel='$\sigma_E$ (MeV)')
    sp_ctr += 1
    sp_curr = sp.twinx()
    #sp_curr.set_ylabel('I (kA)')
    sp_curr.step(curr_time*1e15, curr/1e3, color='black')
    sp_curr.set_yticks([])

    pspread = slices.get_slice_func('get_beamsize', 'p')
    espread = pspread*m_e_eV
    sp.plot(tt*1e15, espread/1e6)

    sp = subplot(sp_ctr, title='Energy chirp', xlabel='t (fs)', ylabel='Chirp (MeV/fs)')
    sp_ctr += 1
    sp_curr = sp.twinx()
    #sp_curr.set_ylabel('I (kA)')
    sp_curr.step(curr_time*1e15, curr/1e3, color='black')
    sp_curr.set_yticks([])


    pmean = slices.get_slice_func('get_mean', 'p')
    emean = pmean*m_e_eV
    chirp = np.diff(emean)/np.diff(tt)
    tt_plot = tt[:-1] + (tt[1] - tt[0])/2
    sp.plot(tt_plot*1e15, chirp/1e6/1e15)

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

