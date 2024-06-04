import os
import numpy as np
from scipy.constants import m_e, e, c
import matplotlib.pyplot as plt

import PassiveWFMeasurement.myplotstyle as ms

m_e_eV = m_e*c**2/e

def plot(sim, watchplot='long'):
    fig = ms.figure(sim.filename)
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    subplot = ms.subplot_factory(3, 3, True)
    sp_ctr = 1

    sp_beta = subplot(sp_ctr, title='Beta functions', xlabel='$s$ (m)', ylabel=r'$\beta$ (m)')
    sp_ctr += 1

    sp_ene = subplot(sp_ctr, title='Energy', xlabel='$s$ (m)', ylabel='$E$ (MeV)')
    sp_ctr += 1

    sp_disp = subplot(sp_ctr, title='Dispersion', xlabel='$s$ (m)', ylabel='$D$ (m)')
    sp_ctr += 1

    sp_blen = subplot(sp_ctr, title='Bunch duration', xlabel='$s$ (m)', ylabel='sig s7 (fs)')
    sp_ctr += 1
    sp_r56 = sp_blen.twinx()
    sp_r56.set_ylabel('R56 (m)')

    for dim in ['x', 'y']:
        color = sp_beta.plot(sim.twi['s'], sim.twi['beta'+dim], label=dim)[0].get_color()
        sp_beta.plot(sim.sig['s'], sim.sig['beta'+dim+'Beam'], ls='--', color=color)
        sp_disp.plot(sim.twi['s'], sim.twi['eta'+dim], label=dim)
    sp_beta.legend()
    sp_disp.legend()

    sp_ene.plot(sim.twi['s'], sim.twi['pCentral0']/1e6*m_e_eV, label='twi')
    sp_ene.plot(sim.cen['s'], sim.cen['pCentral']/1e6*m_e_eV, label='cen')
    sp_ene.legend()

    sp_blen.semilogy(sim.sig['s'], sim.sig['s7']*1e15)

    n_watch_plots = 10 - sp_ctr

    if watchplot == 'long':
        xlabel = '$t$ (fs)'
        ylabel = '$\Delta E$ (MeV)'
        xdim = 't'
        ydim = 'p'
        xfactor = 1e15
        yfactor = m_e_eV/1e6
    elif watchplot == 'trans':
        xlabel = '$x$ (mm)'
        ylabel = '$y$ (mm)'
        xdim = 'x'
        ydim = 'y'
        xfactor = 1e-3
        yfactor = 1e-3
    else:
        raise ValueError(watchplot)

    for w in sim.watch[:n_watch_plots]:
        if len(w['x']) == 1:
            continue
        sp_w = subplot(sp_ctr, title=os.path.basename(w.filename), xlabel=xlabel, ylabel=ylabel, grid=False)
        sp_ctr += 1

        xx = w[xdim] - w[xdim].mean()
        yy = w[ydim] - w[xdim].mean()
        hist, xedges, yedges = np.histogram2d(xx, yy, bins=(100, 100))
        x_axis = xedges[:-1] + xedges[1] - xedges[0]
        y_axis= yedges[:-1] + yedges[1] - yedges[0]
        extent = [x_axis[0]*xfactor, x_axis[-1]*xfactor, y_axis[0]*yfactor, y_axis[-1]*yfactor]
        sp_w.imshow(hist, aspect='auto', extent=extent, origin='lower', cmap=plt.get_cmap('hot'))


    if sim.mat:
        sp_r56.plot(sim.mat['s'], sim.mat['R56'])

