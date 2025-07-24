import os
import numpy as np
from scipy.constants import m_e, e, c
#import matplotlib.pyplot as plt

from PassiveWFMeasurement import image_analysis
import PassiveWFMeasurement.myplotstyle as ms

m_e_eV = m_e*c**2/e

def plot(sim, watchplot='long', void_cut=(2e-3, 2e-3)):
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
    sp_r56.set_ylabel('R56 (mm)')

    for dim in ['x', 'y']:
        color = sp_beta.plot(sim.twi['s'], sim.twi['beta'+dim], label=dim)[0].get_color()
        sp_beta.plot(sim.sig['s'], sim.sig['beta'+dim+'Beam'], ls='--', color=color)
        sp_disp.plot(sim.twi['s'], sim.twi['eta'+dim], label=dim)
    sp_beta.legend()
    sp_disp.legend()

    sp_ene.plot(sim.twi['s'], sim.twi['pCentral0']/1e6*m_e_eV, label='twi')
    sp_ene.plot(sim.cen['s'], sim.cen['pCentral']/1e6*m_e_eV, label='cen')
    sp_ene.legend()

    if np.any(sim.sig['s7'] > 0):
        sp_blen.semilogy(sim.sig['s'], sim.sig['s7']*1e15, label='s7')
    if sim.mat:
        sp_r56.plot(sim.mat['s'], sim.mat['R56']*1e3, color='tab:orange', label='mat')
    ms.comb_legend(sp_blen, sp_r56)

    n_watch_plots = 10 - sp_ctr

    if watchplot == 'long':
        xlabel = '$t$ (fs)'
        ylabel = '$\Delta E$ (MeV)'
        xdim = 't'
        ydim = 'p'
        xfactor = 1
        yfactor = m_e_eV
        x_unit = 's'
        x_unit2 = 'fs'
        y_unit = 'eV'
        y_unit2 = 'MeV'
    elif watchplot == 'trans':
        xlabel = '$x$ (mm)'
        ylabel = '$y$ (mm)'
        xdim = 'x'
        ydim = 'y'
        xfactor = 1
        yfactor = 1
        x_unit = y_unit = 'm'
        x_unit2 = y_unit2 = 'mm'
    else:
        raise ValueError(watchplot)

    for w in sim.watch[:n_watch_plots]:
        if len(w['x']) == 1:
            continue
        sp_w = subplot(sp_ctr, title='%s at %.1f m' % (os.path.basename(w.filename), w.s), xlabel=xlabel, ylabel=ylabel, grid=False)
        sp_ctr += 1

        xx0, yy0 = w[xdim], w[ydim]
        hist, xedges, yedges = np.histogram2d(xx0, yy0, bins=(100, 100))
        x_axis = xedges[:-1] + xedges[1] - xedges[0]
        y_axis= yedges[:-1] + yedges[1] - yedges[0]
        img0 = image_analysis.Image(hist.T, x_axis, y_axis)
        img0 = img0.cut_voids(*void_cut)

        mask_x = np.logical_and(xx0 > img0.x_axis.min(), xx0 < img0.x_axis.max())
        mask_y = np.logical_and(yy0 > img0.y_axis.min(), yy0 < img0.y_axis.max())
        mask = np.logical_and(mask_x, mask_y)

        xx1 = xx0[mask]
        yy1 = yy0[mask]
        xx = (xx1 - xx1.mean())*xfactor
        yy = (yy1 - yy1.mean())*yfactor
        hist, xedges, yedges = np.histogram2d(xx, yy, bins=(100, 100))
        x_axis = xedges[:-1] + xedges[1] - xedges[0]
        y_axis= yedges[:-1] + yedges[1] - yedges[0]
        img = image_analysis.Image(hist.T, x_axis, y_axis, charge=w['Charge'], energy_eV=w['p'].mean()*m_e_eV, x_unit=x_unit, y_unit=y_unit)
        img = img.cut_voids(*void_cut)

        img.plot_img_and_proj(sp_w, plot_gauss=False)

        textbbox = {'boxstyle': 'square', 'alpha': 0.75, 'facecolor': 'white', 'edgecolor': 'gray'}
        x_factor = image_analysis.unit_to_factor(img.x_unit)
        y_factor = image_analysis.unit_to_factor(img.y_unit)
        textstr = '%s: %.2f %s; %s: %.2f %s (rms)' % (xdim, xx.std()*x_factor, x_unit2, ydim, yy.std()*y_factor, y_unit2)
        sp_w.text(0.95, 0.05, textstr, transform=sp_w.transAxes, verticalalignment='bottom', horizontalalignment='right', bbox=textbbox)


        #extent = [x_axis[0]*xfactor, x_axis[-1]*xfactor, y_axis[0]*yfactor, y_axis[-1]*yfactor]
        #sp_w.imshow(hist, aspect='auto', extent=extent, origin='lower', cmap=plt.get_cmap('hot'))


