import os
#import functools
import itertools
import h5py
import numpy as np
import scipy.optimize as optimize
from scipy.constants import c, m_e, e

m_e_eV = m_e*c**2/e


def sdds2hdf(file_in, file_out):
    raw_file = file_in
    processed_file = file_out

    if not os.path.isfile(raw_file):
        try:
            import simulation
        except:
            from . import simulation
        raise simulation.ElegantWrapperError('File %s does not exist!' % raw_file)

    if not os.path.isfile(processed_file) or \
            os.path.getmtime(processed_file) < os.path.getmtime(raw_file):
        cmd = 'sdds2hdf %s %s 2>/dev/null 1>/dev/null' % (raw_file, processed_file)
        status = os.system(cmd)
        if status != 0:
            raise SystemError('Status %i for file %s\nCommand: %s' % (status, raw_file, cmd))

        #print('Generated %s' % processed_file)

def mu_fit_func(z, *parameters):
    output = 0
    for order, par in enumerate(parameters):
        zz = z**(order+1)
        output += par*(zz - np.mean(zz))
    return output

class FileViewer:
    def __init__(self, filename, sim=None, no_page1=False):
        if not filename.endswith('.h5'):
            raise ValueError('File is not an .h5')
        self.filename = filename
        self.sim = sim
        self._dict = {}
        self.no_page1 = no_page1

        if not os.path.isfile(self.filename):
            raise OSError('File %s not found!' % self.filename)
        with h5py.File(self.filename, 'r') as ff:
            if no_page1:
                self.columns = list(ff.keys())
                self.parameters = []
                self.pages = self.columns
            else:
                self.pages = list(ff.keys())
                keys0 = ff['page1'].keys()
                if 'columns' in keys0:
                    self.columns = list(ff['page1/columns'].keys())
                else:
                    self.columns = []
                if 'parameters' in keys0:
                    self.parameters = list(ff['page1/parameters'].keys())
                else:
                    self.parameters = []

    def __lt__(self, other):
        return self.filename < other.filename

    def __le__(self, other):
        return self.filename <= other.filename

    def __eq__(self, other):
        return hasattr(other, 'filename') and self.filename == other.filename

    def __ne__(self, other):
        return hasattr(other, 'filename') and self.filename != other.filename

    def __gt__(self, other):
        return self.filename > other.filename

    def __ge__(self, other):
        return self.filename >= other.filename

    def __str__(self):
        return self.__class__.__name__+': '+os.path.basename(self.filename)
    __repr__ = __str__

    def __getitem__(self, key, page=1):

        if self.no_page1:
            realkey = key
            if realkey not in self._dict:
                with h5py.File(self.filename, 'r') as ff:
                    try:
                        out = np.array(ff[realkey])
                    except KeyError:
                        raise KeyError('Key %s not found' % realkey)
                    if out.dtype == 'O':
                        out = np.array(out, dtype='U')
                    self._dict[realkey] = out
            return self._dict[realkey]
        else:
            if key in self.columns and key not in self.parameters:
                realkey = 'columns/' + key
            elif key in self.parameters and key not in self.columns:
                realkey = 'parameters/' + key
            else:
                realkey = key

            if realkey not in self._dict:
                with h5py.File(self.filename, 'r') as ff:
                    out = np.array(ff['page%i/' % page + realkey])
                    if out.dtype == 'O':
                        out = np.array(out, dtype='U')
                    self._dict[realkey] = out
            return self._dict[realkey]

    def print_tree(self):
        def name_and_size(key, ff):
            try:
                print((key, ff[key].shape, ff[key].dtype))
            except:
                print(key)

        with h5py.File(self.filename, 'r') as ff:
            ff.visit(lambda x: name_and_size(x, ff))

    def get_row(self, key, value):
        print('\t'.join(self.columns))

        indices = np.argwhere(self[key] == value)
        if len(indices) == 0:
            print('No results for %s in %s' % (key, self.filename))
            return

        for i_ctr, index in enumerate(indices):
            output = []
            for column in self.columns:
                output.append(str(self[column][index].squeeze()))
            print('\t'.join(output))


class Watcher(FileViewer):
    def __init__(self, filename, *args, **kwargs):
        FileViewer.__init__(self, filename, *args, **kwargs)
        self.s = kwargs.pop('s', None)

        zz_0 = -self['t']*c
        self.zz = zz_0 - np.mean(zz_0)
        if 'z' not in self._dict:
            self._dict['z'] = self.zz
        try:
            self.s = self['parameters/s'].squeeze()
        except KeyError:
            pass

    # Look at paper from Guetg et al: PRSTAB 18, 2015, p. 030701
    # to understand what mu, zij, mn, mnstar means

    def get_zz(self):
        zz_0 = self['columns/t']*c
        self.zz = zz_0 - np.mean(zz_0)
        if 'z' not in self._dict:
            self._dict['z'] = self.zz
        return self.zz

    #@functools.lru_cache()
    def get_mu(self, dimension, order=1):
        assert dimension in ('x', 'y', 'xp', 'yp')

        trans_0 = self['columns/'+dimension]
        # <x> = <x'> = 0
        trans = trans_0 - np.mean(trans_0)

        fit, _ = optimize.curve_fit(mu_fit_func, self.zz, trans, [1]*order)
        return fit

    #@functools.lru_cache()
    def get_zij(self, i, j):
        return np.mean(self.zz**(i+j)) - np.mean(self.zz**i)*np.mean(self.zz**j)

    #@functools.lru_cache()
    def get_mn(self, dimension, n, twi0, sig0_or_e0, debug=False):
        assert dimension in ('x', 'y')

        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]

        try:
            e0 = float(sig0_or_e0)
        except:
            e0 = sig0_or_e0['columns/e%s' % dimension][index]
        b0 = twi0['columns/beta%s' % dimension][index]
        a0 = twi0['columns/alpha%s' % dimension][index]
        g0 = (1. + a0**2)/b0

        output = 0.
        mu_vect = self.get_mu(dimension, n)
        mup_vect = self.get_mu(dimension+'p', n)
        mu = {i: mu_vect[i-1] for i in range(1, n+1)}
        mup = {i: mup_vect[i-1] for i in range(1, n+1)}

        for i,j in itertools.product(range(1,n+1), repeat=2):
            contribution = (b0*mup[i]*mup[j] + g0*mu[i]*mu[j] + 2*a0*mu[i]*mup[j]) * self.get_zij(i,j)
            output += contribution

            if debug:
                import pdb; pdb.set_trace()
        return output/e0

    def get_m1tilde(self, dimension, twi0, sig0):
        assert dimension in ('x', 'y')
        mu = self.get_mu(dimension, 1)
        mup = self.get_mu(dimension+'p', 1)
        etilde = self.get_etilde(dimension, 1, twi0, sig0)
        zz = self.get_zij(1,1)
        betatilde = self.get_beta_tilde(dimension, 1, twi0, sig0)
        alphatilde = self.get_alpha_tilde(dimension, 1, twi0, sig0)
        gammatilde = (1. + alphatilde**2)/betatilde
        return float((gammatilde*mu**2 + betatilde*mup**2 + 2*alphatilde*mu*mup)/etilde*zz)

    def get_mnstar(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')
        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]
        e0 = sig0['columns/e%s' % dimension][index]

        output = 0.
        mu_vect = self.get_mu(dimension, n)
        mup_vect = self.get_mu(dimension+'p', n)
        mu = {i: mu_vect[i-1] for i in range(1, n+1)}
        mup = {i: mup_vect[i-1] for i in range(1, n+1)}

        for i,j,k,l in itertools.product(range(1,n+1), repeat=4):
            tmp = mu[i]*mu[j]*mup[k]*mup[l] - mu[i]*mup[j]*mu[k]*mup[l]
            output += tmp * self.get_zij(i,j) * self.get_zij(k,l)

        return output/e0**2

    def get_etilde(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')
        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]
        e0 = sig0['columns/e%s' % dimension][index]

        mn = self.get_mn(dimension, n, twi0, sig0)
        mnstar = self.get_mnstar(dimension, n, twi0, sig0)

        return e0 * np.sqrt(1+mn+mnstar)

    def get_beta_tilde(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')
        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]
        e0 = sig0['columns/e%s' % dimension][index]
        b0 = twi0['columns/beta%s' % dimension][index]

        factor = 0
        mu = self.get_mu(dimension, n)
        for i,j in itertools.product(range(1,n+1), repeat=2):
            factor += mu[i-1]*mu[j-1]*self.get_zij(i,j)

        return (b0*e0 + factor) / self.get_etilde(dimension, n, twi0, sig0)

    def get_gamma_tilde(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')
        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]
        e0 = sig0['columns/e%s' % dimension][index]
        try:
            g0 = twi0['columns/gamma%s' % dimension][index]
        except:
            twi0.print_tree()
            raise

        factor = 0
        mup = self.get_mu(dimension+'p', n)
        for i,j in itertools.product(range(1,n+1), repeat=2):
            factor += mup[i-1]*mup[j-1]*self.get_zij(i,j)

        return (g0*e0 + factor) / self.get_etilde(dimension, n, twi0, sig0)

    def get_alpha_tilde(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')
        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]
        e0 = sig0['columns/e%s' % dimension][index]
        a0 = twi0['columns/alpha%s' % dimension][index]

        factor = 0
        mu = self.get_mu(dimension, n)
        mup = self.get_mu(dimension+'p', n)
        for i,j in itertools.product(range(1,n+1), repeat=2):
            factor += mu[i-1]*mup[j-1]*self.get_zij(i,j)

        return (a0*e0 - factor) / self.get_etilde(dimension, n, twi0, sig0)

    @staticmethod
    def get_emittance_from_points(x, xp):
        x = x - np.mean(x)
        xp = xp - np.mean(xp)
        return np.sqrt(np.mean(x**2)*np.mean(xp**2) - np.mean(x*xp)**2)

    def get_emittance_from_beam(self, dimension, normalized=False):
        assert dimension in ('x', 'y')

        space = self['%s' % dimension]
        angle = self['%sp' % dimension]
        if normalized:
            norm_factor = np.mean(self['p'])
        else:
            norm_factor = 1
        return self.get_emittance_from_points(space, angle)*norm_factor

    def get_beta_from_beam(self, dimension):
        assert dimension in ('x', 'y')
        e = self.get_emittance_from_beam(dimension)
        beamsize_squared = self.get_beamsize(dimension)**2
        return beamsize_squared/e

    def get_gamma_from_beam(self, dimension):
        assert dimension in ('x', 'y')
        e = self.get_emittance_from_beam(dimension)
        beamsize_squared = self.get_beamsize(dimension+'p')**2
        return beamsize_squared/e

    def get_alpha_from_beam(self, dimension):
        assert dimension in ('x', 'y')
        e = self.get_emittance_from_beam(dimension)
        arr_x0 = self[dimension]
        arr_xp0 = self[dimension+'p']

        arr_x = arr_x0 - arr_x0.mean()
        arr_xp = arr_xp0 - arr_xp0.mean()

        return -np.mean(arr_x*arr_xp)/e

    def get_optics_from_beam(self, dimension):
        """
        Returns beta, alpha, gamma
        """
        return (
                self.get_beta_from_beam(dimension),
                self.get_alpha_from_beam(dimension),
                self.get_gamma_from_beam(dimension),
                )

    def get_mismatch(self, dimension, beta_ref, alpha_ref):
        beta, alpha, gamma = self.get_optics_from_beam(dimension)
        gamma_ref = (1+alpha_ref**2)/beta_ref
        mismatch = (beta*gamma_ref - 2*alpha*alpha_ref + gamma*beta_ref)/2.
        #import pdb; pdb.set_trace()
        return mismatch

    def get_current(self, t_or_z, charge=None, bins=50, center_time=True):
        """
        Output: bin_edges, bin_values. Use with with plt.step.
        """
        if t_or_z == 't':
            zz = self['t']
        elif t_or_z == 'z':
            zz = -self['t']*c
        else:
            raise ValueError('t_or_z must be either t or z, is:', t_or_z)
        if center_time:
            zz -= zz.mean()
        #print('Max:', zz.max())
        hist, bin_edges = np.histogram(zz, bins=bins)
        #sp.step(bin_edges[:-1], hist)
        #sp.hist(zz, bins=50)
        if charge is None:
            factor = 1
        elif t_or_z == 't':
            factor = charge/(np.diff(bin_edges)[0]*np.sum(hist))
        elif t_or_z == 'z':
            factor = c*charge/(np.diff(bin_edges)[0]*np.sum(hist))
        return bin_edges, np.concatenate(([0,],hist*factor))

    def get_beamsize(self, dimension):
        assert dimension in ('x', 'y', 'xp', 'yp', 't', 'p')
        arr = self[dimension]-np.mean(self[dimension])
        return np.sqrt(np.mean(arr**2))

    def get_mean(self, dimension):
        assert dimension in ('x', 'y', 'xp', 'yp', 't', 'p')
        return np.mean(self[dimension])

    def gaussianBeamsizeFit(self, dimension):
        assert dimension in ('x', 'y', 'xp', 'yp')

        def gauss_function(x, a, x0, sigma):
            return a*np.exp(-(x-x0)**2/(2*sigma**2))

        xx = self[dimension]
        hist, bin_edges = np.histogram(xx, bins=50)
        xx = bin_edges[:-1]+np.diff(bin_edges).mean()/2

        mean = np.mean(xx)
        sigma = self.get_beamsize(dimension)/np.sqrt(2)

        popt, pcov = optimize.curve_fit(gauss_function, xx, hist, p0=[hist.max(), mean, sigma])
        #from visualize_fit import visualize_fit
        #visualize_fit(gauss_function, xx, hist, popt)

        return popt[-1]

    def slice_beam(self, n_bins, axis='z', method='const_delta'):
        """
        Method may be const_delta or const_size
        """
        if axis == 'z':
            xx = -self['t']*c
        else:
            xx = self[axis]
        xx = xx - xx.mean()

        if method == 'const_delta':
            slices = []
            xx_min, xx_max = np.min(xx), np.max(xx)
            delta = (xx_max-xx_min)/n_bins
            for n_bin in range(n_bins):
                if n_bin == n_bins-1:
                    this_mask = np.logical_and(xx_min+n_bin*delta <= xx, xx_min+(n_bin+1)*delta >= xx)
                else:
                    this_mask = np.logical_and(xx_min+n_bin*delta <= xx, xx_min+(n_bin+1)*delta > xx)
                columns_dict = {column: self[column][this_mask] for column in self.columns}
                parameters_dict = {parameter: self[parameter] for parameter in self.parameters}
                new_slice = Watcher2(parameters_dict, columns_dict, s=self.s)
                slices.append(new_slice)
                new_slice.delta = delta

            return slices

        elif method == 'const_size':
            slices = []
            n_slice = len(xx)//n_bins
            #sort_arr = np.sort(np.array([xx, np.arange(len(xx))]).T, axis=0)
            sort_arr = np.array(list(zip(xx, np.arange(len(xx)))), dtype=[('value', float), ('counter', int)])
            sorted_arr = np.sort(sort_arr, order='value')
            indices = sorted_arr['counter']

            for n_bin in range(n_bins):
                this_indices = indices[n_bin*n_slice:(n_bin+1)*n_slice]
                columns_dict = {column: np.take(self[column], this_indices) for column in self.columns}
                parameters_dict = {parameter: self[parameter] for parameter in self.parameters}
                new_slice = Watcher2(parameters_dict, columns_dict, s=self.s)
                slices.append(new_slice)

            return slices

    def toGenesis(self, n_slices, filename=None):
        if filename is None:
            filename = self.filename.replace('.h5','')
        new_name = filename+'.genesis%i' % n_slices
        if not os.path.isfile(new_name) or \
                os.path.getmtime(new_name) < os.path.getmtime(filename):

            cmd = 'elegant2genesis %s %s -slices=%i' % (filename, new_name, n_slices)
            print(cmd)
            os.system(cmd)
        new_name_h5 = new_name+'.h5'
        sdds2hdf(new_name, new_name_h5)
        return FileViewer(new_name_h5)

    def toSDDS(self, filename):
        with open(filename, 'w') as fid:
            fid.write('SDDS1\n')
            #fid.write('&column name=z,    units=m,    type=double,    &end\n')
            fid.write('&column name=t,    units=s,    type=double,    &end\n')
            fid.write('&column name=x,    units=m,  type=double,    &end\n')
            fid.write('&column name=y,    units=m,  type=double,    &end\n')
            fid.write('&column name=xp,  type=double,    &end\n')
            fid.write('&column name=yp,  type=double,    &end\n')
            fid.write('&column name=p,    units=m$be$nc,  type=double,    &end\n')
            fid.write('&data mode=ascii, &end\n')
            #fid.write('! page number 1\n')
            fid.write('%i\n' % len(self['t']))
            for t, x, y, xp, yp, p in zip(
                    self['t'], self['x'], self['y'], self['xp'], self['yp'], self['p']):
                fid.write('  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e  %12.6e\n' % (t, x, y, xp, yp, p))


class Watcher2(Watcher):
    """
    Pseudo Watcher that returns values from a dictionary but otherweise behaves like a Watcher.
    """

    def __init__(self, parameters_dict, columns_dict, s=None, filename='None'):
        self.s = s
        self.filename = filename
        self.columns = list(columns_dict.keys())
        self.parameters = list(parameters_dict.keys())
        self.columns_dict = columns_dict
        self.parameters_dict = parameters_dict


        zz_0 = self['t']*c
        self.zz = zz_0 - np.mean(zz_0)

    def __getitem__(self, key):
        if key in self.columns and key in self.parameters:
            raise ValueError('Key %s in both columns and parameters!' % key)
        elif key in self.columns:
            return self.columns_dict[key]
        elif key in self.parameters:
            return self.parameters_dict[key]
        elif key.startswith('columns/'):
            return self.columns_dict[key[8:]]
        elif key.startswith('parameters/'):
            return self.parameters_dict[key[11:]]
        else:
            raise KeyError(key)

class SliceCollection:
    def __init__(self, slices, parent):
        self.slices = slices
        self.n_slices = len(slices)
        self.parent = parent

        s_list = []

        mean_t = parent['t'].mean()

        for s in slices:
            s_list.append(-(np.mean(s['t'])-mean_t)*c)

        self.s_arr = np.array(s_list)

    def get_mismatch(self, dimension, reference):
        if reference == 'proj':
            reference = self.parent
        if hasattr(reference, '__len__') and len(reference) == 2:
            beta_ref, alpha_ref = reference
        elif hasattr(reference, 'get_beta_from_beam'):
            beta_ref = reference.get_beta_from_beam(dimension)
            alpha_ref = reference.get_alpha_from_beam(dimension)
        else:
            raise ValueError(reference)
        gamma_ref = (1+alpha_ref**2)/beta_ref

        mm_list = []
        z_list = []
        mean_t = self.parent['t'].mean()
        for s in self.slices:
            z_list.append(-(s['t'].mean()-mean_t)*c)
            alpha = s.get_alpha_from_beam(dimension)
            beta = s.get_beta_from_beam(dimension)
            gamma = s.get_gamma_from_beam(dimension)
            mismatch = (beta*gamma_ref - 2*alpha*alpha_ref + gamma*beta_ref)/2.
            mm_list.append(mismatch)
        return np.array(z_list), np.array(mm_list)

    def get_misplacement(self, dimension, reference_optics, reference_orbit):
        if hasattr(reference_optics, '__len__') and len(reference_optics) == 3:
            beta_ref, alpha_ref, eps_ref = reference_optics
        else:
            if type(reference_optics) is str and reference_optics == 'proj':
                reference = self.parent
            else:
                reference = reference_optics
            beta_ref = reference.get_beta_from_beam(dimension)
            alpha_ref = reference.get_alpha_from_beam(dimension)
            eps_ref = reference.get_emittance_from_beam(dimension)
        gamma_ref = (1+alpha_ref**2)/beta_ref


        if reference_orbit == 'ref':
            x_ref = reference[dimension].mean()
            xp_ref = reference[dimension+'p'].mean()
        else:
            x_ref, xp_ref = reference_orbit

        m = []
        for s in self.slices:
            x = np.mean(s[dimension]) - x_ref
            xp = np.mean(s[dimension+'p']) - xp_ref

            inv = (beta_ref*xp**2 + 2*alpha_ref*xp*x + gamma_ref*x**2)/eps_ref
            m.append(inv)

        return np.array(m)

    def get_slice_func(self, funcname, *args, **kwargs):
        output = []
        for s in self.slices:
            func = getattr(s, funcname)
            output.append(func(*args, **kwargs))
        return np.array(output)


def watcher_to_parray(w, charge):
    nn = w['x'].size
    rparticles = np.zeros([6, nn])
    rparticles[0] = w['x']
    rparticles[1] = w['xp']
    rparticles[2] = w['y']
    rparticles[3] = w['yp']
    rparticles[4] = w['t']*c
    rparticles[5] = w['p'] / w['p'].mean() - 1
    E = w['p'].mean()/1e9
    q_array = np.ones(nn, float) * charge / nn
    return rparticles, E, q_array

def parray_to_watcher(pa):
    dd = {}
    rparticles = pa.rparticles
    energy_eV = pa.E*1e9
    dd['x'] = rparticles[0].copy()
    dd['xp'] = rparticles[1].copy()
    dd['y'] = rparticles[2].copy()
    dd['yp'] = rparticles[3].copy()
    dd['t'] = rparticles[4].copy()/c
    dd['p'] = energy_eV * (1 + rparticles[5].copy())
    return Watcher2({}, dd)

