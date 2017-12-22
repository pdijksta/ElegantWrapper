import functools
import itertools
import h5py
import numpy as np
import scipy.optimize as optimize
from scipy.constants import c

class FileViewer:
    def __init__(self, filename, sim=None):
        if not filename.endswith('.h5'):
            raise ValueError('File is not an .h5')
        self.filename = filename
        self.sim = sim
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


class Watcher(FileViewer):
    def __init__(self, *args, **kwargs):
        if 's' in kwargs:
            self.s = kwargs['s']
            del kwargs['s']
        super().__init__(*args, **kwargs)

        zz_0 = self['columns/t']*c
        self.zz = zz_0 - np.mean(zz_0)
        try:
            self.s = self['parameters/s']
        except KeyError:
            pass

    # Look at paper from Guetg et al: PRSTAB 18, 2015 p. 030701
    # to understand what mu, zij, mn, mnstar means

    @functools.lru_cache()
    def get_mu(self, dimension, order=1):
        assert dimension in ('x', 'y', 'xp', 'yp')

        trans_0 = self['columns/'+dimension]
        if dimension[-1] == 'p':
            trans = trans_0
        else:
            trans = trans_0 - np.mean(trans_0)

        def fit_func(z, *parameters):
            output = 0
            for order, par in enumerate(parameters):
                output += par*z**(order+1)
            return output

        fit, _ = optimize.curve_fit(fit_func, self.zz, trans, [1]*order)

        return fit[-1]

    @functools.lru_cache()
    def get_zij(self, i, j):
        return np.mean(self.zz**(i+j)) - np.mean(self.zz**i)*np.mean(self.zz**j)

    @functools.lru_cache()
    def get_mn(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')

        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]

        e0 = sig0['columns/e%s' % dimension][index]
        b0 = twi0['columns/beta%s' % dimension][index]
        a0 = twi0['columns/alpha%s' % dimension][index]
        g0 = (1 + a0**2)/b0

        output = 0
        mu = {i: self.get_mu(dimension, i) for i in range(1, n+1)}
        mup = {i: self.get_mu(dimension+'p', i) for i in range(1, n+1)}

        for i,j in itertools.product(range(1,n+1), repeat=2):
            output += (b0*mup[i]*mup[j] + g0*mu[i]*mu[j] + 2*a0*mu[i]*mup[j]) * self.get_zij(i,j)

        return output/e0

    @functools.lru_cache()
    def get_mnstar(self, dimension, n, twi0, sig0):
        assert dimension in ('x', 'y')
        index = np.argwhere(twi0['columns/s'] == self.s)[0,0]
        e0 = sig0['columns/e%s' % dimension][index]

        output = 0
        mu = {i: self.get_mu(dimension, i) for i in range(1, n+1)}
        mup = {i: self.get_mu(dimension+'p', i) for i in range(1, n+1)}

        for i,j,k,l in itertools.product(range(1,n+1), repeat=4):
            tmp = (mu[i]*mu[j]*mup[k]*mup[l] - mu[i]*mup[j]*mu[k]*mup[l])
            output += tmp * self.get_zij(i,j) * self.get_zij(k,j)

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
        for i,j in itertools.product(range(1,n+1), repeat=2):
            mu_i = self.get_mu(dimension, i)
            mu_j = self.get_mu(dimension, j)
            factor += mu_i*mu_j*self.get_zij(i,j)

        return (b0*e0 + factor) / self.get_etilde(dimension, n, twi0, sig0)

