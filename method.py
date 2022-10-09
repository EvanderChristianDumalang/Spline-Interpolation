dfitpack_int = dfitpack.types.intvar.dtype

_extrap_modes = {0: 0, 'extrapolate': 0,
                 1: 1, 'zeros': 1,
                 2: 2, 'raise': 2,
                 3: 3, 'const': 3}

class UnivariateSpline(object):
    def __init__(self, x, y, w=None, bbox=[None]*2, k=3, s=None,
                 ext=0, check_finite=False):

        x, y, w, bbox, self.ext = self.validate_input(x, y, w, bbox, k, s, ext,
                                                      check_finite)
        data = dfitpack.fpcurf0(x, y, k, w=w, xb=bbox[0],
                                xe=bbox[1], s=s)
        if data[-1] == 1:
            data = self._reset_nest(data)
        self._data = data
        self._reset_class()

    @staticmethod
    def validate_input(x, y, w, bbox, k, s, ext, check_finite):
        x, y, bbox = np.asarray(x), np.asarray(y), np.asarray(bbox)
        if w is not None:
            w = np.asarray(w)
        if check_finite:
            w_finite = np.isfinite(w).all() if w is not None else True
            if (not np.isfinite(x).all() or not np.isfinite(y).all() or
                    not w_finite):
                raise ValueError("x and y array must not contain "
                                 "NaNs or infs.")
        if s is None or s > 0:
            if not np.all(diff(x) >= 0.0):
                raise ValueError("x must be increasing if s > 0")
        else:
            if not np.all(diff(x) > 0.0):
                raise ValueError("x must be strictly increasing if s = 0")
        if x.size != y.size:
            raise ValueError("x and y should have a same length")
        elif w is not None and not x.size == y.size == w.size:
            raise ValueError("x, y, and w should have a same length")
        elif bbox.shape != (2,):
            raise ValueError("bbox shape should be (2,)")
        elif not (1 <= k <= 5):
            raise ValueError("k should be 1 <= k <= 5")
        elif s is not None and not s >= 0.0:
            raise ValueError("s should be s >= 0.0")

        try:
            ext = _extrap_modes[ext]
        except KeyError:
            raise ValueError("Unknown extrapolation mode %s." % ext)

        return x, y, w, bbox, ext

    @classmethod
    def _from_tck(cls, tck, ext=0):
        self = cls.__new__(cls)
        t, c, k = tck
        self._eval_args = tck
        self._data = (None, None, None, None, None, k, None, len(t), t,
                      c, None, None, None, None)
        self.ext = ext
        return self

    def _reset_class(self):
        data = self._data
        n, t, c, k, ier = data[7], data[8], data[9], data[5], data[-1]
        self._eval_args = t[:n], c[:n], k
        if ier == 0:
            pass
        elif ier == -1:
            self._set_class(InterpolatedUnivariateSpline)

    def _set_class(self, cls):
        self._spline_class = cls
        if self.__class__ in (UnivariateSpline, InterpolatedUnivariateSpline):
            self.__class__ = cls
        else:
            pass

    def _reset_nest(self, data, nest=None):
        n = data[10]
        if nest is None:
            k, m = data[5], len(data[0])
            nest = m+k+1
        else:
            if not n <= nest:
                raise ValueError("`nest` can only be increased")
        t, c, fpint, nrdata = [np.resize(data[j], nest) for j in
                               [8, 9, 11, 12]]

        args = data[:8] + (t, c, n, fpint, nrdata, data[13])
        data = dfitpack.fpcurf1(*args)
        return data

    def set_smoothing_factor(self, s):
        data = self._data
        if data[6] == -1:
            warnings.warn('smoothing factor unchanged for'
                          'spline with fixed knots')
            return
        args = data[:6] + (s,) + data[7:]
        data = dfitpack.fpcurf1(*args)
        if data[-1] == 1:
            data = self._reset_nest(data)
        self._data = data
        self._reset_class()

    def __call__(self, x, nu=0, ext=None):
        x = np.asarray(x)
        if x.size == 0:
            return array([])
        if ext is None:
            ext = self.ext
        else:
            try:
                ext = _extrap_modes[ext]
            except KeyError:
                raise ValueError("Unknown extrapolation mode %s." % ext)
        return fitpack.splev(x, self._eval_args, der=nu, ext=ext)

    def get_knots(self):
        data = self._data
        k, n = data[5], data[7]
        return data[8][k:n-k]

    def get_coeffs(self):
        data = self._data
        k, n = data[5], data[7]
        return data[9][:n-k-1]

    def get_residual(self):
        return self._data[10]

    def integral(self, a, b):
        return dfitpack.splint(*(self._eval_args+(a, b)))

    def derivatives(self, x):
        d, ier = dfitpack.spalde(*(self._eval_args+(x,)))
        if not ier == 0:
            raise ValueError("Error code returned by spalde: %s" % ier)
        return d

    def roots(self):
        k = self._data[5]
        if k == 3:
            z, m, ier = dfitpack.sproot(*self._eval_args[:2])
            if not ier == 0:
                raise ValueError("Error code returned by spalde: %s" % ier)
            return z[:m]
        raise NotImplementedError('finding roots unsupported for '
                                  'non-cubic splines')

    def derivative(self, n=1):
        tck = fitpack.splder(self._eval_args, n)
        # if self.ext is 'const', derivative.ext will be 'zeros'
        ext = 1 if self.ext == 3 else self.ext
        return UnivariateSpline._from_tck(tck, ext=ext)

    def antiderivative(self, n=1):
        tck = fitpack.splantider(self._eval_args, n)
        return UnivariateSpline._from_tck(tck, self.ext)


class InterpolatedUnivariateSpline(UnivariateSpline):
    def __init__(self, x, y, w=None, bbox=[None]*2, k=3,
                 ext=0, check_finite=False):

        x, y, w, bbox, self.ext = self.validate_input(x, y, w, bbox, k, None,
                                            ext, check_finite)
        if not np.all(diff(x) > 0.0):
            raise ValueError('x must be strictly increasing')

        # _data == x,y,w,xb,xe,k,s,n,t,c,fp,fpint,nrdata,ier
        self._data = dfitpack.fpcurf0(x, y, k, w=w, xb=bbox[0],
                                      xe=bbox[1], s=0)
        self._reset_class()