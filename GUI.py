import warnings
import tkinter as tk
import matplotlib.pyplot as pl
from numpy import zeros, concatenate, ravel, diff, array, ones
import numpy as np
from scipy.interpolate import fitpack
from scipy.interpolate import dfitpack
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)

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


# Windows
Window = tk.Tk()

# Judul & Ukuran
Window.title("Spline Interpolation")
Window.geometry('1300x650')

# Label
Array1 = []  # List 1
Array2 = []  # List 2


def Boxes():
    Data2 = tk.Frame(Data1)
    Data2.pack()

    xx = size.get()
    Array1.clear()
    Array2.clear()
    for i in range(xx):
        listA = tk.Entry(Data2)
        listA.grid(row=2+i, column=0, sticky="w")
        Array1.append(listA)
        listB = tk.Entry(Data2)
        listB.grid(row=2+i, column=1, sticky="w")
        Array2.append(listB)

    def Clear():
        Data2.destroy()

    ClearArray = tk.Button(Array, text="Clear", command=Clear)
    ClearArray.grid(row=0, column=3, sticky="w")

# Plot1


def Plot1():
    xx = size.get()
    xnew = float(xs.get())
    x = []
    y = []
    try:
        for i in range(xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning(
            "Warning", f"There's Something Wrong, Element of Array Empty.")

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    spl1 = InterpolatedUnivariateSpline(x, y, k=1)
    Spline1 = tk.Label(Number, text="Spline Orde 1: " + str(spl1(xnew)))
    Spline1.grid(row=2, column=0, sticky="w")

    fig = Figure(figsize=(3, 3),
                 dpi=100)

    pl1 = fig.add_subplot(111)

    x.insert(i, xnew)
    y.insert(i, spl1(xnew))

    pl1.plot(x, y, 'bo-')
    pl1.set_title('Spline Orde 1')

    canvas1 = FigureCanvasTkAgg(fig, master=Plot)
    canvas1.get_tk_widget().grid(row=6, column=0, sticky="w")

# Plot2


def Plot2():
    xx = size.get()
    xnew = float(xs.get())
    x = []
    y = []
    try:
        for i in range(xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning(
            "Warning", f"There's Something Wrong, Element of Array Empty.")

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    spl2 = InterpolatedUnivariateSpline(x, y, k=2)
    Spline2 = tk.Label(Number, text="Spline Orde 2: " + str(spl2(xnew)))
    Spline2.grid(row=3, column=0, sticky="w")

    fig = Figure(figsize=(3, 3),
                 dpi=100)

    pl2 = fig.add_subplot(111)

    x.insert(i, xnew)
    y.insert(i, spl2(xnew))

    pl2.plot(x, y, 'ro-')
    pl2.set_title('Spline Orde 2')

    canvas2 = FigureCanvasTkAgg(fig, master=Plot)
    canvas2.get_tk_widget().grid(row=6, column=1, sticky="w")

# Plot3


def Plot3():
    xx = size.get()
    xnew = float(xs.get())
    x = []
    y = []
    try:
        for i in range(xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning(
            "Warning", f"There's Something Wrong, Element of Array Empty.")

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    spl3 = InterpolatedUnivariateSpline(x, y, k=3)
    Spline3 = tk.Label(Number, text="Spline Orde 3: " + str(spl3(xnew)))
    Spline3.grid(row=4, column=0, sticky="w")

    fig = Figure(figsize=(3, 3),
                 dpi=100)

    pl3 = fig.add_subplot(111)

    x.insert(i, xnew)
    y.insert(i, spl3(xnew))

    pl3.plot(x, y, 'go--')
    pl3.set_title('Spline Orde 3')

    canvas3 = FigureCanvasTkAgg(fig, master=Plot)
    canvas3.get_tk_widget().grid(row=6, column=2, sticky="w")

# Plot1&Plot2&Plot3


def Plot4():
    a = []
    b = []
    x = []
    y = []
    xx = size.get()
    xnew = float(xs.get())
    try:
        for i in range(xx):
            x.append(float(Array1[i].get()))
            y.append(float(Array2[i].get()))
            a.append(float(Array2[i].get()))
            b.append(float(Array2[i].get()))
    except:
        tk.messagebox.showwarning(
            "Warning", f"There's Something Wrong, Element of Array Empty.")

    spl1 = InterpolatedUnivariateSpline(x, y, k=1)
    spl2 = InterpolatedUnivariateSpline(x, y, k=2)
    spl3 = InterpolatedUnivariateSpline(x, y, k=3)

    start = 0
    i = 0
    while x[start] < xnew:
        start += 1
        i += 1

    fig = Figure(figsize=(3, 3),
                 dpi=100)

    pl4 = fig.add_subplot(111)

    x.insert(i, xnew)
    y.insert(i, spl1(xnew))
    a.insert(i, spl2(xnew))
    b.insert(i, spl3(xnew))

    pl4.set_title('Differences')
    pl4.plot(x, y, 'bo-')
    pl4.plot(x, a, 'ro-')
    pl4.plot(x, b, 'go--')

    canvas4 = FigureCanvasTkAgg(fig, master=Plot)
    canvas4.get_tk_widget().grid(row=6, column=3, sticky="w")


# Frame1
Array = tk.Frame(Window)
Array.pack()

text1 = tk.Label(Array, text="Enter the Size of Array:", font="Arial 10 bold")
text1.grid(row=0, column=0, sticky="w")

size = tk.IntVar()

ArraySize = tk.Entry(Array, textvariable=size)
ArraySize.grid(row=0, column=1, sticky="w")

SizeofArray = tk.Button(Array, text="Submit", command=Boxes)
SizeofArray.grid(row=0, column=2, sticky="w")

# Frame2
Data1 = tk.Frame(Window)
Data1.pack()

# Frame3
Data2 = tk.Frame(Data1)
Data2.pack()

# Frame4
Number = tk.Frame(Window)
Number.pack()

text2 = tk.Label(
    Number, text="Enter New Number to Interpolate:", font="Arial 10 bold")
text2.grid(row=1, column=0, sticky="w")

xs = tk.DoubleVar()

NNumber = tk.Entry(Number, textvariable=xs)
NNumber.grid(row=1, column=1, sticky="w")

NNumberB = tk.Button(Number, text="Submit", command=lambda: [
                     Plot1(), Plot2(), Plot3(), Plot4()])
NNumberB.grid(row=1, column=2, sticky="w")

# Frame5
Plot = tk.Frame(Window)
Plot.pack()

Window.mainloop()
