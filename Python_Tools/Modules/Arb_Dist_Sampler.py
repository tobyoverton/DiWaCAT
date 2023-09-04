import numpy as np
import scipy.interpolate as interpolate

try:
    from pygam import LinearGAM,s
except ImportError:
    print("not imported pygam - only needed for specific testing")

class Distribution(object):
    """
    draws samples from a one dimensional probability distribution,
    by means of inversion of a discrete inverstion of a cumulative density function

    the pdf can be sorted first to prevent numerical error in the cumulative sum
    this is set as default; for big density functions with high contrast,
    it is absolutely necessary, and for small density functions,
    the overhead is minimal

    a call to this distibution object returns indices into density array

    from:
    https://stackoverflow.com/questions/21100716/fast-arbitrary-distribution-random-sampling-inverse-transform-sampling

    """
    def __init__(self, pdf, sort = True, interpolation = True, transform = lambda x: x):
        self.shape          = pdf.shape
        self.pdf            = pdf.ravel()
        self.sort           = sort
        self.interpolation  = interpolation
        self.transform      = transform

        #a pdf can not be negative
        assert(np.all(pdf>=0))

        #sort the pdf by magnitude
        if self.sort:
            self.sortindex = np.argsort(self.pdf, axis=None)
            self.pdf = self.pdf[self.sortindex]
        #construct the cumulative distribution function
        self.cdf = np.cumsum(self.pdf)
    @property
    def ndim(self):
        return len(self.shape)
    @property
    def sum(self):
        """cached sum of all pdf values; the pdf need not sum to one, and is imlpicitly normalized"""
        return self.cdf[-1]
    def __call__(self, N):
        """draw """
        #pick numbers which are uniformly random over the cumulative distribution function
        choice = np.random.uniform(high = self.sum, size = N)
        #find the indices corresponding to this point on the CDF
        index = np.searchsorted(self.cdf, choice)
        #if necessary, map the indices back to their original ordering
        if self.sort:
            index = self.sortindex[index]
        #map back to multi-dimensional indexing
        index = np.unravel_index(index, self.shape)
        index = np.vstack(index)
        #is this a discrete or piecewise continuous distribution?


        if self.interpolation:
            index = index + np.random.uniform(size=index.shape)
        return self.transform(index)

class Distribution_V2(object):
    """
    draws samples from a one dimensional probability distribution,
    by means of inversion of a discrete inverstion of a cumulative density function

    the pdf can be sorted first to prevent numerical error in the cumulative sum
    this is set as default; for big density functions with high contrast,
    it is absolutely necessary, and for small density functions,
    the overhead is minimal

    a call to this distibution object returns indices into density array

    from:
    https://stackoverflow.com/questions/21100716/fast-arbitrary-distribution-random-sampling-inverse-transform-sampling

    """
    def __init__(self, pdf, sort = True, interpolation = True, transform = lambda x: x):
        self.shape          = pdf.shape
        self.pdf            = pdf.ravel()
        self.sort           = sort
        self.interpolation  = interpolation
        self.transform      = transform

        #a pdf can not be negative
        assert(np.all(pdf>=0))

        #sort the pdf by magnitude
        if self.sort:
            self.sortindex = np.argsort(self.pdf, axis=None)
            self.pdf = self.pdf[self.sortindex]
        #construct the cumulative distribution function
        self.cdf = np.cumsum(self.pdf)
    @property
    def ndim(self):
        return len(self.shape)
    @property
    def sum(self):
        """cached sum of all pdf values; the pdf need not sum to one, and is imlpicitly normalized"""
        return self.cdf[-1]
    def __call__(self, N):
        """draw """
        #pick numbers which are uniformly random over the cumulative distribution function
        choice = np.random.uniform(high = self.sum, size = N)
        #find the indices corresponding to this point on the CDF
        index = np.searchsorted(self.cdf, choice)
        #if necessary, map the indices back to their original ordering
        if self.sort:
            index = self.sortindex[index]
        #map back to multi-dimensional indexing
        index = np.unravel_index(index, self.shape)
        index = np.vstack(index)
        #is this a discrete or piecewise continuous distribution?

        #TODO implement Gaussian blurring here?


        if self.interpolation:
            index = index + np.random.uniform(size=index.shape)
        return self.transform(index)

class Distribution_Interp1D(object):

    def __init__(self,pdf,x,kind = "linear",method=None):

        #TODO add some dimension checking - 1D only

        self.method = method

        #PDF must be a numpy array
        assert isinstance(pdf,np.ndarray)
        assert isinstance(x, np.ndarray)

        #a pdf can not be negative
        assert(np.all(pdf>=0))

        self.kind = kind

        self.summed_up = np.sum(pdf)

        self.pdf = pdf / self.summed_up

        self.shape = pdf.shape
        #self.pdf = np.divide(pdf,np.max(pdf)) # no ravel here, we aren't horrible
        self.x_width = x[1] - x[0]

        #This fixes the offset as it sets the CDF x values to the centre of the PDF bins
        self.cdf_x = x + self.x_width / 2

        self.cdf = np.zeros(self.shape)

        self.pdf_x = x

        #TODO fix offset here by half apparent bin width
        self.cdf[1:] = np.cumsum(self.pdf[1:])



    @property
    def ndim(self):
        return len(self.shape)

    @property
    def sum(self):
        """cached sum of all pdf values; the pdf need not sum to one, and is imlpicitly normalized"""
        return self.cdf[-1]

    @property
    def inv_cdf(self):
        if self.method == "GAM":
            gam = LinearGAM(s(0,constraints="monotonic_inc"),max_iter=200).fit(self.cdf, self.cdf_x)
            return gam.predict

        elif self.method == "PCHIP":
            return interpolate.PchipInterpolator(self.cdf, self.cdf_x)

        elif self.method == "Spline":
            return interpolate.UnivariateSpline(self.cdf, self.cdf_x, k=3)

        elif self.method == "HermiteSpline":
            return interpolate.CubicHermiteSpline(self.cdf, self.cdf_x, np.ones(self.cdf_x.shape))

        elif self.method == "Akima":
            return interpolate.Akima1DInterpolator(self.cdf, self.cdf_x)

        else:
            return interpolate.interp1d(self.cdf, self.cdf_x, self.kind)

    @property
    def min_x(self):
        return np.min(self.cdf_x)

    @property
    def max_x(self):
        return np.max(self.cdf_x)


    def __call__(self, N):
        """draw """
        # pick numbers which are uniformly random over the cumulative distribution function
        choice = np.random.default_rng().uniform(low=0, high=self.sum, size=N)

        choice.sort()

        return self.inv_cdf(choice)

class Distribution_Interp2D(object):

    def __init__(self,pdf,x,y):

        #TODO add some dimesnion checking - 2D

        assert pdf.shape == (len(x),len(y))

        #PDF must be a numpy array
        assert isinstance(pdf,np.ndarray)
        assert isinstance(x, np.ndarray)

        #a pdf can not be negative
        assert(np.all(pdf>=0))

        self.shape = pdf.shape

        self.pdf = pdf

        self.summed_up = np.sum(self.pdf)

        self.pdf = self.pdf / self.summed_up

        #This is for plotting of the generating PDF
        self.pdf_x = x
        self.pdf_y = y

        self.x_width = x[1]-x[0]
        self.y_width = y[1]-y[0]

        self.cdf_x = x + self.x_width/2
        self.cdf_y = y + self.y_width/2

        #Need to goddam build a CDF

        self.cdf = self.buildCDF()

        #And a marginal CDF in X
        self.marg_x_cdf = self.build_marginal_CDF()

        #And this are popualted on call
        self.yi = None
        self.y_data_for_sample = None
        self.y_sample = None
        self.x_sample = None

        self.progress_counter = None
        self.percent_saver = None
        self.ten_percent_point = None

    def buildCDF(self):
        #Could make this Ndimension quite easily
        cdf = np.cumsum(self.pdf,axis=0) # sum across here
        cdf = np.cumsum(cdf,axis=1) # then sum across here
        cdf = cdf / cdf.max()



        return cdf

    def build_marginal_CDF(self):
        cdf = np.sum(self.cdf,axis=0)

        return cdf

    def set_x_values(self,N):
        x_choices = np.random.default_rng().uniform(low=min(self.marg_x_cdf),
                                                   high=max(self.marg_x_cdf),
                                                   size=N)
        x_choices.sort()


        self.x_sample = self.inv_marg_CDF(x_choices)

    def set_y_data_for_sample(self):
        self.y_data_for_sample = interp_pchip_n(self.cdf_x, self.cdf_y,
                                                self.cdf,
                                                self.x_sample, self.yi)

    def progress_reporter(self):
        #Overhead here?
        if self.progress_counter == self.ten_percent_point:
            self.percent_saver += 10
            print('Progress report:', self.percent_saver, "%")
            self.progress_counter = 0

    def set_y_values(self):
        temp_ydata = []
        self.progress_counter = 0

        self.ten_percent_point = int(len(self.x_sample)/10)
        self.percent_saver=0

        self.progress_counter = 0
        for i, set_x in enumerate(self.x_sample):

            self.progress_reporter()
            self.progress_counter += 1

            set_cdf = self.y_data_for_sample[i]
            inv_set_CDF = interpolate.PchipInterpolator(set_cdf, self.yi)
            y_sample = np.random.default_rng().uniform(low=min(set_cdf), high=max(set_cdf))

            y_val = inv_set_CDF(y_sample)

            temp_ydata.append(y_val)

        print("Progress report: 100 %")
        self.y_sample = np.asarray(temp_ydata)

    @property
    def inv_marg_CDF(self):
        print("debugging here with x len {}, y len {}".format(len(self.marg_x_cdf),len(self.cdf_x)))
        is_sorted = lambda a: np.all(a[:-1] <= a[1:])
        print("The marginal CDF is sorted: ",is_sorted(self.marg_x_cdf))

        return interpolate.PchipInterpolator(self.marg_x_cdf, self.cdf_x)


    def __call__(self, N,y_interp_points = 300):
        """draw """
        # pick numbers which are uniformly random over the cumulative distribution function
        assert type(N) == int
        assert type(y_interp_points) == int

        print("Sampling",N," x values")
        self.set_x_values(N)


        #y interpoaltion values
        print("Using", y_interp_points,"Y values for interpolation")
        self.yi = np.linspace(self.cdf_y.min(), self.cdf_y.max(), y_interp_points)

        print("Retrieving interpolated y sample")
        self.set_y_data_for_sample()

        print("Evaluating y sample for each x point")
        self.set_y_values()

        return np.array((self.x_sample,self.y_sample))


#My special N Pchip function
def interp_pchip_n(*args, **kw):
    """Interpolation on N-D.

    ai = interpn(x, y, z, ..., a, xi, yi, zi, ...)
    where the arrays x, y, z, ... define a rectangular grid
    and a.shape == (len(x), len(y), len(z), ...)
    Modified by t pacey from:
    https://github.com/scipy/scipy/issues/2246
    """
    if kw:
        raise ValueError("Unknown arguments: " % kw.keys())
    nd = (len(args)-1)//2
    if len(args) != 2*nd+1:
        raise ValueError("Wrong number of arguments")
    q = args[:nd]
    qi = args[nd+1:]
    a = args[nd]
    # TODO method control here similar to the other one?
    for j in range(nd):
        a = interpolate.PchipInterpolator(q[j], a, axis=j)(qi[j])
    return a


if __name__=='__main__':
    shape = 3,3
    pdf = np.ones(shape)
    pdf[1]=0
    dist = Distribution(pdf, transform=lambda i:i-1.5)
    print(dist(10))
    import matplotlib.pyplot as pp
    pp.scatter(*dist(1000))
    pp.show()