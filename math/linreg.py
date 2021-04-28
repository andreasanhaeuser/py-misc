# standard
from collections import Iterable

# PyPI
import numpy as np
from scipy import odr
from scipy.special import erf
from scipy.special import erfinv

class LinearRegression(object):
    """
        Reference
        ---------
        p. 16 in
        https://www.bpa.gov/EE/Utility/Evaluation/Evaluation/
                UncertaintyMethodsComparisonsFinal.pdf
    """
    def __init__(self, x, y, sigma_level=1):
        self.sigma_level = sigma_level

        # input check
        # =================================================
        x = np.array(x)
        y = np.array(y)

        if len(y) != len(x):
            raise ValueError(
                    'Lengths of x (%i) and y (%i) are not equal.'
                    % (len(x), len(y))
                    )
        # =================================================

        # remove nan's
        # =================================================
        idx = np.isfinite(x) & np.isfinite(y)
        x = x[idx]
        y = y[idx]
        N = len(x)
        # =================================================

        # input check
        # =================================================
        if N < 3:
            raise ValueError('Need at least 3 finite data points, got %i' % N)
        # =================================================

        a, b, Da, Db, N = lin_reg(x, y)
        y_predicted = a * x + b
        variances = (y - y_predicted)**2
        sigma_y = np.sqrt(np.sum(variances)/N)

        self.sigma_y = sigma_y
        self.mean_x = np.mean(x)
        self.sum_x_deviations = np.sum((x - self.mean_x)**2)
        self.slope = a
        self.intercept = b
        self.sigma_slope = Da
        self.sigma_intercept = Db
        self.N = N

    def __call__(self, x):
        return self.model(x)

    ############################################################
    # sigma level & confidence level                           #
    ############################################################
    def set_sigma_level(self, sigma_level):
        """Set sigma level (e. g. 1.96 for 95%-confidence interval."""
        self.sigma_level = sigma_level
        return self

    def set_confidence_level(self, clevel):
        """Set confidence level (e. g. 0.95 for sigma-level 1.96)."""
        if not (0 <= clevel < 1):
            expect = 'Confidence level must be >= 0 and < 1'
            whatis = 'got %s' % str(clevel)
            message = expect + ', ' + whatis
            raise ValueError(message)

        self.sigma_level = np.sqrt(2) * erfinv(clevel)
        return self

    def get_sigma_level(self):
        """Return sigma level as float."""
        return self.sigma_level

    def get_confidence_level(self, clevel):
        """Return sigma level as float between [0..1)."""
        a = self.sigma_level / np.sqrt(2)
        return erf(a)

    ############################################################
    # model & confidence intervals                             #
    ############################################################
    def model(self, x):
        """Return best prediction (linear regression line)."""
        a = self.slope
        b = self.intercept
        return a * x + b

    def model_confidence_interval(self, x):
        """Return lower and upper model confidence bounds."""
        return self._get_interval(x, 'conf')

    def prediction_interval(self, x):
        """Return lower and upper prediction confidence bounds."""
        return self._get_interval(x, 'pred')

    # helpers ----------------------------------------
    def _get_interval(self, x, kind='conf'):
        """Helper function, return confidence interval."""
        # best
        y = self.model(x)

        # spread
        t = self.get_sigma_level()
        sigma = self._get_sigma_value(x, kind)
        Dy = t * sigma

        # lower and upper bounds
        ylo = y - Dy
        yhi = y + Dy

        return ylo, yhi

    def _get_sigma_value(self, x, kind='conf'):
        sigma_y = self.sigma_y
        sqrt = self._get_sqrt(x, kind)

        sigma = sigma_y * sqrt
        return sigma

    def _get_sqrt(self, x, kind='conf'):
        """Helper function, return square root term."""
        if kind.startswith('conf'):
            first_term = 0
        elif kind.startswith('pred'):
            first_term = 1
        else:
            raise ValueError(kind)

        N = self.N
        mean_x = self.mean_x
        numerator = (x - mean_x)**2
        denominator = self.sum_x_deviations

        radicant = first_term + 1/N + numerator / denominator
        return np.sqrt(radicant)

    ############################################################
    # confidence                                               #
    ############################################################
    def get_model_confidence(self, x, ylo, yhi):
        return _get_confidence(self, x, ylo, yhi, 'conf')

    def get_prediction_confidence(self, x, ylo, yhi):
        return _get_confidence(self, x, ylo, yhi, 'pred')

    def _get_confidence(self, x, ylo, yhi, kind):
        if not ylo <= yhi:
            raise ValueError('ylo must be smaller than yhi.')

        mean = self.model(x)
        sigma = self._get_sigma_value(self, x, kind)

        p_lo = probability_than_smaller_than(ylo, mean, sigma)
        p_hi = probability_than_larger_than(yhi, mean, sigma)
        return 1 - (p_lo + p_hi)

def probability_that_larger_than(value, mean, sigma):
    """Return probability that truth is larger than y."""
    a = (value - mean) / (sigma * np.sqrt(2))

    P_inside = erf(np.abs(a))
    P_outside = 1 - P_inside

    if a > 0:
        return P_outside / 2
    return 1 - P_outside / 2

def probability_that_smaller_than(value, mean, sigma):
    """Return probability that truth is larger than y."""
    a = (value - mean) / (sigma * np.sqrt(2))

    P_inside = erf(np.abs(a))
    P_outside = 1 - P_inside

    if a < 0:
        return P_outside / 2
    return 1 - P_outside / 2

def lin_reg(x, y):
    """Linear least square fit of point cloud.

        Parameters
        ----------
        x : array_like
        y : array_like of the same shape as x

        Returns
        -------
        a : slope
        b : y-axis intercept
        Da : uncertainty in a
        Db : uncertainty in b
        N : number of non-NaN data pairs

        Notes
        -----
        The result y = a * x + b is the least square error fit to the points.

        The x-values are considered to be precise. The y-values are con-
        sidered to have all the same uncertainty.
    """
    #### INPUT CHECK ###
    if not isinstance(x, Iterable):
        raise TypeError('x must be iterable.')
    if not isinstance(y, Iterable):
        raise TypeError('y must be iterable.')
    if np.shape(x) != np.shape(y):
        raise IndexError('x and y must be of same length.')

    #### CONVERT DATA ###
    # convert to arrays:
    x = np.array(x[:])
    y = np.array(y[:])

    # exclude nan's:
    idx = np.isfinite(x) & np.isfinite(y)
    x = x[idx]
    y = y[idx]

    #### STATISTICS ###
    N = len(x)  # number of data points
    if N < 3:
        raise ValueError('Need at least 3 data points, got %i' % N)
    f = N - 2   # degrees of freedom

    # slope and intercept:
    Delta = N * np.sum(x*x) - np.sum(x) * np.sum(x)
    a = 1./Delta * (N * np.sum(x*y) - np.sum(x) * np.sum(y))
    b = 1./Delta * (np.sum(x*x) * np.sum(y) - np.sum(x) * np.sum(x*y))

    # y-error:
    dy = (a * x + b) - y                # array
    Dy = np.sqrt(1./f * np.sum(dy*dy))  # ensemble error

    # errors in a and b:
    Da = Dy * np.sqrt(1./Delta * N)
    Db = Dy * np.sqrt(1./Delta * np.sum(x*x))

    return (a, b, Da, Db, N)

def orthogonal_regression(x, y, contains_nans=True):
    """Minimize square distances perpendicular to the regression line.

        A line of the shape
            y = m * x + t                                 (1)
        is returned that reduces the squared distance of the data points to the
        line.

        Conventional linear regression minmizied square distances to the
        regression line in y-direction. This assumes that the x-values are free
        of errors.  However, if x and y are equal in terms of certainty, this
        method yield more "geometric" regression.

        Parameters
        ----------
        x, y : arrays of equal shape
            coordinates of the data
        contains_nans : bool, optional
            set this to True, if you are not sure.

        Returns
        -------
        reg : scipy.odr.odrpack.Output object
            Useful fields are:
            beta : [m, t]
                m, t (floats) are the linear regression parameters in Eq. (1).
            sd_beta : [sd_m, sd_t]
                standard deviations of m and t

        Author
        ------
        Written in 2014
        by Andreas Anhaeuser
        Insitute for Geophysics and Meteorology
        University of Cologne
        Germany
        <anhaeus@meteo.uni-koeln.de>
    """
    # Step 1: filter nan's
    # ====================
    #
    # Step 2: find a first guess
    # ==========================
    # the result of the conventional linear regression is used as first guess
    # (beta0).
    #
    # Step 3: perform regression
    # ==========================
    assert np.shape(x) == np.shape(y)

    # Step 1: filter nan's
    if contains_nans:
        idx = np.isnan(x) | np.isnan(y)  # abort these entries
        x = x[-idx]
        y = y[-idx]

    # Step 2: first guess
    # vertical least square regression:
    a, b, Da, Db, N = lin_reg(x, y)

    # Step 3: regression
    f = lambda x: beta[0] * x + beta[1]    # fitting function
    beta0 = [a, b]                          # first guess
    model = odr.Model(f)
    data = odr.Data(x, y)
    myodr = odr.ODR(mydata, model, beta0)
    return myodr.run()
