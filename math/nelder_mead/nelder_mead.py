#!/usr/bin/python
"""Minimization of scalar functions using Nelder-Mead algorithm.

    This solver is intended to be an improvement of the Nelder-Mead function
    in `scipy.optimize.minizme`.

    Todo
    ====
    add max_time


    Improvements are:
    =================
    - The output now also includes the history of the function values on the
      vertices of each iteration step.
    - The codes has been completely re-written in order to make it more
      readable.
    - Documentation has been added.


    Terms
    =====
        Simplex
        -------
        A simplex is the higher-dimensional generalisation of a triangle on a
        plane or a tetrahedron in space. It exists in N-dimensional
        (hyper-)space and has (N+1) "corners".

        Vertex
        ------
        A "corner" of a simplex.

        Centroid
        --------
        The arithmetic mean of a simplex's vertices (i. e. its "center" point).


    Method
    ======
    The algorithm tries to minimize a scale-valued function in N-dimensional
    argument space. An N-dim simplex is initialized an the iteratively moved
    and reshaped depending on the function values on its vertices. If
    successful, this will eventually lead to a contraction of the simplex
    around the function minimum. For details see docstring of the main function.


    Author
    ======
    Andreas Anhaeuser (AA)
    University of Cologne
    <andreas.anhaeuser@posteo.net>
"""

# standard modules
from copy import deepcopy as copy
import datetime as dt

# PyPI modules
import numpy as np
from scipy.optimize import OptimizeResult

# local modules
from . import messages

_alpha = 1.     # reflection factor
_gamma = 2.     # expansion factor
_rho = 0.5      # contraction factor
_sigma = 0.5    # shrinking factor

_fatol = 1e-4   # function difference tolerance
_xatol = 1e-4   # parameter difference tolerance
_maxiter = 200  # maximum number of iterations
_maxfev = 400   # maximum number of function evaluation
_maxtime_sec = 5 * 60   # maximum real time to use

# simplex : ndarray, length N
# vertex : ndarray, shape (N+1) x N

###################################################
# MAIN                                            #
###################################################
def minimize(fun, x0, args=(), method='nelder-mead', jac=None, hess=None,
        hessp=None, bounds=None, constraints=(), tol=None, callback=None,
        options=None):
    if method.lower() != 'nelder-mead':
        raise NotImplementedError()
    
    return minimize_nelder_mead(fun, x0, args=args, callback=callback,
            **options)

def minimize_nelder_mead(
        fun,
        x0,
        args=(),
        callback=None,
        maxiter=None,
        maxfev=None,
        maxtime_sec=None,
        disp=False,
        return_all=False,
        initial_simplex=None,
        xatol=None,
        fatol=None,
        adaptive=False,
        x0_uncert=None,
        alpha=_alpha,
        gamma=_gamma,
        rho=_rho,
        sigma=_sigma,
        parameter_names=None,
    ):
    """Return a dict.

        Parameters
        ----------
        fun : callable
            Function that is to be minimized.
        x0 : array of shape (N,)
            Initial guess of the parameters to be optimized.
        args : tuble, optional
            Additional parameters (other than `x0`) passed to the objective
            function `fun`.
        callback : callable, optional
            Called after each iteration. If callback returns True the algorithm
            execution is terminated. The signature is:
            
                ``callback(xk)``
            
            where ``xk`` is the current parameter vector.
        maxiter : int, optional
            Maximum number of iterations.
        maxfev : int, optional
            Maximum number of function evaluations.
        maxtime_sec : int, optional
            Maximum number of seconds to use.
        disp : bool, optional
            Print convergence messages. Default: False
        return_all : bool, optional
            Return complete iteration record. Default: False
        initial_simplex : array of shape (N + 1, N), optional
            If given, overrides `x0` and `x0_uncert`. ``initial_simplex[j,:]``
            contains the coordinates of the j-th vertex of the ``N+1`` vertices
            in the simplex, where ``N`` is the dimension.
        xatol : array of length N, optional
            Maximum difference of the vertex coordinates upon termination
        fatol : float, optional
            Maximum difference of function value to the previous one upon
            termination
        adaptive : bool, optional
            Unused remnant of the original numpy implementation.
            Caution: If True, function will exit with ``NotImplementedError``.
            Default: False.
        x0_uncert : array of shape (N,), optional
            Initial uncertaintainties of the parameters in `x0`. If
            `initial_simplex` is given, `x0_uncert` remains unused. Otherwise,
            it is used to construct `initial_simplex`. It neither of the two is
            given, the function will still run (by just making up some
            `x0_uncert` itself), but this is discouraged as the function has no
            means to estimate `x0_uncert` is a sensible way.

        Returns
        -------
        scipy.optimize.OptimizeResult

        History
        -------
        2018-04-28 (AA): Created
        2018-05-14 (AA): Conformed to scipy.optimize._minimize_neldermead
    """
    ###################################################
    # DEFAULTS                                        #
    ###################################################
    if fatol is None:
        fatol = _fatol
    if xatol is None:
        xatol = _xatol * np.ones_like(x0)
    if maxiter is None:
        maxiter = _maxiter
    if maxfev is None:
        maxfev = _maxfev
    if maxtime_sec is None:
        maxtime_sec = _maxtime_sec

    if x0_uncert is None:
        x0_uncert = make_up_uncertainties(x0)

    if initial_simplex is None:
        initial_simplex = get_initial_simplex(x0, x0_uncert)

    if adaptive:
        raise NotImplementedError('Adaptive parameters not implemented. ' + 
                'Set `adaptive` to False.')

    ###################################################
    # CHECKS                                          #
    ###################################################
    if any(x0_uncert==0):
        raise ValueError('Uncertaities must not be 0.')
    
    ###################################################
    # INITIALIZE                                      #
    ###################################################
    options = {
            'Nparam' : len(x0),
            'maxiter' : maxiter,
            'maxfev' : maxfev,
            'maxtime' : dt.timedelta(seconds=maxtime_sec),
            'fatol' : fatol,
            'xatol' : xatol,
            'alpha' : alpha,
            'gamma' : gamma,
            'rho' : rho,
            'sigma' : sigma,
            }

    M, N = np.shape(initial_simplex)
    table = initialize_lookup_table(initial_simplex[0])
    record = OptimizeResult(
            {
            'time_beg' : dt.datetime.now(),
            'all_simplices' : [],
            'best_vertices' : [],
            'worst_vertices' : [],
            'all_fun_values': [],
            'best_fun_values' : [],
            'worst_fun_values' : [],
            'lookup_table' : table,
            'action' : [],              # mode in which simplex is manipulated
            }
            )

    ###################################################
    # ABBREVIATIONS                                   #
    ###################################################
    f = lambda vertex: get_function_value(vertex, fun, args, table)
    F = lambda simplex: get_function_values(simplex, fun, args, table)

    # action flags
    _INITIALIZE = 0
    _SHRINK = 1
    _CONTRACT = 2
    _REFLECT = 3
    _EXPAND = 4

    ###################################################
    # ITERATE                                         #
    ################################################### 
    simplex = initial_simplex
    record['action'].append(_INITIALIZE)

    while True:
        ###################################################
        # EVALUATE                                        #
        ###################################################
        fun_values = F(simplex)

        # best
        m_best = np.argmin(fun_values)
        x_best = simplex[m_best]
        f_best = fun_values[m_best]                    

        # worst
        m_worst = np.argmax(fun_values)
        x_worst = simplex[m_worst]
        f_worst = fun_values[m_worst]

        ###################################################
        # RECORD                                          #
        ###################################################
        record['all_simplices'].append(copy(simplex))
        record['best_vertices'].append(copy(x_best))
        record['worst_vertices'].append(copy(x_worst))
        record['all_fun_values'].append(copy(fun_values))
        record['best_fun_values'].append(copy(f_best))
        record['worst_fun_values'].append(copy(f_worst))
        record['x0'] = x0
        if x0_uncert is not None:
            record['x0_unc'] = x0_uncert
        else:
            record['x0_unc'] = np.ptp(initial_simplex, 0)
        record['f0'] = get_function_value(x0, fun, args, table)
        if 'f0_unc' not in record.keys():
            record['f0_unc'] = f_worst - f_best

        if parameter_names is None:
            parameter_names = ['parameter %i' % n for n in range(N)]
        record['parameter_names'] = parameter_names
        
        ###################################################
        # DISPLAY                                         #
        ###################################################
        if disp:
            update_record(record, options)
            messages.display_intermediate_message(record, options)

        ###################################################
        # TERMINATE                                       #
        ###################################################
        if termination_reached(record, options):
            break

        ###################################################
        # FIND A BETTER SIMPLEX                           #
        ###################################################
        # centroid
        x_centroid = get_centroid_vertex(simplex, m_worst)

        # reflected
        x_reflected = get_reflected_vertex(x_worst, x_centroid,
                options['alpha'])
        f_reflected = f(x_reflected)

        # check quality of reflected vertex
        if f_reflected < f_best:
            # x_reflected is best
            # -> try expansion
            x_expanded = get_expanded_vertex(x_reflected, x_centroid,
                    options['gamma'])
            f_expanded = f(x_expanded)
            if f_expanded < f_reflected:
                # expanded is better
                simplex[m_worst] = x_expanded
                record['action'].append(_EXPAND)
            else:
                # reflected is better
                simplex[m_worst] = x_reflected
                record['action'].append(_REFLECT)

        elif f_reflected > f_worst:
            # x_reflected is worst
            # -> try contraction
            x_contracted = get_contracted_vertex(x_worst, x_centroid,
                    options['rho'])
            f_contracted = f(x_contracted)
            if f_contracted < f_worst:
                # contracted is better
                # -> use it
                simplex[m_worst] = x_contracted
                record['action'].append(_CONTRACT)
            else:
                # contracted is even worse (this can happen close to minima)
                # -> shrink
                simplex = get_shrunk_simplex(simplex, x_centroid, m_best,
                        options['sigma'])
                record['action'].append(_SHRINK)

        else:
            # x_reflected is neither worst nor best
            # -> use it
            simplex[m_worst] = x_reflected
            record['action'].append(_REFLECT)

        if callback is not None:
            terminate = callback(x_best)
            if terminate is True:
                break

    ###################################################
    # CREATE OptimizeResult OBJECT                    #
    ###################################################
    record = update_record(record, options)

    # list -> array
    for key in record.keys():
        if key in ['lookup_table']:
            continue
        record[key] = np.array(record[key])
    record['success'] = convergence_reached(record, options)
    record['message'] = messages.get_termination_message(record, options)

    return record

###################################################
# MANIPULATE SIMPLEX                              #
###################################################
def get_centroid_vertex(simplex, exclude_index):
    """Return a vertex.

        Parameters
        ----------
        simplex : array of shape (N + 1, N)
        exclude_index : int (0...N)
            index of the vertex that is not taken into account

        Returns
        -------
        centroid : array of shape (N,)
    """
    # ========== input check  ============================ #
    shape = np.shape(simplex)
    assert len(shape) == 2
    M, N = shape
    assert M == N + 1
    assert exclude_index < M
    # ==================================================== #

    # exclude one vertex
    m = exclude_index
    vertices = np.concatenate((simplex[:m], simplex[m+1:]), 0)

    # compute
    centroid = np.mean(vertices, 0)
    return centroid

def get_reflected_vertex(vertex, centroid, alpha):
    """Return a vertex.

        Parameters
        ----------
        vertex : array of shape (N,)
        centroid : array of shape (N,)
        alpha: float
            reflection factor (1 for a 'true' reflection)

        Returns
        -------
        reflected_vertex : array of shape (N,)
    """
    assert alpha > 0
    return centroid + alpha * (centroid - vertex)

def get_expanded_vertex(vertex, centroid, gamma):
    """Return a vertex.

        Parameters
        ----------
        vertex : array of shape (N,)
        centroid : array of shape (N,)
        gamma: float
            expansion factor

        Returns
        -------
        expanded_vertex : array of shape (N,)
    """
    assert gamma > 1
    return centroid + gamma * (vertex - centroid)

def get_contracted_vertex(vertex, centroid, rho):
    """Return a vertex.

        Parameters
        ----------
        vertex : array of shape (N,)
        centroid : array of shape (N,)
        rho: float
            contraction factor

        Returns
        -------
        contracted_vertex : array of shape (N,)
    """
    assert 0 < rho < 1
    return centroid + rho * (vertex - centroid)

def get_shrunk_simplex(simplex, centroid, best_index, sigma):
    """Return a simplex.

        Parameters
        ----------
        simplex : array of shape (N + 1, N)
        centroid : array of shape (N,)
        best_index : int (0...N)
            index of the vertex that is to be kept unchanged
        sigma: float
            shrink factor

        Returns
        -------
        shrunk_simplex : array of shape (N + 1, N)
    """
    # ========== input check  ============================ #
    shape = np.shape(simplex)
    assert len(shape) == 2
    M, N = shape
    assert M == N + 1
    assert 0 <= best_index < M

    assert 0 < sigma < 1

    # ========== shrink  ================================= #
    best_vertex = simplex[best_index]
    new_simplex = np.nan * np.ones_like(simplex)

    for m in range(M):
        vertex = simplex[m]

        if m != best_index:
            # non-best vertex: shrink towards best_vertex
            new_vertex = best_vertex + sigma * (vertex - best_vertex)
        else:
            # best vertex
            new_vertex = copy(vertex)

        new_simplex[m] = new_vertex

    # ========== output check  =========================== #
    assert np.sum(np.isnan(new_simplex)) == 0

    return new_simplex

###################################################
# FORWARD FUNCTION                                #
###################################################
def get_function_values(simplex, fun, args, lookup_table):
    """Return an array."""
    # input check
    shape = np.shape(simplex)
    assert len(shape) == 2
    M, N = shape
    assert M == N + 1

    # initialize
    fvals = np.nan * np.ones(M)

    # fill array
    for m in range(M):
        vertex = simplex[m]
        fvals[m] = get_function_value(vertex, fun, args, lookup_table)

    return fvals

def get_function_value(vertex, fun, args, lookup_table):
    """Return a float."""
    # input check
    shape = np.shape(vertex)
    assert len(shape) == 1

    # check loopup_table
    result = get_function_value_from_lookup_table(lookup_table, vertex)

    # call function
    if result is None:
        result = fun(vertex, *args)
        add_function_value_to_loopkup_table(lookup_table, vertex, result)

    assert np.shape(result) == ()
    return result

###################################################
# LOOKUP TABLE                                    #
###################################################
def initialize_lookup_table(vertex):
    """Return a dict of function values on vertices."""
    N = len(vertex)
    return {
            'vertices' : np.ones((0, N)),
            'fun_values' : np.ones(0),
            }

def get_function_value_from_lookup_table(lookup_table, vertex):
    """Return result if vertex is in lookup_table, else None."""
    index = find_index_in_lookup_table(lookup_table, vertex)
    if index is None:
        result = None
    else:
        result = lookup_table['fun_values'][index]

    return result

def find_index_in_lookup_table(lookup_table, vertex):
    """Return index if vertex is in lookup_table, else None."""
    vertices = lookup_table['vertices']
    K = len(vertices)
    index = None

    # try to find vertex in list of vertices
    for k in range(K):
        is_equal = np.sum(vertices[k] != vertex) == 0
        if is_equal:
            index = k
            break

    return index

def add_function_value_to_loopkup_table(lookup_table, vertex, result):
    index = find_index_in_lookup_table(lookup_table, vertex)
    if index is not None:
        return 

    table = lookup_table
    K = len(table['vertices'])
    table['vertices'] = np.insert(table['vertices'], K, vertex, 0)
    table['fun_values'] = np.insert(table['fun_values'], K, result, 0)

###################################################
# TERMINATION                                     #
###################################################
def termination_reached(record, options):
    """Return a bool."""
    iteration_criterion = iteration_limit_reached(record, options)[0]
    convergence_criterion = convergence_reached(record, options)
    terminate = iteration_criterion or convergence_criterion
    return terminate

def iteration_limit_reached(record, options):
    """Return (bool, str)."""
    update_record(record, options)

    it_lim_reached = False
    message = ''

    # maxiter
    if record['Niter'] > options['maxiter']:
        it_lim_reached = True
        message = 'Maximum number of iterations reached.'

    # maxfev
    if record['Nfev'] > options['maxfev']:
        it_lim_reached = True
        message = 'Maximum number of function evaluations reached.'

    # too little iterations
    if record['Niter'] < 2:
        it_lim_reached = False
        message = ''

    # maxtime
    if record['time_passed'] > options['maxtime']:
        it_lim_reached = True
        message = 'Time limit surpassed.'

    record['it_lim_reached'] = it_lim_reached, message

    return it_lim_reached, message

def convergence_reached(record, options):
    """Return a bool."""
    update_record(record, options)

    conv_reached = True

    if record['Niter'] < 2:
        conv_reached = False

    # fatol
    span = function_uncertainty(record, options)
    fatol = options['fatol']
    if not span < fatol:
        conv_reached = False

    # xatol
    spans = parameter_uncertainty(record, options)
    xatol = options['xatol']
    if np.sum(spans >= xatol) > 0:
        conv_reached = False

    record['convergence_reached'] = conv_reached

    return conv_reached

def function_uncertainty(record, options):
    """Return a float.

        Return maximum difference in function values of all points of the last
        two simplices. Return infinity if history is shorter than 2.
    """
    last_fun_values = record['all_fun_values'][-2:]

    # history too short:
    if len(last_fun_values) < 2:
        return np.inf

    return np.ptp(last_fun_values)

def parameter_uncertainty(record, options):
    """Return an array of shape (N,).

        Return maximum difference in parameter values of all points of the last
        two simplices. Return infinity if history is shorter than 2.
    """
    all_simplices = record['all_simplices']
    K, M, N = np.shape(all_simplices)

    # history too short:
    if K < 2:
        return np.ones(N) * np.inf

    # retrieve and reshape last two simplices
    last_simplices = np.array(all_simplices[-2:])
    last_simplices_reshaped = last_simplices.reshape((M*2, N))

    # return maximum differences
    return np.ptp(last_simplices_reshaped, 0)

###################################################
# MISC                                            #
###################################################
def update_record(record, options):
    table = record['lookup_table']
    record['Nfev'] = len(table['fun_values'])
    record['Niter'] = len(record['all_simplices'])
    record['time_passed'] = dt.datetime.now() - record['time_beg']
    record['x_unc'] = parameter_uncertainty(record, options)
    record['f_unc'] = function_uncertainty(record, options)
    record['f'] = record['best_fun_values'][-1]
    record['x'] = record['best_vertices'][-1]
    return record

def get_initial_simplex(x0, x0_uncert):
    """Return a (N+1 x N)-array."""
    # input check
    shape = np.shape(x0)
    assert len(shape) == 1
    assert np.shape(x0_uncert) == shape

    # initialize
    N = shape[0]
    simplex = np.nan * np.zeros((N+1, N))

    # ========== create vertices  ======================== #
    for n in range(-1, N):
        vertex = 1. * copy(x0)

        # manipulate n-th coordinate
        if n >= 0:
            vertex[n] += 0.5 * x0_uncert[n]

        # add vertex to simplex
        simplex[n] = vertex
    # ==================================================== #

    # output check
    assert np.sum(np.isnan(simplex)) == 0

    return simplex

def make_up_uncertainties(x0):
    unc = 0.1 * np.abs(x0) 
    unc[unc==0] = 0.1
    return unc

###################################################
# TEST FUNCTIONS                                  #
###################################################
def himmelblau(a):
    x = a[0]
    y = a[1]
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def styblinski(a):
    N = len(a)
    return np.sum(a**4 - 16 * a**2 + 5 * a) + N * 2 * 39.16617

def styblinski_slow(a):
    N = len(a)
    sleep(0.03)
    return np.sum(a**4 - 16 * a**2 + 5 * a) + N * 2 * 39.16617


###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    from time import sleep
    reload(messages)
    # x0 = - np.array([2.2, 1.5, 2.0, 16., 1., 4., 5.])
    # x0_uncert = np.array([5., 5., 5., 10., 1.,  1., 1.])
    maxtime_sec = 5
    N = 9
    x0 = - np.array([2.2]*N)
    x0_uncert = np.array([5.]*N)
    # f = himmelblau
    f = styblinski
    # f = styblinski_slow

    record = minimize_nelder_mead(f, x0, x0_uncert=x0_uncert, disp=True,
            maxfev=np.inf, maxiter=np.inf, maxtime_sec=maxtime_sec,
            fatol=10e-10)
