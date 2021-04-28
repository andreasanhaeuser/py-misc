#!/usr/bin/python
"""Plot result of nelder_mead.minimize_nelder_mead

    Author
    ======
    Andreas Anhaeuser (AA)
    University of Cologne
    <andreas.anhaeuser@posteo.net>
"""

# PyPI modules
import numpy as np
from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt

_actionlabels = ('init', 'shr.', 'contr.', 'refl.', 'exp.')
_actionmarkers = ('d', '.', 'o', '^', 's')

def plot_all_parameters(record, filename, fun_log=False, fun_min=None,
        fun_max=None, figsize=None, parameter_names=None):
    ###################################################
    # SETUP                                           #
    ###################################################
    # figure size
    if figsize is None:
        figsize = (160/25.4, 240/25.4)

    # colors
    col_f = (0.5,) * 3
    col_param = 'r'
    col_param_thin = (1, 0.5, 0.5)
    col_action = (0.3, 0.3, 1.)

    # linestyles
    ls_grid = ':'

    # linewidths
    lw_thick = 1
    lw_thin = 0.5
    lw_grid = 0.25

    # markersizes
    ms_action = 1
    ms_param = 1

    # fontsizes
    fs_fun = 8

    ###################################################
    # RETRIEVE VALUES                                 #
    ###################################################
    best_params = record['best_vertices']
    best_fvals = record['best_fun_values']
    worst_params = record['worst_vertices']
    all_simplices = record['all_simplices']
    worst_fvals = record['worst_fun_values']
    actions = record['action']

    ###################################################
    # FUNCTION                                        #
    ###################################################
    # min/max
    if fun_min is None:
        fun_min = np.nanmin(best_fvals)
    if fun_max is None:
        fun_max = np.nanmax(worst_fvals)
        
    # plot values and parameters
    if fun_log:
        # logarithmic
        plot_best = np.log10(best_fvals)
        plot_worst = np.log10(worst_fvals)
        ylabel_fun = 'log10(function value)'
        ymin_fun = np.log10(fun_min)
        ymax_fun = np.log10(fun_max)
    else:
        # linear
        plot_best = best_fvals
        plot_worst = worst_fvals
        ylabel_fun = 'function value'
        ymin_fun = fun_min
        ymax_fun = fun_max

    ###################################################
    # CREATE PLOT                                     #
    ###################################################
    K, M, N = np.shape(all_simplices)
    Nplots = N + 1
    Ncols = int(np.ceil(np.sqrt(Nplots)))
    if Ncols > 2:
        Ncols -= 1
    Nrows = int(np.ceil(1. * Nplots / Ncols))

    plt.figure(figsize=figsize)
    for nplot in range(1, Nplots+1):
        ax1 = plt.subplot(Nrows, Ncols, nplot)
        ncol = (nplot - 1) % Ncols
        nrow = (nplot - 1) // Ncols

        ###################################################
        # SIMPLEX TRANSFORMATION                          #
        ###################################################
        if nplot == Nplots:
            # color
            color1 = col_action

            # line
            ax1.plot(actions, '.-', lw=lw_thick, color=color1, ms=ms_action)

            # labels
            ax1.set_yticks(range(len(_actionlabels)))
            ax1.set_yticklabels(_actionlabels)

            # title
            title = 'simplex transformation'

        ###################################################
        # PARAMETER VALUE                                 #
        ###################################################
        else:
            nparam = nplot - 1

            # color
            color1 = col_param


            # all
            for m in range(M):
                param = all_simplices[:, m, nparam]
                ax1.plot(param, '-', color=col_param_thin,
                    lw=lw_thin, zorder=10)

            # best
            ax1.plot(best_params[:, nparam], '.-', color=color1,
                    lw=lw_thick, zorder=11, ms=ms_param)

            # title
            if parameter_names is None:
                title = 'parameter %i' % nparam 
            else:
                title = parameter_names[nparam]

        # title
        plt.title(title)

        # grid
        ax1.tick_params('y', colors=color1)
        ax1.grid(color=color1, lw=lw_grid, ls=ls_grid)

        # x-axis
        ax1.set_xlim(0, K-1)

        # zorder
        ax1.set_zorder(10)
        ax1.patch.set_visible(False)

        ###################################################
        # FUNCTION                                        #
        ###################################################
        ax2 = ax1.twinx()
        ax2.fill_between(np.arange(K), plot_best, ymin_fun, color=col_f)
        ax2.plot(plot_worst, '--', color=col_f, lw=lw_thin)
        ax2.set_ylim(ymin_fun, ymax_fun)
        ax2.tick_params('y', colors=col_f)
        # ax2.grid(color=col_f, lw=lw_grid/2, ls=ls_grid)
        ax2.set_zorder(0)

        if ncol == Ncols - 1:
            ax2.set_ylabel(ylabel_fun, fontsize=fs_fun, color=col_f)
        else:
            ax2.set_yticks([])
            ax2.set_yticklabels([])

    plt.tight_layout()
    plt.savefig(filename)
    print('Saved plot to %s' % filename)
    plt.close()

def plot_all_simplices(record, filename, axes=(0, 1), figsize=None, vmin=None,
        vmax=None, parameter_names=None):
    a0, a1 = axes
    all_simplices = record['all_simplices']
    K, M, N = np.shape(all_simplices)
    fvals = record['all_fun_values']
    if vmin is None:
        vmin = np.nanmin(record['best_fun_values'])
    if vmax is None:
        vmax = np.nanmax(record['worst_fun_values'])
    vinc = (vmax - vmin) / 50.
    vlevels = np.arange(vmin, vmax + vinc / 2., vinc)
    cmap = plt.get_cmap('gist_rainbow_r', 2**10)

    if figsize is None:
        figsize = (80/25.4, 80/25.4)
    plt.figure(figsize=figsize)

    for k in range(K):
        f = 0.9 * (1 - 1. * k/K)
        col_edge = 'k'

        x = all_simplices[k, :, a0]
        y = all_simplices[k, :, a1]
        f = fvals[k]

        for m1 in range(M):
            for m2 in range(m1 + 1, M):
                add_color_gradient_line((x[m1], x[m2]), (y[m1], y[m2]),
                        (f[m1], f[m2]), vmin, vmax, cmap)
            # plt.plot(x[m1], y[m1], 'k.')

    # axes labels
    if parameter_names is None:
        xlabel = 'parameter %i' % a0 
        ylabel = 'parameter %i' % a1 
    else:
        xlabel = parameter_names[a0]
        ylabel = parameter_names[a1]

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()

    ###################################################
    # COLORBAR                                        #
    ###################################################
    ax = plt.gca()
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    cb = plt.colorbar(sm)
    cb.set_label('function value')

    plt.tight_layout()
    print('Saving to file %s ...' % filename)
    plt.savefig(filename)
    print('Done.')
    plt.close()

###################################################
# HELPERS                                         #
###################################################
def add_color_gradient_line(xx, yy, ff, vmin, vmax, cmap, N=3, lw=1):
    # interpolate
    x1, x2 = xx
    y1, y2 = yy
    f1, f2 = ff
    inc = 1./N
    a = np.arange(0, 1 + inc/2., inc)
    x = x1 + a * (x2 - x1)
    y = y1 + a * (y2 - y1)
    f = f1 + a * (f2 - f1)
    for n in range(len(a)-1):
        fmean = np.mean((f[n], f[n+1]))
        cval = (fmean - vmin) / (vmax - vmin)
        color = cmap(cval)
        plt.plot((x[n], x[n+1]), (y[n], y[n+1]), '-', color=color, lw=lw)



###################################################
# TESTING                                         #
###################################################
if __name__ == '__main__':
    filename = '/home/anhaeus/tmp/nmplot.eps'
    plot_all_parameters(record, filename, fun_log=True)
    filename = '/home/anhaeus/tmp/nm_simplices.eps'
    plot_all_simplices(record, filename)
