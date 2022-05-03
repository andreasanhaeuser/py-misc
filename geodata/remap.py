#!/usr/bin/python3

# PyPI
import numpy as np
from scipy import interpolate

# misc
from misc.date_time.timer import Timer

def gdf_to_grid(
        gdf, value_name, target_x, target_y, method='linear', fill_value=0,
        dist_fwhm=5., max_dist=25.,
        ):
    """Return a 2d-array.

        Parameters
        ----------
        gdf : GeoDataFrame
        value_name : str
            field of `gdf` to be interpolated
        target_x : nd-array
            x-coordinates of target grid
        target_y : nd-array, same size as `target_x`
            y-coordinates of target grid
        method : 'linear' or 'cubic' or ...
            as in scipy.interpolate.griddata
        fill_value : object
            as in scipy.interpolate.griddata

        Returns
        -------
        remapped : nd-array, same size as `target_x`
            remapped data
    """
    ############################################################
    # input check                                              #
    ############################################################
    shape_x = np.shape(target_x)
    shape_y = np.shape(target_y)
    if shape_x != shape_y:
        i_need = "target coords must have same shape"
        this_is = "got %s and %s" % (shape_x, shape_y)
        message = (i_need + ', but ' + this_is).capitalize()
        raise ValueError(message)

    ############################################################
    # main                                                     #
    ############################################################
    # source values
    source_values = np.array(gdf[value_name])

    # source coordinates
    source_x = np.array(gdf.geometry.x)
    source_y = np.array(gdf.geometry.y)

    # cast coordinates to suitable format
    source_xy = np.array([source_x, source_y]).T
    target_xy = np.array([target_x, target_y]).T

    # interpolate
    if method in ('linear', 'cubic'):
        remapped = interpolate.griddata(
            source_xy,
            source_values,
            target_xy,
            method=method,
            fill_value=fill_value,
            ).T
    elif method == 'rbf':
        args = (source_x, source_y, source_values)
        kwargs = {
                'epsilon' : 1,
                'function' : 'cubic',
                }
        interpolator = interpolate.Rbf(*args, **kwargs)
        remapped = interpolator(target_x, target_y)
    elif method in ('calpuff_devel', 'crea'):
        N = len(source_values)
        weighted_sum = np.zeros_like(target_x)
        weights = np.zeros_like(target_x)
        # zmin = np.nanmax(source_values) / 1000.
        # dist_fwhm2 = dist_fwhm**2
        maxdist2 = max_dist**2
        for n in range(N):
            z = source_values[n]

            if not np.isfinite(z):
                continue

            if n % 1000 == 0:
                print('- interpolate: %1.1f%% (%i/%i)' % (n/N*100, n, N))

            x = source_x[n]
            y = source_y[n]

            dx = x - target_x
            dy = y - target_y
            dist2 = dx**2 + dy**2

            idx = dist2 <= maxdist2

            dist2 = dist2[idx]

            idx0 = dist2 == 0
            if np.any(idx0):
                mindist = np.min(dist2)
                dist2[idx0] = mindist2 #* 1e-6


            weight = 1/dist2
            # weight = 1/dist2**2
            # weight = np.exp(- dist2/dist_fwhm2)
            # weight[np.isnan(weight)] = 0.
            # weight[dist2>maxdist2] = 0.

            weighted_sum[idx] += (z * weight)
            weights[idx] += weight

        weights[weighted_sum==0] = 1.
        remapped = weighted_sum / weights

    ############################################################
    # output check                                             #
    ############################################################
    assert np.shape(remapped) == shape_x

    return remapped


################################################################
# devel                                                        #
################################################################
def load(setup, data):
    time = dt.datetime(2016, 12, 31, 12)
    setup['times_current'] = time

    clusters = setup['clusters']
    scenario = setup['scenario']
    loader = calpuff.get_spatial_one_instant

    key = 'pollution_instant_raw' 
    if key in data:
        raw_all = data[key]
    else:
        raw_all = {}

    kwargs = {
            'project_name' : setup['project_name'],
            'substance' : setup['substance'],
            'setup' : setup,
            'time' : setup['times_current'],
            }

    for n, cluster in enumerate(clusters):
        run_name = '%s%s' % (scenario, cluster)
        kwargs['run_name'] = run_name
        if run_name in raw_all:
            continue

        print('- load  %s' % run_name)
        raw = loader(**kwargs)
        raw_all[run_name] = raw

    data[key] = raw_all
    return setup, data

def remap(setup, data):
    grid_xm = setup['grid_x_mesh']
    grid_ym = setup['grid_y_mesh']
    method = 'devel'
    # method = 'linear'

    clusters = setup['clusters']
    scenario = setup['scenario']
    which = 'instant'

    cluster_combiner = get_cluster_combiner(setup, which)

    for n, cluster in enumerate(clusters):
        run_name = '%s%s' % (scenario, cluster)

        # remap
        print('- remap %s' % run_name)
        gdf = data['pollution_instant_raw'][run_name]
        value_name = calpuff_utils.get_data_column_name(gdf)
        args = gdf, value_name, grid_xm, grid_ym, method
        remapped_one = gdf_to_grid(*args)
        
        # combine
        # --------------------------------------------
        if n == 0:
            remapped_total = np.zeros_like(remapped_one)

        to_combine = (remapped_total, remapped_one)
        remapped_total = cluster_combiner(to_combine, axis=0)
        # --------------------------------------------

    data['pollution_instant'] = remapped_total
    return setup, data

def plot(setup, data):
    set_up_figure(setup, data)
    plot_pollution_instant(setup, data)
    return setup, data


################################################################
# helpers                                                      #
################################################################
def get_cluster_combiner(setup, which):
    """Return appropriate function to combine pollution from clusters."""
    source = setup['pollution_source']

    # Calpuff =============================================
    if source == 'calpuff':
        if which == 'accumulated':
            return np.nansum
        if which == 'instant':
            return np.nansum
        if which == 'statistics':
            return get_statistics_combiner(setup)
        raise NotImplementedError('Unknown pollution type: %s' % which)

    # Tropomi =============================================
    if source == 'tropomi':
        return np.nanmax

    # unknwon source ======================================
    raise NotImplementedError('Unknown pollution source: %s' % source)


################################################################
# main                                                         #
################################################################
def main(setup_file=None, data=None):
    args = parse_command_line_arguments(setup_file)
    setup = load_setup_file(args)
    process_setup(setup)

    if data is None:
        data = {}
    if 'artists' not in data:
        data['artists'] = {}

    timer = Timer()

    timer.start(False)
    load(setup, data)
    timer.stop(False)
    timer.show()
    timer.start(False)
    remap(setup, data)
    timer.stop(False)
    timer.show()

    plot(setup, data)

    return setup, data


if __name__ == '__main__':
    # standard
    import os
    import datetime as dt

    # PyPI
    import matplotlib.pyplot as plt

    # GAPU
    from gapu.vis.utils.figure import set_up_figure
    from gapu.vis.utils.pollution import plot_instant
    from gapu.in_out.pollution import calpuff
    from gapu.in_out.pollution import calpuff_utils

    # setup
    from gapu.vis.pollution_map import helpers
    parse_command_line_arguments = helpers.setup.parse_args
    load_setup_file = helpers.setup.load_setup_file
    process_setup = helpers.setup.process

    plt.close()
    dirname = '~/gp/gapu_code/vis/pollution_map/korea/coal_free'
    # basename = 'timeseries.txt'
    basename = 'devel.txt'
    setup_file = os.path.expanduser(dirname + '/' + basename)

    if 'data' not in globals():
        data = None

    setup, data = main(data=data, setup_file=setup_file)
    plt.show()
