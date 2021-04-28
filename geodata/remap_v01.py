import numpy as np
from scipy import interpolate

def gdf_to_grid(
        gdf, value_name, target_x, target_y, method='linear', fill_value=0,
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
        raise ValueError(i_need + ', but ' + this_is)

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
    remapped = interpolate.griddata(
        source_xy,
        source_values,
        target_xy,
        method=method,
        fill_value=fill_value,
        ).T

    ############################################################
    # output check                                             #
    ############################################################
    assert np.shape(remapped) == shape_x

    return remapped
