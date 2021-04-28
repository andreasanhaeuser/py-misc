def get_filename_out(filename_in):
    fn = filename_in.replace('/mnt/', '/data/')
    fn = fn.replace('/data/air_pollution_unit/projects/', '/data/')
    return prefix_bucket + fn
