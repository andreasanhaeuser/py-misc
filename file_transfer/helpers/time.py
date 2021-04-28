def get_utc_offset():
    utc_offset_exact = dt.datetime.now() - dt.datetime.utcnow()
    secs = utc_offset_exact.total_seconds()
    hours = round(secs/3600.)
    utc_offset = dt.timedelta(hours=hours)
    return utc_offset

def get_times_from_hours_key(jobs, hours_key):
    """Return a list of datetime.datetime."""
    assert hours_key in jobs
    hours_all = jobs[hours_key]
    times = []
    for n in range(len(hours_all)):
        hours = float(hours_all[n])
        time = dt.datetime.now() - dt.timedelta(hours=hours)
        times.append(time)
    return times

def get_times_copy(jobs):
    return get_times_from_hours_key(jobs, 'hours_copy')

def get_times_delete(jobs):
    return get_times_from_hours_key(jobs, 'hours_delete')

def get_time_file(filename):
    """Return last file modification time as datetime.datetime."""
    if not os.path.isfile(filename):
        return None
    utc_offset = get_utc_offset()
    mtime_utc_sec = os.path.getmtime(filename)
    mtime_utc = dt_utils.seconds_to_datetime(mtime_utc_sec)
    mtime = mtime_utc + utc_offset
    return mtime
