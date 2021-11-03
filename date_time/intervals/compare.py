def contains(self, other):
    if isinstance(self, Interval):
        return interval_contains(other)
    if isinstance(self, DaytimePeriod):
        return daytime_period_contains(other)
    if isinstance(self, Season):
        return season_contains(other)
    raise TypeError(type(self))

################################################################
# x contains                                                   #
################################################################
def interval_contains(other):
    assert isinstance(self, Interval)
    if isinstance(other, DaytimePeriod):
        return interval_contains_daytime_period(self, other)
    if isinstance(other, Season):
        return interval_contains_season(self, other)
    raise TypeError(type(other))

def daytime_period_contains(other):
    assert isinstance(self, DaytimePeriod)
    if isinstance(other, Interval):
        return daytime_period_contains_interval(self, other)
    if isinstance(other, Season):
        return daytime_period_contains_season(self, other)
    raise TypeError(type(other))

def season_contains(other):
    assert isinstance(self, Season)
    if isinstance(other, Interval):
        return season_contains_interval(self, other)
    if isinstance(other, DaytimePeriod):
        return season_contains_daytime_period(self, other)
    raise TypeError(type(other))

################################################################
# x contains y                                                 #
################################################################
def interval_contains_daytime_period(self, other):
    assert isinstance(self, Interval)
    assert isinstance(other, DaytimePeriod)
    raise NotImplementedError()

def interval_contains_season(self, other):
    assert isinstance(self, Interval)
    assert isinstance(other, Season)
    raise NotImplementedError()

def daytime_period_contains_interval(self, other):
    assert isinstance(self, DaytimePeriod)
    assert isinstance(other, Interval)
    raise NotImplementedError()

def daytime_period_contains_season(self, other):
    assert isinstance(self, DaytimePeriod)
    assert isinstance(other, Season)
    raise NotImplementedError()

def season_contains_interval(self, other):
    assert isinstance(self, Season)
    assert isinstance(other, Interval)
    raise NotImplementedError()

def season_contains_daytime_period(self, other):
    assert isinstance(self, Season)
    assert isinstance(other, DaytimePeriod)
    raise NotImplementedError()
