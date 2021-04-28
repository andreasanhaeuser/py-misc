"""A type that measures elapsed time."""

# standard modules
import datetime as dt

# misc
from misc.text.string_utils import human_format

################################################################
# Main                                                         #
################################################################
class Timer(object):
    """Initialize stopped timer.

        Parameters
        ----------
        <none>

        Methods
        -------
        start() : (re-)start counting
        stop() : stop counting
        get() : show elapsed time
        reset() : set to 0. and stop

        Methods may be concatenated:
        >>> timer = Timer().start()
        >>> timer.stop().show().reset().start()

        Examples
        --------
        >>> from time import sleep

        >>> # initialize and start timer
        >>> timer = Timer().start()
        >>> sleep(0.3)

        >>> # show status while running
        >>> print(timer)
        >>> sleep(0.7)

        >>> # show status after stopping
        >>> timer.stop()
        >>> print(timer)

        >>> # show status of stopped timer later
        >>> sleep(0.5)
        >>> print(timer)
    """
    # class variables
    autostart = False
    autoshow = False

    ############################################################
    # INIT                                                     #
    ############################################################
    def __init__(self, name=None):
        """Return self.
            
            A Timer has two states:
            - running
            - stopped
            When stopped, its time count (field `elapsed`) may be non-zero if
            it has been running previously.

            Fields
            ------
            elapsed : float
                (s) elapsed time while running
            first_started : dt.datetime
                time when it was first started
            last_started : None of dt.datetime
                time when it was last started; None if currently stopped.
        """
        self.elapsed = dt.timedelta()
        self.first_started = dt.datetime.now()
        self.last_started = None
        self.name = name
        if self.autostart:
            self.start()

    ############################################################
    # CONTEXT MANAGEMENT                                       #
    ############################################################
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Clean up."""
        if self.autoshow:
            self.show()

    def __del__(self):
        pass

    ############################################################
    # USER FUNCTIONS                                           #
    ############################################################
    def __str__(self):
        """Return elapsed time as human readable string."""
        secs = self.get('s')
        return str(human_format(secs)) + 's'

    def __repr__(self):
        """Return a str."""
        elapsed = self.get('s')
        if self.is_running():
            state = 'running'
        else:
            state = 'stopped'

        if self.name is None:
            name_str = ''
        else:
            name_str = ' "%s"' % self.name
        return '%s Timer%s at %f seconds' % (state, name_str, elapsed)

    def reset(self):
        """Return a timer with zero elapsed time."""
        self.elapsed = dt.timedelta()
        now = dt.datetime.now()
        self.first_started = now
        if self.is_running():
            self.last_started = now
        return self

    def start(self, ignore_running=False):
        """Start the timer and return self.
            
            Parameters
            ----------
            ignore_running : bool, optional
                (default: False) If False, an error is thrown, if the timer is
                already running. If True, an already running Timer is left
                unchanged.

            Throws
            ------
            AlreadyRunningError
        """
        if self.is_running():
            if ignore_running:
                return self
            else:
                message = 'Attempt to start already running Timer.'
                raise AlreadyRunningError(message)

        # regular case
        self.last_started = dt.datetime.now()
        return self

    def stop(self, ignore_stopped=False):
        """Stop the timer and return self.
            
            Parameters
            ----------
            ignore_stopped : bool, optional
                (default: False) If False, an error is thrown, if the timer is
                already stopped. If True, an already stopped Timer is left
                unchanged.

            Throws
            ------
            AlreadyStoppedError
        """
        if not self.is_stopped():
            pass

        elif ignore_stopped:
            return self

        else:
            raise AlreadyStoppedError('Attempt to stop already stopped Timer.')

        # regular case
        now = dt.datetime.now()
        diff = now - self.last_started
        self.elapsed += diff
        self.last_started = None
        return self

    def is_running(self):
        """Return a bool."""
        return not self.is_stopped()

    def is_stopped(self):
        """Return a bool."""
        return self.last_started is None

    def get(self, type='t'):
        """Return elapsed time in seconds as float."""
        if self.is_running():
            self.stop()
            self.start()

        if type[:1] == 't':
            return self.elapsed
        elif type[:1] == 's':
            return self.elapsed.total_seconds()
        else:
            raise ValueError(
                    'Output type must be "(s)econds" or "(t)imedelta".'
                    )

    def get_start(self):
        """Return a datetime.datetime."""
        return self.first_started

    def show(self, verbose=True):
        """Print elapsed time and return self."""
        if verbose:
            if self.name is None:
                name_str = ''
            else:
                name_str = ' [%s]' % self.name
            print(str(self.get()) + name_str)
        else:
            print(self.get())
        return self

    def set(self, elapsed):
        """Set elapsed time and return self.

            Running/stopping state remains unchanged.

            Parameters
            ----------
            elapsed : float or datetime.timedelta
        """
        if not isinstance(elapsed, dt.timedelta):
            elapsed = dt.timedelta(seconds=float(elapsed))

        self.elapsed = elapsed
        if self.is_running():
            self.last_started = dt.datetime.now()
        return self


################################################################
# Exceptions                                                   #
################################################################
class TimerError(Exception):
    pass

class AlreadyRunningError(TimerError):
    pass

class AlreadyStoppedError(TimerError):
    pass
