from time import sleep

def call(
        f, *args, n=3, ignore_error=True, error_return_value=None,
        report_error=True, reporter=print, sleep_time=0, **kwargs
        ):
    """Try to call function for n times.

        If the function call raises an error, it is tried again; for `n` times.

        The intended purpose is for file and network operations, which may
        occasionally fail and then work just fine if you simply try it again
        with the same arguments.
    
        Parameters
        ----------
        f : callable
            the actual function to be called
        *args : non-keyword arguments for `f`
        **kwargs : keyword arguments for `f`
        n : int, optional
            (default: 3) number of times to try to call the function
        ignore_error : bool, optional
            (default: True) if this is False, it effectively disables the
            repitition (intended for debugging purposes)
        error_return_value : object, optional
            (default: None) if the function call remains unsuccessful after `n`
            trials, this value is returned
        report_error : bool, optional
            (default: True) report errors on function calls, even if they are
            handled.
        reporter : callable, optional
            (default: print [built-in]) a function that takes str as argument.
            The default action is to simply print it to STDOUT, but other
            strategies can be though of (eg. write it to a log file)

        Returns
        -------
        On success:
            Whatever `f` should return
        Otherwise:
            `error_return_value`

        Raises
        ------
        Exception
    """
    ntry = 0
    assert isinstance(n, int)
    assert n > 0
    assert callable(f)

    success = False
    while ntry < n:
        try:
            result = f(*args, **kwargs)
            success = True
            break
        except KeyboardInterrupt as ki:
            raise ki
        except Exception as exception:
            ntry += 1
            if ntry >= n:
                if report_error:
                    reporter('ERROR: %s' % exception)
                if not ignore_error:
                    raise exception
        sleep(sleep_time)

    if success:
        return result

    if ignore_error:
        return error_return_value

    raise Exception('Function execution failed %i times' % n)
