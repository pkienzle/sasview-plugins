def profile(fn, *args, **kw):
    """
    Profile a function called with the given arguments.
    """
    import cProfile, pstats, os
    global call_result
    def call():
        global call_result
        call_result = fn(*args, **kw)
    datafile = 'profile.out'
    cProfile.runctx('call()', dict(call=call), {}, datafile)
    stats = pstats.Stats(datafile)
    #order='calls'
    order='cumulative'
    #order='pcalls'
    #order='time'
    stats.sort_stats(order)
    stats.print_stats()
    os.unlink(datafile)
    return call_result