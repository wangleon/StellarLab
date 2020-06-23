import collections
import functools

class memoized(object):
    """
    Memoized class for Decorator.
    
    Notes
    ------
    This class caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    From https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize

    Examples
    --------

    .. code-block:: python

        from stella.utils.memoize import memoized

        @memoized
        def fibonacci(n):
            "Return the nth fibonacci number."
            if n in (0, 1):
                return n
            return fibonacci(n-1) + fibonacci(n-2)

        print fibonacci(12)

    """
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value
    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__
    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)
