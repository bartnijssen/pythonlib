import collections

def flatten(l):
    """Flatten a list or tuple which may contain lists or tuples"""
    """For example, list(flatten([1, [2, 3, [4, 5]]]) returns [1, 2, 3, 4, 5]"""
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el