import collections
import errno
import os

def flatten(l):
    """Flatten a list or tuple which may contain lists or tuples"""
    """For example, list(flatten([1, [2, 3, [4, 5]]]) returns [1, 2, 3, 4, 5]"""
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def make_sure_path_exists(path):
    """Create path, but don't complain if it already exists"""
    """Note in python 3.x (x>2), this is os.makedirs(path, exist_ok=True)"""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise