import fileinput
from itertools import islice

def sliding_window(seq, n=3):
    """Returns an iterator for a sliding window over a given sequence.

    Each element is a tuple of length `n` containing the window.

    Parameters
    ----------
    seq : string
        The sequence to slide a window over
    n : int
        Size of the sliding window.

    Example
    -------

    s -> (s0, s1, ..., sn-1), (s1, s2, sn)

    Notes
    -----

    See also this stackoverflow function:
    https://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
    """

    it = iter(seq)
    result = tuple(islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + (elem,)
        yield result

