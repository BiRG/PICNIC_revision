"""
This module provides splitBy function which splits a sequence into two
subsequences, in single pass and preserving order.

Some tests:

>>> def odd(x): return x%2 != 0
>>> odds, evens = splitBy(odd, range(10))
>>> odds.next()
1
>>> odds.next()
3
>>> list(evens)
[0, 2, 4, 6, 8]
>>> list(odds)
[5, 7, 9]

"""

from collections import deque

__all__ = [ "splitBy" ]
__author__ = "Sergey Astanin"
__license__ = "BSD3"
__version__ = "0.1"

class SplitSeq:
    """
    Lazily process a sequence in single pass and split into two.

    Computes both output sequences even if only one of them is consumed.
    """
    def __init__(self, condition, sequence):
        self.cond = condition
        self.goods = deque([])
        self.bads = deque([])
        self.seq = iter(sequence)
    def getNext(self, getGood=True):
        if getGood:
            these, those, cond = self.goods, self.bads, self.cond
        else:
            these, those, cond = self.bads, self.goods, lambda x: not self.cond(x)
        if these:
            return these.popleft()
        else:
            while 1: # exit on StopIteration
                n = self.seq.next()
                if cond(n):
                    return n
                else:
                    those.append(n)

def splitBy(condition, sequence):
    """
    Split a sequence into two subsequences, in single-pass and preserving order.

    Arguments:

    condition   a function; if condition is None, split true and false items
    sequence    an iterable object

    Return a pair of generators (seq_true, seq_false). The first one
    builds a subsequence for which the condition holds, the second one
    builds a subsequence for which the condition doesn't hold.

    As the function works in single pass, it leads to build-up of both
    subsequences even if only one of them is consumed.
    """
    cond = condition if condition else bool  # evaluate as bool if condition == None
    ss = SplitSeq(cond, sequence)
    def goods():
        while 1:
            yield ss.getNext(getGood=True)
    def bads():
        while 1:
            yield ss.getNext(getGood=False)
    return goods(), bads()

if __name__ == "__main__":
    import doctest
    doctest.testmod()