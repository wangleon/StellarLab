from itertools import tee
import numpy as np

def get_edge_bin(array):
    '''
    Detect the edge indcies of a binary 1-D array.

    Args:
        array (:class:`numpy.array`): A list or Numpy 1d array, with binary
            (0/1) or boolean (True/False) values.
    Returns:
        list: A list containing starting and ending indices of the non-zero
            blocks.
    Examples:
        .. code-block:: python

            >>> a = [0,1,1,0,0,0,1,0,1]
            >>> get_edge_bin(a)
            [(1, 3), (6, 7), (8, 9)]
            >>> b = [True, False, True, True, False, False]
            >>> get_edge_bin(b)
            [(0, 1), (2, 4)]

    '''
    array1 = np.int64(array)
    array1 = np.insert(array1, 0, 0)
    array1 = np.append(array1, 0)
    tmp = array1 - np.roll(array1, 1)
    i1_lst = np.nonzero(tmp == 1)[0] - 1
    i2_lst = np.nonzero(tmp ==-1)[0] - 1
    return list(zip(i1_lst, i2_lst))

def get_local_minima(x, window=None):
    '''
    Get the local minima of a 1-d array in a window.
    
    Args:
        x (:class:`numpy.array`): A list or Numpy 1d array.
        window (int): An odd integer as the length of searching window.
    Returns:
        tuple: A tuple containing:

            * **index** (:class:`numpy.array`): a numpy 1d array containing the indices of all local minima.
            * **x[index]** (:class:`numpy.array`): a numpy 1d array containing the values of all local minima.

    '''
    x = np.array(x)
    dif = np.diff(x)
    ind = dif > 0
    tmp = np.logical_xor(ind, np.roll(ind,1))
    idx = np.logical_and(tmp,ind)
    index = np.where(idx)[0]
    if window is None:
        return index, x[index]
    else:
        # window must be an odd integer
        if window%2 != 1:
            raise ValueError
        halfwin = int(round((window-1)/2.))
        index_lst = []
        for i in index:
            i1 = max(0, i-halfwin)
            i2 = min(i+halfwin+1, x.size)
            if i == x[i1:i2].argmin() + i1:
                index_lst.append(i)
        index_lst = np.array(index_lst)
        return index_lst, x[index_lst]

def pairwise(iterable_array):
    a, b = tee(iterable_array)
    next(b, None)
    return zip(a, b)
