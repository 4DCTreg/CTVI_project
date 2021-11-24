import numpy as np
import array
def test():
    print('Hello, Matlab!')

def add(a, b):
    a = np.array(a,dtype='float')
    b = np.array(b,dtype='float')
    c = a + b
    d = array.array('d',c)
    e = a - b
    e = array.array('d', e)

    return d, e