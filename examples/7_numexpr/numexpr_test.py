import numpy as np
import numexpr as ne
import time


def compare():
    a = np.random.rand(1000000)
    b = np.random.rand(1000000)

    start = time.time()
    result = 2 * a + 3 * b
    end = time.time()
    print('numpy:', end - start)

    start = time.time()
    result = ne.evaluate("2*a + 3*b")
    end = time.time()
    print('numexpr:', end - start)

def calcualte():
    v = np.random.rand(1000, 1)
    m = np.random.rand(1000, 1)
    m1 = m.reshape(1, len(m))
    m2 = m.reshape(len(m), 1)

    v1 = v.reshape(1, len(v))
    v2 = v.reshape(len(v), 1)

    start = time.time()
    res = (np.log(v1) - np.log(v2)) /2 + ((m1 - m2)**2 + v2 ) / ( 2 * v1 ) - 0.5
    end = time.time()
    print('numpy:', end - start)

    start = time.time()
    res = ne.evaluate(
        '(log(v1) - log(v2))/2 + ( (m1 - m2)**2 + v2 ) / ( 2 * v1 ) - half',
        {
            'v1': v1,
            'v2': v2,
            'm1': m1,
            'm2': m2,
            # numexpr rules are that mixing with floats causes
            # conversion to float64
            # Note that this does not happen for integers
            'half': np.float32(0.5),
        })
    end = time.time()
    print('numexpr:', end - start)

if __name__ == '__main__':
    calcualte()
