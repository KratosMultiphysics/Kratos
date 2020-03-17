import numpy as np
import scipy
import scipy.sparse
np.set_printoptions(linewidth=150)
from scipy.optimize import minimize
import scipy.sparse.linalg as sla
import matplotlib.pyplot as plt

#adapted from consistent_aat.ipynb

columns = True

_a = np.matrix([[3., 2., 1., 0., 0., 0., 0., 0., 0., 0.],
                [2., 3., 2., 1., 0., 0., 0., 0., 0., 0.],
                [1., 2., 3., 2., 1., 0., 0., 0., 0., 0.],
                [0., 1., 2., 3., 2., 1., 0., 0., 0., 0.],
                [0., 0., 1., 2., 3., 2., 1., 0., 0., 0.],
                [0., 0., 0., 1., 2., 3., 2., 1., 0., 0.],
                [0., 0., 0., 0., 1., 2., 3., 2., 1., 0.],
                [0., 0., 0., 0., 0., 1., 2., 3., 2., 1.],
                [0., 0., 0., 0., 0., 0., 1., 2., 3., 2.],
                [0., 0., 0., 0., 0., 0., 0., 1., 2., 3.]])

_a = scipy.sparse.csr_matrix(_a)

for i in range(_a.shape[0]):
  _a[i,:]/=_a[i,:].sum()

m = None
a = None
if columns:
  x = np.array([1.]*_a.shape[0]*2)
else:
  x = np.array([1.]*_a.shape[0])

if columns:
  bounds = [(0.1,None)]*_a.shape[0]*2
else:
  bounds = [(0.1,None)]*_a.shape[0]

def test(matrix, mod_matrix, vec):
  print(vec)
  r1 = matrix @ vec
  r2 = mod_matrix @ vec
  print(r1)
  print(r2)
  coords = np.linspace(1, r1.shape[0], num=r1.shape[0])
  fig = plt.figure()
  plt.plot(coords, r1, label="1")
  plt.plot(coords, r2, label="2")
  plt.legend()
  plt.show()

def _MakeConsistent_NR(matrix):
    from newton_raphson import newton_raphson_solve

    x = np.array([1.]*matrix.shape[0])

    def calculate_system(_x):
        _m = matrix.copy()
        a = _m.dot(scipy.sparse.diags(_x))

        m = a.dot(a.transpose())

        RHS = np.asarray((m.sum(1)-1).flatten())[0]

        LHS = matrix.copy()
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[0]):
                v = (_m[i,j] * _m.getcol(j)).sum() * 2 * _x[i]
                if v != 0:
                    LHS[i,j] = v

        return LHS, RHS

    x, n_iter = newton_raphson_solve(calculate_system, x, 100, 1e-6)

    _m = matrix.copy()
    a = _m.dot(scipy.sparse.diags(x))
    return a.dot(a.transpose())

def _MakeConsistent(matrix):
    LHS = matrix.copy()
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            v = (matrix[i,j] * matrix.getcol(j)).sum()
            if v != 0:
                LHS[i,j] = v

    x = sla.spsolve(LHS, [1]*matrix.shape[0])

    x = np.sqrt(x)

    _m = matrix.copy()
    a = _m.dot(scipy.sparse.diags(x))
    return a.dot(a.transpose())

m = _MakeConsistent(_a)

test1 = np.array([1.]*_a.shape[0])
test2 = np.array([0.]*_a.shape[0])
test2[3] = 1.
test3 = np.linspace(-1, 1, num=_a.shape[0])
test4 = np.array([x if x > 0 else -x for x in test3 ])
_m = _a @ _a.transpose()
print("test1")
test(_m, m, test1)
print("test2")
test(_m, m, test2)
print("test3")
test(_m, m, test3)
print("test4")
test(_m, m, test4)