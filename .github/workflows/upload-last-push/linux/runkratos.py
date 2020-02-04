import scipy
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import matplotlib

import numpy
import sys
print("scipy version = ",scipy.version)
print("runkratos file ", sys.argv[1])
exec(open(sys.argv[1]).read())