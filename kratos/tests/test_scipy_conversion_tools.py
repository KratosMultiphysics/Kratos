import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import numpy as np

try:
    import scipy
    import scipy.sparse
    import KratosMultiphysics.scipy_conversion_tools
    scipy_available = True
except ImportError:
    scipy_available = False

class TestScipyConversionTools(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(not scipy_available,"Missing python libraries (scipy)")
    def test_scipy_conversion_tools(self):
        # Create a kratos matrix
        A = KratosMultiphysics.CompressedMatrix(3,3)
        A[0,0] = 1.0
        A[1,1] = 2.0
        A[2,2] = 3.0
        A[0,2] = 4.0

        # Convert it to scipy. Note that Ascipy is A COPY of A
        Ascipy = KratosMultiphysics.scipy_conversion_tools.to_csr(A)

        for i, j in np.nditer(Ascipy.nonzero()):
            self.assertAlmostEqual(A[int(i), int(j)], Ascipy[int(i), int(j)], 1e-12)

if __name__ == '__main__':
    KratosUnittest.main()
