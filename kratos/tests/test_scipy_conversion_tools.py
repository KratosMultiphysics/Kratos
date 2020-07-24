from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.scipy_conversion_tools
    import numpy as np
    missing_scipy = False
except ImportError as e:
    missing_scipy = True

class TestScipyConversionTools(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(missing_scipy,"Missing python libraries (scipy)")
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
