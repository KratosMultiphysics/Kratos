from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.factories import convergence_criterion_factory


class TestConvergenceCriteria(KratosUnittest.TestCase):
    def test_abc(self):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
