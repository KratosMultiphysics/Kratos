import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

import numpy as np
import os
from copy import deepcopy


class TestMapperInterpolator(KratosUnittest.TestCase):
    def test_mapper_interpolator(self):
        parameter_file_name = os.path.join(os.path.dirname(__file__), 'test_interpolator.json')
        cs_data_structure = ImportDataStructure(parameter_file_name)
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = cs_data_structure.Parameters(parameter_file.read())
        par_mapper_0 = parameters['mapper']

        # test if directions are set correctly in __init__
        par_mapper = deepcopy(par_mapper_0)

        par_mapper['settings'].SetArray('directions', ['Z'])
        mapper = cs_tools.CreateInstance(par_mapper)
        self.assertListEqual(mapper.directions, ['Z0'])

        par_mapper['settings'].SetArray('directions', ['Z', 'Y', 'X'])
        mapper = cs_tools.CreateInstance(par_mapper)
        self.assertListEqual(mapper.directions, ['Z0', 'Y0', 'X0'])

        par_mapper['settings'].SetArray('directions', ['z'])
        self.assertRaises(ValueError, cs_tools.CreateInstance, par_mapper)

        par_mapper['settings'].SetArray('directions', ['Z0'])
        self.assertRaises(ValueError, cs_tools.CreateInstance, par_mapper)

        par_mapper['settings'].SetArray('directions', ['Z', 'Y', 'X', 'Z'])
        self.assertRaises(ValueError, cs_tools.CreateInstance, par_mapper)

        par_mapper['settings'].SetString('directions', 'Z')
        self.assertRaises(TypeError, cs_tools.CreateInstance, par_mapper)

        # test check_bounding_box method
        if True:
            par_mapper = deepcopy(par_mapper_0)

            # check 1D errors and warnings
            par_mapper['settings'].SetArray('directions', ['Z'])
            mapper = cs_tools.CreateInstance(par_mapper)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')

            mp_from.CreateNewNode(0, 0., 0., 0.)
            mp_from.CreateNewNode(1, 0., 0., 1.)

            mp_to.CreateNewNode(0, 0., 0., 0.)
            mp_to.CreateNewNode(1, 0., 0., 1.)

            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(2, 0., 0., 1.01)
            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(3, 0., 0., -.01)
            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(11, 0., 0., 1.1)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(12, 0., 0., 1.25)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(13, 0., 0., -.25)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(14, 0., 0., 2.)
            mp_to.CreateNewNode(15, 0., 0., -1.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            # check 2D errors and warnings
            par_mapper['settings'].SetArray('directions', ['Z', 'X'])
            mapper = cs_tools.CreateInstance(par_mapper)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')

            mp_from.CreateNewNode(0, 0., 0., 0.)
            mp_from.CreateNewNode(1, 1., 0., 1.)

            mp_to.CreateNewNode(0, 0., 0., 0.)
            mp_to.CreateNewNode(1, 1., 0., 1.)

            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(2, 1.01, 0., 1.01)
            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(3, -.01, 0., -.01)
            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(11, 1.1, 0., 1.1)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(12, 1.25, 0., 1.25)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(13, -.25, 0., -.25)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(14, 2., 0., 2.)
            mp_to.CreateNewNode(15, -1., 0., -1.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            # check 3D errors and warnings
            par_mapper['settings'].SetArray('directions', ['Z', 'X', 'Y'])
            mapper = cs_tools.CreateInstance(par_mapper)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')

            mp_from.CreateNewNode(0, 0., 0., 0.)
            mp_from.CreateNewNode(1, 1., 1., 1.)

            mp_to.CreateNewNode(0, 0., 0., 0.)
            mp_to.CreateNewNode(1, 1., 1., 1.)

            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(2, 1.01, 1.01, 1.01)
            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(3, -.01, -.01, -.01)
            mapper.Initialize(mp_from, mp_to)

            mp_to.CreateNewNode(11, 1.1, 1.1, 1.1)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(12, 1.25, 1.25, 1.25)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(13, -.25, -.25, -.25)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            mp_to.CreateNewNode(14, 2., 2., 2.)
            mp_to.CreateNewNode(15, -1., -1., -1.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            # check if method works for lines aligned with coordinate axes in 2D
            par_mapper['settings'].SetArray('directions', ['X', 'Y'])
            mapper = cs_tools.CreateInstance(par_mapper)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1., 0.)
            mp_to.CreateNewNode(1, 1.01, 1.01, 0.)
            mapper.Initialize(mp_from, mp_to)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.05, 0.)
            mp_to.CreateNewNode(1, 1., 1.05, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.25, 0.)
            mp_to.CreateNewNode(1, 1., 1.25, 0.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))
            # mapper.Initialize(mp_from, mp_to)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., .85, 0.)
            mp_to.CreateNewNode(1, 1., 1.15, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            # check if method works for planes aligned with coordinate axes in 3D
            par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])
            mapper = cs_tools.CreateInstance(par_mapper)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1., 0.)
            mp_to.CreateNewNode(1, 1.01, 1.01, 0.)
            mapper.Initialize(mp_from, mp_to)

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.05, 0.)
            mp_to.CreateNewNode(1, 1., 1.05, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., 1.25, 0.)
            mp_to.CreateNewNode(1, 1., 1.25, 0.)
            self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

            model = cs_data_structure.Model()
            mp_from = model.CreateModelPart('mp_from')
            mp_to = model.CreateModelPart('mp_to')
            mp_from.CreateNewNode(0, 0., 1., 0.)
            mp_from.CreateNewNode(1, 1., 1., 0.)
            mp_to.CreateNewNode(0, 0., .85, 0.)
            mp_to.CreateNewNode(1, 1., 1.15, 0.)
            self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

        # test check_duplicate_points method
        par_mapper = deepcopy(par_mapper_0)
        par_mapper['settings'].SetArray('directions', ['X', 'Y', 'Z'])

        mapper = cs_tools.CreateInstance(par_mapper)

        model = cs_data_structure.Model()
        mp_from = model.CreateModelPart('mp_from')
        mp_to = model.CreateModelPart('mp_to')
        mp_from.CreateNewNode(0, 0., 0., 0.)
        mp_from.CreateNewNode(1, 1., 0., 0.)
        mp_to.CreateNewNode(0, 0., 0., 0.)
        mp_to.CreateNewNode(1, 1., 0., 0.)

        mp_from.CreateNewNode(2, 1e-10, 0., 0.)
        self.assertRaises(Warning, mapper.Initialize, *(mp_from, mp_to))

        mp_from.CreateNewNode(3, 1e-14, 0., 0.)
        self.assertRaises(ValueError, mapper.Initialize, *(mp_from, mp_to))

        # to do: check tree? check __call__ method?


if __name__ == '__main__':
    KratosUnittest.main()
