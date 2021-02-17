import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils
try:
    import KratosMultiphysics.ShallowWaterApplication.output.convergence_output_process as convergence_output
    import h5py
    h5py_available = True
except ImportError:
    h5py_available = False

@KratosUnittest.skipIfApplicationsNotAvailable("StatisticsApplication")
class TestConvergenceOutputProcess(KratosUnittest.TestCase):

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("output_file")

    def _CheckValues(self, file_name, dset_name, run_id, variables_list, values_list):
        h5file = h5py.File(file_name, 'r')
        dset = h5file[dset_name]
        for (variable, value) in zip (variables_list, values_list):
            self.assertAlmostEqual(dset[variable.Name()][run_id], value)

    def _CheckTwoDatasets(self, file_name, dset_name_a, dset_name_b, attribute_key, value_a, value_b):
        h5file = h5py.File(file_name, 'r')

        self.assertEqual(len(h5file.keys()), 2)

        self.assertTrue(list(h5file.keys())[0] == dset_name_a)
        self.assertTrue(list(h5file.keys())[1] == dset_name_b)

        self.assertEqual(h5file[dset_name_a].attrs[attribute_key], value_a)
        self.assertEqual(h5file[dset_name_b].attrs[attribute_key], value_b)

    def testSingleConvergenceOutputProcessCheck(self):
        if not h5py_available:
            self.skipTest('h5py is not available')
        settings = KM.Parameters('''{
            "Parameters"         : {
                "model_part_name"            : "model_part",
                "file_name"                  : "output_file",
                "analysis_label"             : "label",
                "analysis_attributes"        : {
                    "density"                     : 1.0
                },
                "convergence_variables_list" : ["ERROR_RATIO","NODAL_ERROR"]
            }
        }''')
        variables = [KM.ERROR_RATIO, KM.NODAL_ERROR]

        values = [1e-3, 1e-4]
        self._AuxiliaryAnalysis(settings, variables, values)
        self._CheckValues(settings["Parameters"]["file_name"].GetString() + ".hdf5", "analysis_000", 0, variables, values)

        values = [2e-3, 2e-4]
        self._AuxiliaryAnalysis(settings, variables, values)
        self._CheckValues(settings["Parameters"]["file_name"].GetString() + ".hdf5", "analysis_000", 1, variables, values)

        kratos_utils.DeleteFileIfExisting(settings["Parameters"]["file_name"].GetString() + ".hdf5")

    def testMultipleConvergenceOutputProcessCheck(self):
        if not h5py_available:
            self.skipTest('h5py is not available')
        settings = KM.Parameters('''{
            "Parameters"         : {
                "model_part_name"            : "model_part",
                "file_name"                  : "output_file",
                "analysis_label"             : "label",
                "analysis_attributes"        : {
                    "density"                    : 1.0,
                    "dummy_flag"                 : true,
                    "dummy_string"               : "string"
                },
                "convergence_variables_list" : ["ERROR_RATIO"]
            }
        }''')
        variables = [KM.ERROR_RATIO]
        values = [1e-3]
        self._AuxiliaryAnalysis(settings, variables, values)

        settings["Parameters"]["analysis_attributes"]["density"].SetDouble(13.6)
        self._AuxiliaryAnalysis(settings, variables, values)

        self._CheckTwoDatasets(settings["Parameters"]["file_name"].GetString() + ".hdf5", "analysis_000", "analysis_001", "density", 1.0, 13.6)

        kratos_utils.DeleteFileIfExisting(settings["Parameters"]["file_name"].GetString() + ".hdf5")

    @staticmethod
    def _AuxiliaryAnalysis(settings, variables_list, values_list):
        m = KM.Model()
        mp = m.CreateModelPart("model_part")

        mp.AddNodalSolutionStepVariable(KM.DISTANCE)

        mp.ProcessInfo[KM.DOMAIN_SIZE] = 2
        mp.ProcessInfo[KM.GRAVITY_Z] = 9.81

        KM.ModelPartIO("processes_tests/model_part").ReadModelPart(mp)

        process = convergence_output.Factory(settings, m)
        process.Check()

        step = 0
        time = 0.0
        end_time = 2.0
        time_step = 0.1

        mp.ProcessInfo[KM.STEP] = step
        mp.ProcessInfo[KM.TIME] = time

        process.ExecuteBeforeSolutionLoop()

        while time < end_time:
            step += 1
            time += time_step

            mp.ProcessInfo[KM.STEP] = step
            mp.ProcessInfo[KM.TIME] = time

            for (variable, value) in zip(variables_list, values_list):
                KM.VariableUtils().SetNonHistoricalVariable(variable, value, mp.Nodes)

        process.ExecuteFinalize()

if __name__ == '__main__':
    KratosUnittest.main()
