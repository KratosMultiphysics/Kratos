import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
import KratosMultiphysics.vtk_output_process as vtk_output_process
import math
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestFluxBasedRedistance(KratosUnittest.TestCase):

    def _ExpectedDistance(self,x,y,z):
        d = x
        if( d > 0.2):
            d = 0.2
        if( d < -0.2):
            d = -0.2
        return x

    def test_model_part_sub_model_parts(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VOLUME)
        KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)


        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0, self._ExpectedDistance(node.X,node.Y,node.Z)  )

        linear_solver = linear_solver_factory.ConstructSolver( KratosMultiphysics.Parameters( """ { "solver_type" : "skyline_lu_factorization" } """ ) )

        model_part.CloneTimeStep(1.0)

        settings = KratosMultiphysics.Parameters()
        KratosMultiphysics.FluxBasedRedistanceProcess3D(model_part, linear_solver, settings).Execute()

        max_distance = -1.0
        min_distance = +1.0
        for node in model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        
        vtk_output_parameters = KratosMultiphysics.Parameters("""{
            "Parameters" : {
                "model_part_name"                    : "Main",
                "file_format"                        : "ascii",
                "output_precision"                   : 8,
                "output_interval"                    : 2,
                "output_sub_model_parts"             : true,
                "output_path"                        : "test_vtk_output",
                "nodal_solution_step_data_variables" : ["DISTANCE"],
                "nodal_data_value_variables"         : ["POTENTIAL_GRADIENT","ADJOINT_SCALAR_1"],
                "nodal_flags"                        : [],
                "element_data_value_variables"       : [],
                "condition_data_value_variables"     : [],
                "condition_flags"                    : []
            }
        }""")

        vtk_output_parameters["Parameters"]["model_part_name"].SetString("Main")
        vtk_output_parameters["Parameters"]["file_format"].SetString("ascii")
        #output_process = vtk_output_process.Factory(vtk_output_parameters, current_model)
        #output_process.ExecuteInitialize()
        #output_process.ExecuteInitializeSolutionStep()
        #output_process.ExecuteFinalizeSolutionStep()
        #output_process.PrintOutput()


        self.assertAlmostEqual(max_distance, 0.49554221316850783)
        self.assertAlmostEqual(min_distance, -0.4950998362886954)


if __name__ == '__main__':
    KratosUnittest.main()
