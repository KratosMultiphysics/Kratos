import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
def GetFilePath(fileName):
    file_path = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), 
            "../../../kratos/tests/auxiliar_files_for_python_unittest/mdpa_files", 
            fileName
        )
    )
    return file_path 


class NavierStokesFractionalVectorialConvectionTest(KratosUnittest.TestCase):

    def tearDown(self):
        # Remove the .time file
        try:
            os.remove('levelset_convection_process_mesh.time')
        except :
            pass
    def _straight_velocity_field(self,x,y,z):
      return[0.1, 0, 0]
    def _gradient_velocity_field(self,x,y,z):
      return[50-x, 0, 0]
    def _acceleration_velocity_field_n(self,x,y,z):
      return[x**2-0.1, 0, 0]
    def _acceleration_velocity_field_n_1(self,x,y,z):
      return[x**2-0.2, 0, 0]
    def _acceleration_velocity_field_n_2(self,x,y,z):
      return[x**2-0.3, 0, 0]
    def _set_and_fill_buffer(self,mp,buffer_size,delta_time):
        # Set buffer size
        mp.SetBufferSize(buffer_size)

        # Fill buffer
        time = mp.ProcessInfo[KratosMultiphysics.TIME]
        time = time - delta_time * (buffer_size)
        mp.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        mp.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, delta_time)
        for size in range(0, buffer_size):
            step = size - (buffer_size -1)
            mp.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            #delta_time is computed from previous time in process_info
            mp.CloneTimeStep(time)

    def test_navier_stokes_fractional_convection_straight_field_2D(self):
        # create model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        # Add convection problem required historical variables
        model_part.AddNodalSolutionStepVariable(KratosCFD.FRACTIONAL_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        # Define the mesh. 
        # NOTE:To avoid having another MDPA file in the tests, # the MDPA used in the classic convection problem in the Kratos core is being used.
        KratosMultiphysics.ModelPartIO(GetFilePath("levelset_convection_process_mesh")).ReadModelPart(model_part)

        # Define the delta time and clone accordingly based on the buffer size.
        delta_time =0.1
        buffer_size = 3
        self._set_and_fill_buffer(model_part,buffer_size,delta_time)

        # Add corresponding DoFs and set the domain. 
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_Y, model_part)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
   
  
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY,0, self._straight_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, self._straight_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, self._straight_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 2, self._straight_velocity_field(node.X, node.Y, node.Z))

        # Fix the inlet boundary condition since it is a hyperbolic problem.
        for node in model_part.Nodes:
            if node.X<0.001:
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_X)
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_Y)
            if node.Y>0.99 and node.Y<0.001:
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_Y)

        # Compute the BDF
        KratosMultiphysics.TimeDiscretization.BDF(2).ComputeAndSaveBDFCoefficients(model_part.ProcessInfo)

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))

        KratosMultiphysics.FindGlobalNodalNeighboursProcess(model_part).Execute()
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "element_type": "ns_fractional_velocity_convection"
            }""")
        KratosCFD.TwoFLuidNavierStokesVectorialFractionalConvectionProcess2D(model_part,
                linear_solver,
                levelset_convection_settings).Execute()
        
        # Expected Results: The expected outcome is a velocity field that remains constant throughout the domain,
        # matching the initial condition. 
        for node in model_part.Nodes:
            if node.Id == 54:
               velocity  = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
               velocity_fractional  = node.GetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY_X,0)
               ref_value = velocity_fractional-velocity

        self.assertAlmostEqual(ref_value, 1e-10)
    
    def test_navier_stokes_fractional_convection_gradient_field_2D(self):
        # create model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        # Add convection problem required historical variables
        model_part.AddNodalSolutionStepVariable(KratosCFD.FRACTIONAL_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        # Define the mesh. 
        # NOTE:To avoid having another MDPA file in the tests, # the MDPA used in the classic convection problem in the Kratos core is being used.
        KratosMultiphysics.ModelPartIO(GetFilePath("levelset_convection_process_mesh")).ReadModelPart(model_part)

        # Define the delta time and clone accordingly based on the buffer size.
        delta_time =0.1
        buffer_size = 3
        self._set_and_fill_buffer(model_part,buffer_size,delta_time)

        # Add corresponding DoFs and set the domain. 
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_Y, model_part)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
   
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY,0, self._gradient_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0,  self._gradient_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1,  self._gradient_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 2,  self._gradient_velocity_field(node.X, node.Y, node.Z))

        # Fix the inlet boundary condition since it is a hyperbolic problem.
        for node in model_part.Nodes:
            if node.X<0.001:
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_X)
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_Y)
            if node.Y>0.99 and node.Y<0.001:
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_Y)

        # Compute the BDF
        KratosMultiphysics.TimeDiscretization.BDF(2).ComputeAndSaveBDFCoefficients(model_part.ProcessInfo)

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(model_part).Execute()

        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "element_type": "ns_fractional_velocity_convection"
            }""")
    
        KratosCFD.TwoFluidNavierStokesVectorialFractionalConvectionProcess2D(model_part,
                linear_solver,
                levelset_convection_settings).Execute()
        
        # Expected Results: The expected outcome is the same velocity field as the initial condition.  
        # Despite being a field with a non-zero gradient, and therefore the solution should change,  
        # it converges in one iteration because the convection of the vector field equals  
        # the past acceleration, which, given the data used, matches the convection of the velocity field at \(n\),  
        # resulting in a zero residual.

        for node in model_part.Nodes:
            if node.Id == 54:
               velocity  = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
               velocity_fractional  = node.GetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY_X,0)
               ref_value = velocity_fractional-velocity

        self.assertAlmostEqual(ref_value, 1e-10)

    def test_navier_stokes_fractional_convection_acceleration_field_2D(self):

        # create model part
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")

        # Add convection problem required historical variables
        model_part.AddNodalSolutionStepVariable(KratosCFD.FRACTIONAL_VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        # Define the mesh. 
        # NOTE:To avoid having another MDPA file in the tests, # the MDPA used in the classic convection problem in the Kratos core is being used.
        KratosMultiphysics.ModelPartIO(GetFilePath("levelset_convection_process_mesh")).ReadModelPart(model_part)
        
        # Define the delta time and clone accordingly based on the buffer size.
        delta_time =0.1
        self._set_and_fill_buffer(model_part,3,delta_time)

        # Add corresponding DoFs and set the domain. 
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.FRACTIONAL_VELOCITY_Y, model_part)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,2)

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY, 0, self._gradient_velocity_field(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, self._acceleration_velocity_field_n(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, self._acceleration_velocity_field_n_1(node.X, node.Y, node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 2, self._acceleration_velocity_field_n_2(node.X, node.Y, node.Z))
        
        # Fix the inlet boundary condition since it is a hyperbolic problem.
        for node in model_part.Nodes:
            if node.X<0.001:
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_X)
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_Y)
            if node.Y>0.99 and node.Y<0.001:
                node.Fix(KratosCFD.FRACTIONAL_VELOCITY_Y)
        
        # Compute the BDF
        KratosMultiphysics.TimeDiscretization.BDF(2).ComputeAndSaveBDFCoefficients(model_part.ProcessInfo)

        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            KratosMultiphysics.Parameters("""{"solver_type" : "skyline_lu_factorization"}"""))
        
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(model_part).Execute()
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "element_type": "ns_fractional_velocity_convection"
            }""")
        for node in model_part.Nodes:
            if node.Id == 54:
               velocity  = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X,0)
               velocity_fractional  = node.GetSolutionStepValue(KratosCFD.FRACTIONAL_VELOCITY_X,0)
               ref_value = velocity_fractional-velocity
        

        KratosCFD.TwoFluidNavierStokesVectorialFractionalConvectionProcess2D(model_part,
                linear_solver,
                levelset_convection_settings).Execute()
        self.assertAlmostEqual(ref_value, -651.9)

if __name__ == '__main__':
    KratosUnittest.main()