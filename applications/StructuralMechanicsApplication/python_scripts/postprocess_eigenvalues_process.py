import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return PostProcessEigenvaluesProcess(Model, settings["Parameters"])

class PostProcessEigenvaluesProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    def __init__(self, Model, settings):

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "model_part_name"   : "Structure",
                "dof_variable_name" : "DISPLACEMENT",
                "animation_steps"   :  1
            }
            """
        );

        settings.ValidateAndAssignDefaults(default_settings)

        KratosMultiphysics.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.eigenvaluevariable = getattr(StructuralMechanicsApplication,  "EIGENVALUE_VECTOR")
        self.eigenvectorvariable = getattr(StructuralMechanicsApplication, "EIGENVECTOR_MATRIX")
        self.dof_variable_name = KratosMultiphysics.KratosGlobals.GetVariable( settings["dof_variable_name"].GetString()) 
        self.animation_steps =  settings["animation_steps"].GetInt()
        
        #self.dim = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                                                                              
    def ExecuteInitialize(self):
        pass
    
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass
              
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
        singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions
        
        ## Debug
        #for node in self.model_part.Nodes:
            #EigenMatrix = node.GetValue(self.eigenvectorvariable)
            #print(EigenMatrix)
        
        eigen_values = [ev for ev in self.model_part.ProcessInfo[self.eigenvaluevariable]]
        for evs, count in zip(eigen_values, range(len(eigen_values))):
            # We create a different ouput for each eigenvalue
            output_file = "EigenValue_w="+str(round(evs, 3))
            
            gid_io = KratosMultiphysics.GidIO(output_file, gid_mode, singlefile, deformed_mesh_flag, write_conditions)
            gid_io.InitializeMesh(0.0) 
            gid_io.WriteMesh(self.model_part.GetMesh())
            gid_io.WriteNodeMesh(self.model_part.GetMesh())
            gid_io.FinalizeMesh() 
            
            gid_io.InitializeResults(0.0, self.model_part.GetMesh())
            
            for label in range(self.animation_steps):
                angle = 2.0 * math.pi * label/self.animation_steps
                # NOTE: This is slow, because I am accessing to the same information several times
                for node in self.model_part.Nodes:
                    EigenMatrix = node.GetValue(self.eigenvectorvariable)
                    # NOTE: The DoF stored include the velocity and acceleration, solve this in the future
                    dof_values = KratosMultiphysics.Vector(3)
                    dof_values[0] = math.cos(angle) * EigenMatrix[count, 0]
                    dof_values[1] = math.cos(angle) * EigenMatrix[count, 1]
                    dof_values[2] = math.cos(angle) * EigenMatrix[count, 2]
                    node.SetSolutionStepValue(self.dof_variable_name, 0, dof_values) 
                
                gid_io.WriteNodalResults(self.dof_variable_name, self.model_part.Nodes, label, 0)

            
            
