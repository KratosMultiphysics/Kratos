import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EmbeddedElementDeactivation(Model, settings["Parameters"])

class EmbeddedElementDeactivation(KratosMultiphysics.Process):
    
    def __init__(self, Model, settings):
        
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "mesh_id"                   : 0,
            "model_part_name"           : "CHOOSE_FLUID_MODELPART_NAME",
            "check_at_each_time_step"   : false
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.nnodes = self.fluid_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE) + 1 # Note that this process only works for simplicial elements (tetrahedras and triangles)
        self.check_at_each_time_step = settings["check_at_each_time_step"].GetBool()
                
                
    def ExecuteInitialize(self):
        
        # Obtain the minimum nodal element size (the minimum element size of all the elements surrounding a node)
        # Itialize the nodal area/volume
        for elem in self.fluid_model_part.Elements:
            for node in elem.GetNodes():
                node.SetValue(KratosFluid.NODAL_H, elem.GetArea())  
                
        # Ensure that the minimum value of the surrounding elements is kept
        for elem in self.fluid_model_part.Elements:
            for node in elem.GetNodes():
                node.SetValue(KratosFluid.NODAL_H, min(node.GetValue(KratosFluid.NODAL_H), elem.GetArea()))  
                
        # Avoid the "bad" intersections
        for node in self.fluid_model_part.Nodes:
            
            h = node.GetValue(KratosFluid.NODAL_H)
            if (self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                d_tol = 1e-2*(h**(1/2))
            else:
                d_tol = 1e-2*(h**(1/3))
            
            if((node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0) and (node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < d_tol)):
                print("Modifying distance value! d_tol = ", d_tol)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -1e-5)

        # Deactivate the full negative distance elements
        for elem in self.fluid_model_part.Elements:
            inside = 0
            
            # Check the sign of the distance function at the element nodes
            for node in elem.GetNodes():
                if(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0):
                    inside += 1
                    
            # If all the nodes have negative distance value (structure domain) deactivate the element
            if(inside == self.nnodes):
                elem.Set(KratosMultiphysics.ACTIVE,False)
                
                
    def ExecuteInitializeSolutionStep(self):
        
        if (self.check_at_each_time_step == True):            
            
            self.mod_d_ids = []
            self.mod_d_values = []
            
            # Obtain the minimum nodal element size (the minimum element size of all the elements surrounding a node)
            for elem in self.fluid_model_part.Elements:
                for node in elem.GetNodes():
                    node.SetValue(KratosFluid.NODAL_H, min(node.GetValue(KratosFluid.NODAL_H), elem.GetArea()))
                    
            # Avoid the "bad" intersections
            for node in self.fluid_model_part.Nodes:
                
                h = node.GetValue(KratosFluid.NODAL_H)
                if (self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                    d_tol = 5e-2*(h**(1/2))
                else:
                    d_tol = 5e-2*(h**(1/3))
                
                if((node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0) and (node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < d_tol)):
                    self.mod_d_ids.append(node.Id)
                    self.mod_d_values.append(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE))
                    
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -1e-5)
            
            # Deactivate the full negative distance elements
            for elem in self.fluid_model_part.Elements:
                inside = 0
                
                # Check the sign of the distance function at the element nodes
                for node in elem.GetNodes():
                    if(node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0):
                        inside += 1
                        
                # If all the nodes have negative distance value (structure domain) deactivate the element
                if(inside == self.nnodes):
                    elem.Set(KratosMultiphysics.ACTIVE, False)
                # Otherwise, activate the element again
                elif(inside < self.nnodes):
                    elem.Set(KratosMultiphysics.ACTIVE, True)
                # Security check
                else:
                    raise Exception ("ERROR: The number of inside elements is larger than the element nodes!")
        
        else:
            pass
            
        
        def ExecuteFinalizeSolutionStep(self):
            
            # Recover the original distances in case it has been modified
            if (self.check_at_each_time_step == True):
                aux_counter = 0
                
                for node in self.fluid_model_part.Nodes:
                    if (node.Id() == self.mod_d_ids[aux_counter]):
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, self.mod_d_values[aux_counter])
                        aux_counter += 1
            else:
                pass
