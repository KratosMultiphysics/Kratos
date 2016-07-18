import KratosMultiphysics
import python_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
        
    return ApplyHydrostaticOutletProcess(Model, settings["Parameters"])


class ApplyHydrostaticOutletProcess(python_process.PythonProcess):
    def __init__(self, Model, Parameters ):
        
        self.model_part = Model[Parameters["model_part_name"].GetString()]
        self.parameters = Parameters      
        self.outlet_type = self.parameters["outlet_type"].GetString()
        
        if self.outlet_type == "Constant":
            
            constant_outlet_settings = KratosMultiphysics.Parameters("{}")
            constant_outlet_settings.AddValue("model_part_name",self.parameters["model_part_name"]) 
            constant_outlet_settings.AddValue("mesh_id",self.parameters["mesh_id"]) 
            constant_outlet_settings.AddEmptyValue("is_fixed").SetBool(True)
            constant_outlet_settings.AddValue("value",self.parameters["value"]) 
            constant_outlet_settings.AddValue("variable_name",self.parameters["variable_name"]) 
            
            KratosMultiphysics.Process.__init__(self)
            self.constant_outlet_process = KratosMultiphysics.ApplyConstantScalarValueProcess(self.model_part, constant_outlet_settings)
        
    def ExecuteInitialize(self):          
              
        if self.outlet_type == "Constant":
            self.constant_outlet_process.ExecuteInitialize()
 
                
    def ExecuteBeforeSolutionLoop(self):
        
        # The "hydrostatic" case is set in the ExecuteBeforeSolutionLoop() step to ensure that the body force has been set in all the nodes in the ExecuteInitialize() step.
        if self.outlet_type == "Hydrostatic":
            tol = 1e-5
            body_force = list(self.model_part.Nodes)[0].GetSolutionStepValue(KratosMultiphysics.BODY_FORCE,0)
            #~ body_force_norm = (body_force[0]*body_force[0]+body_force[1]*body_force[1]+body_force[2]*body_force[2])**0.5
            
            h_top = self.parameters["h_top"].GetDouble()
            rho = self.parameters["hyd_rho"].GetDouble()
        
            ##### Body force in x-direction #####
            if abs(body_force[0]) > tol:
                x_coords_list = []
                gravity = abs(body_force[0])
                
                for node in self.model_part.Nodes:
                    x_coords_list.append(node.X)
            
                # Body force in negative direction
                if body_force[0]<0.0:
                    self.ref_coord = max(x_coords_list)
                    
                    for node in self.model_part.Nodes:
                        hyd_pres = rho*gravity*(h_top+(self.ref_coord-node.X))
                        cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0)
                
                        node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0,hyd_pres+cur_pres)
                        #~ node.Fix(KratosMultiphysics.EXTERNAL_PRESSURE)
            
                # Body force in positive directon
                elif body_force[0]>0.0:
                    self.ref_coord = min(x_coords_list)
                    
                    for node in self.model_part.Nodes:
                        hyd_pres = rho*gravity*(h_top+(node.X-self.ref_coord))
                        cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0)
                
                        node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0,hyd_pres+cur_pres)
                        #~ node.Fix(KratosMultiphysics.EXTERNAL_PRESSURE)
            
            ##### Body force in y-direction #####
            if abs(body_force[1])>1e-5:
                y_coords_list = []
                gravity = abs(body_force[1])
            
                for node in self.model_part.Nodes:
                    y_coords_list.append(node.Y)
            
                # Body force in negative direction
                if body_force[1]<0.0:
                    self.ref_coord = max(y_coords_list)
                    
                    for node in self.model_part.Nodes:
                        hyd_pres = rho*gravity*(h_top+(self.ref_coord-node.Y))
                        cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0)
                
                        node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0,hyd_pres+cur_pres)
                        #~ node.Fix(KratosMultiphysics.EXTERNAL_PRESSURE)
            
                # Body force in positive direction
                elif body_force[1]>0.0:
                    self.ref_coord = min(y_coords_list)
                    
                    for node in self.model_part.Nodes:
                        hyd_pres = rho*gravity*(h_top+(node.Y-self.ref_coord))
                        cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0)
                
                        node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0,hyd_pres+cur_pres)
                        #~ node.Fix(KratosMultiphysics.EXTERNAL_PRESSURE)
                
            ##### Body force in z-direction #####
            if abs(body_force[2])>1e-5:
                z_coords_list = []
                gravity = abs(body_force[2])
            
                for node in self.model_part.Nodes:
                    z_coords_list.append(node.Z)
                            
                # Body force in negative direction
                if body_force[2]<0.0:
                    self.ref_coord = max(z_coords_list)
                    
                    for node in self.model_part.Nodes:
                        hyd_pres = rho*gravity*(h_top+(self.ref_coord-node.Z))
                        cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0)
                
                        node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0,hyd_pres+cur_pres)
                        #~ node.Fix(KratosMultiphysics.EXTERNAL_PRESSURE)
            
                # Body force in positive direction
                elif body_force[2]>0.0:
                    self.ref_coord = min(z_coords_list)
                    
                    for node in self.model_part.Nodes:
                        hyd_pres = rho*gravity*(h_top+(node.Z-self.ref_coord))
                        cur_pres = node.GetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0)
                
                        node.SetSolutionStepValue(KratosMultiphysics.EXTERNAL_PRESSURE,0,hyd_pres+cur_pres)
                        #~ node.Fix(KratosMultiphysics.EXTERNAL_PRESSURE)

