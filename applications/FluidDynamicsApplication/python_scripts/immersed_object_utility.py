from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication 

# Other imports
import math

class ImmersedBCsUtility():
    def __init__(self, Model, settings):
        self.settings = settings;
        self.Model = Model

        default_settings = KratosMultiphysics.Parameters(
            """
            {
                "immersed_structure_model_part" : "immersed_structure_model_part",
                "fluid_domain_model_part"       : "fluid_domain_model_part"
            }
            """
        )

        self.settings.ValidateAndAssignDefaults(default_settings)

        self.structure_mp   = self.Model[self.settings["immersed_structure_model_part"].GetString()]
        self.fluid_mp       = self.Model[self.settings["fluid_domain_model_part"].GetString()]
        
        self.domain_size = self.fluid_mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
 
        if(self.domain_size == 2):
            self.fluid_space_database = KratosMultiphysics.BinBasedFastPointLocator2D(self.fluid_mp)
            self.structure_space_database = KratosMultiphysics.BinBasedFastPointLocator2D(self.structure_mp)
        else:
            self.fluid_space_database = KratosMultiphysics.BinBasedFastPointLocator3D(self.fluid_mp)
            self.structure_space_database = KratosMultiphysics.BinBasedFastPointLocator3D(self.structure_mp)
            
    
    def ApplyFluidImmersedBCs(self):
        
        self.immersed_nodes = []
        
        print(self.structure_mp)
        
        
        zero = KratosMultiphysics.Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0
        
        self.structure_space_database.UpdateSearchDatabase()
        
        N = KratosMultiphysics.Vector(self.domain_size+1)
        coords =  KratosMultiphysics.Array3()
        pelem =  KratosMultiphysics.Element(-1) #UGLY! here i create an empty pointer
        grad =  KratosMultiphysics.Vector(3)

        for node in self.fluid_mp.Nodes: #nodes on the skin of the rotor            
            #now find if inside
            coords[0] = node.X
            coords[1] = node.Y
            coords[2] = node.Z
            found = self.structure_space_database.FindPointOnMesh(coords, N, pelem, 1000, 1e-9) 
            print("-----------",node.Id,coords, found)
            if(found): #fluid node falls within the structure
                
                v = 0.0
                v[0] = 0.0
                v[1] = 0.0
                v[2] = 0.0
                k = 0
                for p in pelem.GetNodes():
                    v += N[k]*p.GetSolutionStepValue(KratosMultiphysics.VELOCITY)                    
                    k+=1
                
                
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,v)
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
                self.immersed_nodes.append(node)
        print("cccccccccc")
    
    def ApplyStructureForces(self):
        
        zero = KratosMultiphysics.Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0

        
        self.immersed_nodes = []
                
        self.fluid_space_database.UpdateSearchDatabase()
        
        KratosMultiphysics.CalculateNodalAreaProcess(self.fluid_mp, self.fluid_mp.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]).Execute()
        
        N = KratosMultiphysics.Vector(self.domain_size+1)
        coords =  KratosMultiphysics.Array3()
        pelem =  KratosMultiphysics.Element(-1) #UGLY! here i create an empty pointer
        grad =  KratosMultiphysics.Vector(3)

        for node in self.structure_mp.Nodes: #nodes on the skin of the rotor            
            #now find if inside
            coords[0] = node.X
            coords[1] = node.Y
            coords[2] = node.Z
            found = self.fluid_space_database.FindPointOnMesh(coords, N, pelem, 1000, 1e-9) 
            
            if(found): #structure node falls within the fluid
                t = KratosMultiphysics.Vector(3)
                t[0] = 0.0
                t[1] = 0.0
                t[2] = 0.0
                k = 0
                for p in pelem.GetNodes():
                    A = p.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
                    rho = p.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                    t += (N[k]/(A*rho))*p.GetSolutionStepValue(KratosMultiphysics.REACTION)                    
                    k+=1
                
                node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,t)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,zero)


    def Unfix(self):
        for node in self.immersed_nodes:
            node.Free(KratosMultiphysics.VELOCITY_X)
            node.Free(KratosMultiphysics.VELOCITY_Y)
            node.Free(KratosMultiphysics.VELOCITY_Z)