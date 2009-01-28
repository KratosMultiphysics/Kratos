#importing the Kratos Library
from Kratos import *

class ExplicitCoupling:
    
    def __init__(self,fluid_model_part,structure_model_part,structural_solver,fluid_solver,mesh_solver,mapper,domain_size):
        self.fluid_solver = fluid_solver
        self.structural_solver = structural_solver
        self.mesh_solver = mesh_solver
        self.mapper = mapper
        self.fluid_model_part = fluid_model_part
        self.structure_model_part = structure_model_part

        self.utilities = VariableUtils()

        #initialize the list of interface nodes
        self.interface_fluid_nodes = (self.utilities).SelectNodeList(IS_INTERFACE,1.0,fluid_model_part.Nodes)

    def Solve_fluid(self):
        (self.fluid_solver).Solve()


    def Solve(self):
        print "qui"
        #solve the structure (prediction)
        (self.structural_solver).Solve()
        print "Structural Prediction: OK"
        
        ##map displacements to the structure
        (self.mapper).StructureToFluid_VectorMap(DISPLACEMENT,DISPLACEMENT)
        print "Displacement Map: OK"

        ##move the mesh
        (self.mesh_solver).Solve()
        print "Mesh Movement: OK"
        
        ##set the fluid velocity at the interface to be equal to the corresponding mesh velocity
        (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);

        ##solve the fluid
        (self.fluid_solver).Solve()
        print "Fluid Solution: OK"

        ##map displacements to the structure
        (self.mapper).FluidToStructure_ScalarMap(PRESSURE,POSITIVE_FACE_PRESSURE)
        print "Pressure Map OK"
        
        #solve the structure (correction)
        (self.structural_solver).Solve()
        print "Structural Correction: OK"
        
    
