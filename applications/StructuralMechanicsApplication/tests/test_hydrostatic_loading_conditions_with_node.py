from __future__ import print_function, absolute_import, division
import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math
import numpy as np


class TestLoadingConditionsSurface(KratosUnittest.TestCase):
    def test_HydrostaticLoad3D3N(self):
        dim = 2
        mp = KratosMultiphysics.ModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        
        #create nodes
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,0.0,1.0,0.0)
        mp.CreateNewNode(3,0.0,0.0,1.0)
        mp.CreateNewNode(4,0.0,1.0,1.0)

       

        specific_weight= 10.0

        centre = [0.0,0.0,0.5]
        radius = 10.0
        normal = [0.0,0.0,1.0]
        fluid_volume = 2.0

        # Create nodal concentrated element
        mp.CreateNewNode(5,centre[0],centre[1],centre[2])
        

        
    
        
        #ensure that the property 1 is created
        self.properties = mp.GetProperties()[1]
        

        self.properties.SetValue(StructuralMechanicsApplication.FREE_SURFACE_RADIUS, radius)
        self.properties.SetValue(StructuralMechanicsApplication.FLUID_VOLUME, fluid_volume)

        self.properties.SetValue(StructuralMechanicsApplication.FREE_SURFACE_CENTRE, centre)
        self.properties.SetValue(StructuralMechanicsApplication.FREE_SURFACE_NORMAL, normal) 
        self.properties.SetValue(StructuralMechanicsApplication.SPECIFIC_WEIGHT, specific_weight) 

        mp.SetBufferSize(2)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        
            
        nodeIds1 = [1,3,2]
        nodeIds2 = [2,3,4]
        nodeIDs2 = [5]
        cond_hydro1 = mp.CreateNewCondition("HydrostaticLoadCondition3D3N", 1, nodeIds1, mp.GetProperties()[1])
        cond_hydro2 = mp.CreateNewCondition("HydrostaticLoadCondition3D3N", 2, nodeIds2, mp.GetProperties()[1])
      
        elem1 = mp.CreateNewElement("ShellThinElement3D3N", 1, nodeIds1, mp.GetProperties()[1])
        elem2 = mp.CreateNewElement("ShellThinElement3D3N", 2, nodeIds2, mp.GetProperties()[1])
        nodal_elem = mp.CreateNewElement("NodalConcentratedFluidElement3D1N", 3, nodeIDs2, mp.GetProperties()[1])

        nodal_elem.SetValue(StructuralMechanicsApplication.STIFFNESS_SCALING,2.0)
        nodal_elem.SetValue(StructuralMechanicsApplication.IS_DYNAMIC, True)

        


        for node in mp.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0,0.0)

            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 1,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 1,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 1,0.0)
            

        avg_nodes = 10
        avg_elems = 10

        

        neighbhor_finder = KratosMultiphysics.FindNodalNeighboursProcess(mp,avg_elems, avg_nodes)

        neighbhor_finder.Execute()
        
       


        lhs_cond1 = KratosMultiphysics.Matrix(0,0)
        rhs_cond1 = KratosMultiphysics.Vector(0)
        lhs_cond2 = KratosMultiphysics.Matrix(0,0)
        rhs_cond2 = KratosMultiphysics.Vector(0)

        lhs_elem = KratosMultiphysics.Matrix(0,0)
        rhs_elem = KratosMultiphysics.Vector(0)

        print(mp)

        master_creator = StructuralMechanicsApplication.BuildMasterConditionsForHydrostaticLoadingProcess(mp,nodal_elem)
        master_creator.Execute()

        
        nodal_elem.CalculateLocalSystem(lhs_elem,rhs_elem,mp.ProcessInfo)
        cond_hydro1.CalculateLocalSystem(lhs_cond1,rhs_cond1,mp.ProcessInfo)
        cond_hydro2.CalculateLocalSystem(lhs_cond2,rhs_cond2,mp.ProcessInfo)
        print("Current volume :: ",self.properties.GetValue(StructuralMechanicsApplication.CURRENT_FLUID_VOLUME))
        

        for node in mp.Nodes:
            print("Nodal Normal :: ",node.GetSolutionStepValue(KratosMultiphysics.NORMAL),0)
        #print(lhs_elem)
        #print('##########################################################')
        #print(rhs_elem)
        

        print(lhs_cond1)
        print('#################')
        print(lhs_cond2)
        print('#################')
        print(rhs_cond1)
        print('#################')
        print(rhs_cond2)

        print('########Nodal Element#########')
        print(lhs_elem)
        print('##########################################################')
        print(rhs_elem)
        force = 0.0
        for i in range(0,len(rhs_cond1)):
            force += rhs_cond1[i]
            force += rhs_cond2[i]

        print(force)



        self.assertAlmostEqual(0,0)
        
if __name__ == '__main__':
    KratosUnittest.main()
