from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import math
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

class RVEAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters,
            boundary_mp_name = "Structure.DISPLACEMENT_Displacement_Auto1",
            averaging_mp_name = "Structure.computing_domain"
            ):
        
        #######################################
        #######################################
        #######################################
        ##parameters to be adjusted by the user
        self.averaging_volume = -1.0 #it will be computed in initialize     
        self.boundary_mp_name = boundary_mp_name
        self.averaging_mp_name = averaging_mp_name
        #######################################
        #######################################
        #######################################
                
        domain_size = project_parameters["solver_settings"]["domain_size"].GetInt()
        if(domain_size == 2):
            self.strain_size = 3
        else:
            self.strain_size = 6

        #pseudo time to be used for output
        self.time = 0.0
            
        super(RVEAnalysis, self).__init__(model, project_parameters)


    ##here populate the submodelparts to be used for periodicity
    def ModifyInitialGeometry(self):
        super(RVEAnalysis, self).ModifyInitialGeometry()

        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]  

        #construct auxiliary modelparts
        self.min_corner,self.max_corner = self._DetectBoundingBox(averaging_mp)
        self._ConstructFaceModelParts(self.min_corner, self.max_corner, boundary_mp)

        self.averaging_volume = (self.max_corner[0]-self.min_corner[0]) * (self.max_corner[1]-self.min_corner[1]) * (self.max_corner[2]-self.min_corner[2])

    def InitializeSolutionStep(self):
        raise Exception("should use the _CustomInitializeSolutionStep instead of this")


        

    def __CustomInitializeSolutionStep(self,strain,boundary_mp, averaging_mp):
        self.ApplyBoundaryConditions() #here the processes are called

        #construct MPCs according to the provided strain
        self._ApplyPeriodicity(strain,averaging_mp, boundary_mp) 
        
        ##apply BCs for RVE according to the provided strain
        self._ApplyMinimalConstraints(averaging_mp,strain,self.min_corner,self.max_corner)

        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def RunSolutionLoop(self):

        #super(RVEAnalysis, self).RunSolutionLoop()
        perturbation = 1e-6
        
        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]  

        stress_and_strain = []
        if(self.strain_size == 3): ##2D case - ordering s00 s11 s01
            stress_and_strain.append( self._ComputeEquivalentStress(0,0,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(1,1,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(0,1,perturbation, boundary_mp, averaging_mp) )
        elif(self.strain_size == 6): ##3D case - ordering:  s00 s11 s22 s01 s12 s02
            stress_and_strain.append( self._ComputeEquivalentStress(0,0,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(1,1,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(2,2,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(0,1,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(1,2,perturbation, boundary_mp, averaging_mp) )
            stress_and_strain.append( self._ComputeEquivalentStress(0,2,perturbation, boundary_mp, averaging_mp) )

        C = self._ComputeEquivalentElasticTensor(stress_and_strain, perturbation)
        averaging_mp.SetValue(KratosMultiphysics.StructuralMechanicsApplication.EIGENVECTOR_MATRIX, C)
        self._MatrixOutput(C)


    def _DetectBoundingBox(self,mp):
        min_corner = KratosMultiphysics.Array3()
        min_corner[0] = 0.0
        min_corner[1] = 0.0
        min_corner[2] = 0.0

        max_corner = KratosMultiphysics.Array3()
        max_corner[0] = 0.0
        max_corner[1] = 0.0
        max_corner[2] = 0.0

        for node in mp.Nodes:
            x = node.X 
            min_corner[0] = min(min_corner[0], x)
            max_corner[0] = max(max_corner[0], x)

            y = node.Y 
            min_corner[1] = min(min_corner[1], y)
            max_corner[1] = max(max_corner[1], y)

            z = node.Z 
            min_corner[2] = min(min_corner[2], z)
            max_corner[2] = max(max_corner[2], z)

        print("Boundng box detected")
        print("min corner = ",min_corner)
        print("max corner = ",max_corner)

        return min_corner, max_corner

    def __PopulateMp(self,face_name, coordinate, component, eps, mp):
        if(len(mp.Conditions) == 0):
            raise Exception("boundary_mp is expected to have conditions and has none")

        mp = mp.GetRootModelPart()

        if not mp.HasSubModelPart(face_name):
            mp.CreateSubModelPart(face_name)
        face_mp = mp.GetSubModelPart(face_name)

        for cond in mp.Conditions:
            xc = cond.GetGeometry().Center() 
            if abs(xc[component]-coordinate) < eps:
                face_mp.AddCondition(cond)

        node_ids = set()
        for cond in face_mp.Conditions:
            for node in cond.GetNodes():
                if(not node.Is(KratosMultiphysics.SLAVE)):
                    node_ids.add(node.Id)
                    node.Set(KratosMultiphysics.SLAVE)



        face_mp.AddNodes(list(node_ids))
        return face_mp

    def _ConstructFaceModelParts(self,min_corner, max_corner, mp):
        
        eps = 0.0001*(max_corner[0] - min_corner[0])/len(mp.Nodes)

        for node in mp.Nodes:
            node.Set(KratosMultiphysics.SLAVE,False)
            node.Set(KratosMultiphysics.MASTER,False)

        #populate the slave faces
        print("----------------- populate slave faces")
        self.max_x_face = self.__PopulateMp("max_x_face", max_corner[0], 0, eps, mp)
        self.max_y_face = self.__PopulateMp("max_y_face", max_corner[1], 1, eps, mp)
        self.max_z_face = self.__PopulateMp("max_z_face", max_corner[2], 2, eps, mp)  

        #first populate the master faces (min)
        print("----------------- populate master faces")
        self.min_x_face = self.__PopulateMp("min_x_face", min_corner[0], 0, eps, mp)
        self.min_y_face = self.__PopulateMp("min_y_face", min_corner[1], 1, eps, mp)
        self.min_z_face = self.__PopulateMp("min_z_face", min_corner[2], 2, eps, mp)   

        if len(self.min_x_face.Conditions) == 0:
            raise Exception("min_x_face has 0 conditions")
        if len(self.min_y_face.Conditions) == 0:
            raise Exception("min_y_face has 0 conditions")
        if len(self.min_z_face.Conditions) == 0:
            raise Exception("min_z_face has 0 conditions")

    def _SelectClosestNode(self,mp,coords):
        min_distance = 1e30
        selected_node = 0
        for node in mp.Nodes:
            dx = node.X0 - coords[0]
            dy = node.Y0 - coords[1]
            dz = node.Z0 - coords[2]
            d = dx**2 + dy**2 + dz**2

            if(d < min_distance):
                selected_node = node
                min_distance = d

        return selected_node

    #prescribed conditions to avoid rigid body motions
    def _ApplyMinimalConstraints(self,mp, strain,min_corner,max_corner):
        aux = KratosMultiphysics.Array3()

        #point coinciding with the min_corner
        node = self._SelectClosestNode(mp, min_corner)
        node.Fix(KratosMultiphysics.DISPLACEMENT_X)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
        node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        coords_min_corner = KratosMultiphysics.Array3(node)
        coords_min_corner[0] = node.X0
        coords_min_corner[1] = node.Y0
        coords_min_corner[2] = node.Z0

        disp_min_corner = strain*coords_min_corner 
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,disp_min_corner)

    def _ComputeEquivalentElasticTensor(self,stress_and_strain,perturbation):
        C = KratosMultiphysics.Matrix(self.strain_size, self.strain_size)

        for j in range(len(stress_and_strain)):
            stress = stress_and_strain[j][0]
            for i in range(self.strain_size):
                C[i,j] = stress[i]

        inverse_perturbation = KratosMultiphysics.Matrix(self.strain_size, self.strain_size)
        inverse_perturbation.fill(0.0)
        if(self.strain_size == 3):
            inverse_perturbation[0,0] = 1.0/perturbation
            inverse_perturbation[1,1] = 1.0/perturbation
            inverse_perturbation[2,2] = 0.5/perturbation        
        else:
            inverse_perturbation[0,0] = 1.0/perturbation
            inverse_perturbation[1,1] = 1.0/perturbation
            inverse_perturbation[2,2] = 1.0/perturbation
            inverse_perturbation[3,3] = 0.5/perturbation
            inverse_perturbation[4,4] = 0.5/perturbation
            inverse_perturbation[5,5] = 0.5/perturbation

        C = C*inverse_perturbation

        return C
        
    def _MatrixOutput(self,C,filename="rve_elasticity_tensor.txt"):
        f = open(filename,'w')

        if(self.strain_size == 3): #2D
            f.write( str(C[0,0])+ " " + str(C[0,1])+ " " + str(C[0,2]) + "\n" )
            f.write( str(C[1,0])+ " " + str(C[1,1])+ " " + str(C[1,2]) + "\n" )
            f.write( str(C[2,0])+ " " + str(C[2,1])+ " " + str(C[2,2]) + "\n" )

        elif(self.strain_size == 6):
            for i in range(6):
                f.write( str(C[i,0])+ " " + str(C[i,1])+ " " + str(C[i,2])+ " " + str(C[i,3])+ " " + str(C[i,4])+ " " + str(C[i,5])+ "\n" )

        # if(self.strain_size == 3): #2D
        #     f.write( str(C[0,0])+ " " + str(C[0,1])+ " " + str(C[0,2]) + "\n" )
        #     f.write( str(C[1,0])+ " " + str(C[1,1])+ " " + str(C[1,2]) + "\n" )
        #     f.write( str(C[2,0])+ " " + str(C[2,1])+ " " + str(C[2,2]) + "\n" )

        # elif(self.strain_size == 6):
        #     f.write( str(C[0,0])+ " " + str(C[0,1])+ " " + str(C[0,3])+ " " + str(C[0,2])+ " " + str(C[0,5])+ " " + str(C[0,4])+ "\n" )
        #     f.write( str(C[1,0])+ " " + str(C[1,1])+ " " + str(C[1,3])+ " " + str(C[1,2])+ " " + str(C[1,5])+ " " + str(C[1,4])+ "\n" )
        #     f.write( str(C[3,0])+ " " + str(C[3,1])+ " " + str(C[3,3])+ " " + str(C[3,2])+ " " + str(C[3,5])+ " " + str(C[3,4])+ "\n" )
        #     f.write( str(C[2,0])+ " " + str(C[2,1])+ " " + str(C[2,3])+ " " + str(C[2,2])+ " " + str(C[2,5])+ " " + str(C[2,4])+ "\n" )
        #     f.write( str(C[5,0])+ " " + str(C[5,1])+ " " + str(C[5,3])+ " " + str(C[5,2])+ " " + str(C[5,5])+ " " + str(C[5,4])+ "\n" )
        #     f.write( str(C[4,0])+ " " + str(C[4,1])+ " " + str(C[4,3])+ " " + str(C[4,2])+ " " + str(C[4,5])+ " " + str(C[4,4])+ "\n" )
        
        f.close()
                
        


    def _ComputeEquivalentStress(self,i,j,perturbation, boundary_mp, averaging_mp):

        #here use a pseudotime for output
        self.time = self.time + 1.0
        averaging_mp.GetRootModelPart().CloneTimeStep(self.time)

        strain = KratosMultiphysics.Matrix(3,3)
        strain.fill(0.0)

        strain[i,j] = perturbation
        strain[j,i] = perturbation

        strain_vector = KratosMultiphysics.Vector(self.strain_size)
        if(self.strain_size == 2):
            strain_vector[0] = strain[0,0]
            strain_vector[1] = strain[1,1]
            strain_vector[2] = 2.0*strain[1,2]
        elif(self.strain_size == 6):
            strain_vector[0] = strain[0,0]
            strain_vector[1] = strain[1,1]
            strain_vector[2] = strain[2,2]
            strain_vector[3] = 2.0*strain[0,1]
            strain_vector[4] = 2.0*strain[1,2]
            strain_vector[5] = 2.0*strain[0,2]

        #self._GetSolver().Initialize()
        self.__CustomInitializeSolutionStep(strain,boundary_mp, averaging_mp)

        self._GetSolver().Predict()
        print(self._GetSolver().GetComputingModelPart())

        self._GetSolver().SolveSolutionStep()
        process_info = averaging_mp.ProcessInfo
        avg_stress = KratosMultiphysics.Vector(self.strain_size)
        avg_stress.fill(0.0)
        measured_volume = 0.0

        for elem in averaging_mp.Elements:
            tmp = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PK2_STRESS_VECTOR,process_info)
            ngauss = len(tmp)
            A = elem.GetGeometry().Area() 
            measured_volume += A
            Agauss = A/ngauss #TODO: this is only valid for gauss points with the same weight. should be generalized
            for item in tmp:
                avg_stress = avg_stress + item*Agauss

        self._GetSolver().Clear()            
        
        print("measured volume = ", measured_volume)
     
                
        avg_stress /= self.averaging_volume

        print("avg_stress = ", avg_stress)
       

            

        self.OutputSolutionStep()

        #reset position of nodes
        for node in averaging_mp.Nodes:
            node.X = node.X0
            node.Y = node.Y0
            node.Z = node.Z0
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0,0.0)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0,0.0)



        return avg_stress, strain_vector




    def _ApplyPeriodicity(self,strain, volume_mp, boundary_mp):

        #clear 
        for constraint in volume_mp.GetRootModelPart().MasterSlaveConstraints:
            constraint.Set(KratosMultiphysics.TO_ERASE)
        volume_mp.GetRootModelPart().RemoveMasterSlaveConstraintsFromAllLevels(KratosMultiphysics.TO_ERASE)

        dx = self.max_corner[0] - self.min_corner[0]
        dy = self.max_corner[1] - self.min_corner[1]
        dz = self.max_corner[2] - self.min_corner[2]

        periodicity_utility = KratosMultiphysics.StructuralMechanicsApplication.RVEPeriodicityUtility(self._GetSolver().GetComputingModelPart())

        ##assign periodicity to faces
        periodicity_utility.AssignPeriodicity(self.min_x_face,self.max_x_face,strain, KratosMultiphysics.Vector([dx, 0.0, 0.0]))
        periodicity_utility.AssignPeriodicity(self.min_y_face,self.max_y_face,strain, KratosMultiphysics.Vector([0.0, dy, 0.0]))
        periodicity_utility.AssignPeriodicity(self.min_z_face,self.max_z_face,strain, KratosMultiphysics.Vector([0.0, 0.0, dz]))

        periodicity_utility.Finalize(KratosMultiphysics.DISPLACEMENT)

        #start from the exact solution in the case of a constant strain
        x = KratosMultiphysics.Vector(3)
        for node in volume_mp.Nodes:
            x[0] = node.X0
            x[1] = node.Y0
            x[2] = node.Z0
            d = strain*x 
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,d)



