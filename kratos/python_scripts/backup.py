# print("/n ::TESTING:: START Calculate Jacobian /n")
            # for element in self._GetSolver().GetComputingModelPart().Elements:
            #     J = element.GetGeometry().Jacobian(0)
            #     Jred = KratosMultiphysics.Matrix([[J[0,0], J[0,1]],[J[1,0], J[1,1]]])
            #     print("J: ", element.Id, J)
            #     #print("Jred: ", element.Id, Jred)
            #     # element.GetNodes()[0].X -= 0.3xx
            #     # element.GetNodes()[0].Y -= 0.1
            #     element.GetNodes()[0].Z = 10
            #     element.GetNodes()[1].Z = 10
            #     element.GetNodes()[2].Z = 10
            #     # for node in element.GetNodes():
            #     #     # node.X = node.X0
            #     #     # node.Y = node.Y0
            #     #     node.Z = node.Z0
            #     J0 = element.GetGeometry().Jacobian(0)
            #     J0red = KratosMultiphysics.Matrix([[J0[0,0], J0[0,1]],[J0[1,0], J0[1,1]]])
            #     J0red_inv = self.InvertMatrix(J0red)
            #     print("J0: ", element.Id, J0)
            #     #print("J0red: ", element.Id, J0red)
            #     #print("J0red_inv: ", element.Id, J0red_inv)
            #     F = Jred * J0red_inv
            #     # print("F: ", element.Id, F)
            # print("/n ::TESTING:: FINISH Calculate Jacobian /n")
            # hehehehe

            # DN1DXI = -1.0
            # DN2DXI = 1.0
            # DN3DXI = 0.0
            
            # DNXI = KratosMultiphysics.Matrix(3,9)
            # # DN1X1
            # DNXI[0,0] = DN1DXI
            # DNXI[0,1] = 0.0
            # DNXI[0,2] = 0.0
            # DNXI[1,0] = 0.0
            # DNXI[1,1] = DN1DXI
            # DNXI[1,2] = 0.0
            # DNXI[2,0] = 0.0
            # DNXI[2,1] = 0.0
            # DNXI[2,2] = DN1DXI
            # # DN2XI
            # DNXI[0,3] = DN2DXI
            # DNXI[0,4] = 0.0
            # DNXI[0,5] = 0.0
            # DNXI[1,3] = 0.0
            # DNXI[1,4] = DN2DXI
            # DNXI[1,5] = 0.0
            # DNXI[2,3] = 0.0
            # DNXI[2,4] = 0.0
            # DNXI[2,5] = DN2DXI
            # # DN3XI
            # DNXI[0,6] = DN3DXI
            # DNXI[0,7] = 0.0
            # DNXI[0,8] = 0.0
            # DNXI[1,6] = 0.0
            # DNXI[1,7] = DN3DXI
            # DNXI[1,8] = 0.0
            # DNXI[2,6] = 0.0
            # DNXI[2,7] = 0.0
            # DNXI[2,8] = DN3DXI

            # # DNXI[0] = [DN1DXI,      0,       0,   DN2DXI,      0,       0,   DN3DXI,   0,     0]
            # # DNXI[1] = [     0,  DN1DXI,      0,        0,  DN2DXI,      0,       0,   DN3DXI, 0]
            # # DNXI[2] = [     0,       0, DN1DXI,        0,       0, DN2DXI,       0,   0, DN3DXI]

            # DN1DETA = -1.0
            # DN2DETA = 0.0
            # DN3DETA = 1.0
            
            # DNETA = KratosMultiphysics.Matrix(3,9)
            # # DN1X1
            # DNETA[0,0] = DN1DETA
            # DNETA[0,1] = 0.0
            # DNETA[0,2] = 0.0
            # DNETA[1,0] = 0.0
            # DNETA[1,1] = DN1DETA
            # DNETA[1,2] = 0.0
            # DNETA[2,0] = 0.0
            # DNETA[2,1] = 0.0
            # DNETA[2,2] = DN1DETA
            # # DN2XI
            # DNETA[0,3] = DN2DETA
            # DNETA[0,4] = 0.0
            # DNETA[0,5] = 0.0
            # DNETA[1,3] = 0.0
            # DNETA[1,4] = DN2DETA
            # DNETA[1,5] = 0.0
            # DNETA[2,3] = 0.0
            # DNETA[2,4] = 0.0
            # DNETA[2,5] = DN2DETA
            # # DN3XI
            # DNETA[0,6] = DN3DETA
            # DNETA[0,7] = 0.0
            # DNETA[0,8] = 0.0
            # DNETA[1,6] = 0.0
            # DNETA[1,7] = DN3DETA
            # DNETA[1,8] = 0.0
            # DNETA[2,6] = 0.0
            # DNETA[2,7] = 0.0
            # DNETA[2,8] = DN3DETA
            
            # # DNETA[0] = [DN1DETA,      0,       0,   DN2DETA,      0,       0,   DN3DETA,   0,     0]
            # # DNETA[1] = [     0,  DN1DETA,      0,        0,  DN2DETA,      0,       0,   DN3DETA, 0]
            # # DNETA[2] = [     0,       0, DN1DETA,        0,       0, DN2DETA,       0,   0, DN3DETA]

            # print("DNXI: ", DNXI)
            # print("DNETA: ", DNETA)

            # for element in self._GetSolver().GetComputingModelPart().Elements:
            #     # Jacob = KratosMultiphysics.Matrix([DNXI * element.GetNodes().GetCoordinates()], [DNETA * element.GetNodes().GetCoordinates()])
            #     # for node in element.GetNodes():
            #     #     node.X = node.X0
            #     #     node.Y = node.Y0
            #     #     node.Z = node.Z0
            #     # element.GetNodes()[0].X += 0.1
            #     # element.GetNodes()[0].Y += 0.3
            #     # element.GetNodes()[0].Z += 0.5
            #     coordinates = [element.GetNodes()[0].X, element.GetNodes()[0].Y, element.GetNodes()[0].Z,
            #                    element.GetNodes()[1].X, element.GetNodes()[1].Y, element.GetNodes()[1].Z,
            #                    element.GetNodes()[2].X, element.GetNodes()[2].Y, element.GetNodes()[2].Z,]
            #     # print("ELEMENT:", coordinates)

            #     Jac_xi = DNXI * coordinates
            #     Jac_eta = DNETA * coordinates
            #     J_total = KratosMultiphysics.Matrix([Jac_xi,Jac_eta])
            #     J_total = J_total.transpose()
            #     # print(element.Id, "J_xi: ", Jac_xi)
            #     # print(element.Id, "J_eta: ", Jac_eta)
            #     # print(element.Id, "J_tot: ", J_total)
            #     J_X = element.GetGeometry().Jacobian(0)
            #     print(element.Id, "J_X : ", J_X)
            #     for node in element.GetNodes():
            #         # node.X = node.X0
            #         # node.Y = node.Y0
            #         node.Z = 0.0
            #     J_X0 = element.GetGeometry().Jacobian(0)
            #     # print(element.Id, "J_X : ", J_X)
            #     print(element.Id, "J_X0: ", J_X0)
            #     # J_diff = J_total - J_X
            #     # print(element.Id, "J_diff = ", J_diff)

            # Calculate Jacobian for CURRENT flat configuration
            


            
            # new_imposed_strain = KratosMultiphysics.Vector(3)
            # new_imposed_strain[0] = 0.01
            # new_imposed_strain[1] = 0.01
            # new_imposed_strain[2] = 0.0
            # for element in self._GetSolver().GetComputingModelPart().Elements:
            #     element.SetValue(KratosMultiphysics.INITIAL_STRAIN_VECTOR, new_imposed_strain)

####################################################################################
# import numpy as np
# def ModifyInitialGeometry(self):
#     var_utils = KratosMultiphysics.VariableUtils()
#     self.initial_unmodified_coordinates = var_utils.GetInitialPositionsVector(mp.Nodes, 3)
    
#     self.deformed_config_jacobians = []
#     for elem in elements:
#         J = np.array(elem.GetGeometry().Jacobian())
#         self.deformed_config_jacobians.append(J)
        
#     #calculate normals
    
#     for node in mp.Nodes:
#         node.Z0 = 0.0
#         node.Z = 0.0
#     self.flattened_coordinates = var_utils.GetInitialPositionsVector(mp.Nodes, 3)
#     self.initial_displacements = self.initial_unmodified_coordinates - self.deformed_config_jacobians
    
#     self.flatted_config_jacobians = []
#     for elem in elements:
#         J = np.array(elem.GetGeometry().Jacobian())
#         self.flatted_config_jacobians.append(J)
        
#     Fs = []
#     for J,J0 in zip(self.deformed_config_jacobians,self.flatted_config_jacobians):
#         J0_mod #append a column with h/h0
#         J_mod
#         F3D = J_mod@np.inv(J0_mod)
#         F = F3D[0:2,0:2]
        
#         Fs.append(F)
        


# element.SetValue[VARIABLE, [0,0,0]]
####################################################################################
# "loads_process_list"       : [{
#     "python_module" : "assign_vector_by_direction_to_condition_process",
#     "kratos_module" : "KratosMultiphysics",
#     "check"         : "DirectorVectorNonZero direction",
#     "process_name"  : "AssignVectorByDirectionToConditionProcess",
#     "Parameters"    : {
#         "model_part_name" : "Structure.PointLoad3D_Load_on_points_Auto1",
#         "variable_name"   : "POINT_LOAD",
#         "interval"        : [0.0,"End"],
#         "modulus"         : 1000.0,
#         "direction"       : [0.0,0.0,-1.0]
#     }
# }],
####################################################################################