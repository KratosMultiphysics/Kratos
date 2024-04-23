# function2
from re import X
import sys
import time
import importlib
import numpy as np
import math
# from gid_output_process import GiDOutputProcess

import KratosMultiphysics
from KratosMultiphysics import *


def Create_Fluid_model_part(divisions):
    current_model = KratosMultiphysics.Model()
    model_part = current_model.CreateModelPart("ModelPart")
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X_GRADIENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
    model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

    problem_domain = KratosMultiphysics.Quadrilateral2D4(
        KratosMultiphysics.Node(1, 1.0, 1.0, 0.0),
        KratosMultiphysics.Node(2, 1.0,  2.0, 0.0),
        KratosMultiphysics.Node(3,  2.0,  2.0, 0.0),
        KratosMultiphysics.Node(4,  2.0, 1.0, 0.0))
    parameters = KratosMultiphysics.Parameters("{}")
    parameters.AddEmptyValue("element_name").SetString("Element2D3N")
    # parameters.AddEmptyValue("condition_name").SetString("LineCondition2D2N")
    parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
    parameters.AddEmptyValue("number_of_divisions").SetInt(divisions)

    KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, model_part, parameters).Execute()
    return model_part

def Import_Background_model_part(Model,name_background_mpda):
    # Set skin_model_part geometry
    # current_model = KratosMultiphysics.Model()
    model_part = Model.CreateModelPart("fluid_model_part")
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE_GRADIENT)
    model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
    KratosMultiphysics.ModelPartIO(name_background_mpda).ReadModelPart(model_part)
    return model_part

def Import_Structural_model_part(name_structural_mpda):
    # Set skin_model_part geometry
    current_model = KratosMultiphysics.Model()
    skin_model_part = current_model.CreateModelPart("skin_model_part")
    skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
    skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X)
    skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Y)
    KratosMultiphysics.ModelPartIO(name_structural_mpda).ReadModelPart(skin_model_part)
    return skin_model_part

# def Import_Structural_SUB_model_part(model_part, name_structural_mpda):
#     # Set skin_model_part geometry
#     skin_model_part = model_part.CreateSubModelPart("skin_model_part")
#     skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
#     KratosMultiphysics.ModelPartIO(name_structural_mpda).ReadModelPart(skin_model_part)
#     return skin_model_part, model_part

def FindSurrogateNodes(model_part,iter):
    start_time = time.time()
    # General methods for find the surrogate nodes
    name_surrogate_sub_model_part = "surrogate_sub_model_part" + "_" + str(iter)
    surrogate_sub_model_part = model_part.CreateSubModelPart(name_surrogate_sub_model_part)
    # second_level_surrogate_sub_model_part = model_part.CreateSubModelPart("second_level_surrogate_sub_model_part")
    count2 = 0
    # # Inside-Outside problem
    # for node in model_part.Nodes:
    #     a = - node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
    #     node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, a)
    for elem in model_part.Elements :
        count_pos = 0
        count_neg = 0
        for node in elem.GetGeometry() :
            phi = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if phi > 0 :
                count_pos = count_pos + 1
            else :
                count_neg = count_neg + 1
        # When count_neg*count_pos != 0 --> the element is cut
        if count_neg * count_pos != 0 :
            # The element is cut
            for node in elem.GetGeometry() :
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0 :
                    node.Set(BOUNDARY, True)
        if count_pos == 3 :
            # The element is a fluid element completely ouside the surrogate boundary
            elem.Set(MARKER,True)
            # elem.Set(ACTIVE,True)
            count2 = count2 + 1
        # if count_neg > 0 :
            # elem.Set(ACTIVE,False)
            # model_part.RemoveElement(elem)
    print('Total number of "Fluid" element : ', count2)
    # Count the number of surrogate nodes
    tot_sur_nodes = 0
    for node in model_part.Nodes :
        if node.Is(BOUNDARY):
            surrogate_sub_model_part.AddNode(node,0)
            tot_sur_nodes = tot_sur_nodes +1 
    print('Number of surrogate nodes: ',tot_sur_nodes)
    # Discretize if an element is "boundary" so if it has at least one node that is surrogate
    for elem in model_part.Elements :
        if elem.Is(MARKER) :
            for node in elem.GetGeometry() :
                if node.Is(BOUNDARY) :
                    elem.Set(BOUNDARY, True)
                    break
    print("--> %s seconds for Find_surrogate_nodes" % (time.time() - start_time))
    print('ciao!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PERCHE\' NON POSSO DE-ACTIVATE?')
    return surrogate_sub_model_part,tot_sur_nodes



def FindClosestSkinElement(model_part,skin_model_part,tot_sur_nodes,tot_skin_el) :
    start_time = time.time()
    # if we are interested in the closest skin ELEMENT
    # file_due = open("closest_skin_element.txt", "w")
    # Inizializzo array closest_element
    closest_element = [0] * tot_sur_nodes   # the number of element and nodes of 
                                            # the skin is the same in this case
    i = 0
    for node in model_part.Nodes :
        if node.Is(BOUNDARY):
            # Inizializzo array equation_el
            equation_el = [0] * tot_skin_el
            for j in range(tot_skin_el) :  # Run over the skin ELEMENTS
                node1 = skin_model_part.Conditions[j+1].GetNodes()[0]
                node2 = skin_model_part.Conditions[j+1].GetNodes()[1]
                # Distance squared from the element
                equation_el[j] = (node1.X-node.X)**2  + (node1.Y-node.Y)**2 + (node2.X-node.X)**2 + (node2.Y-node.Y)**2
            index_min = np.argmin(equation_el) # j-esimo nodo
            closest_element[i] = skin_model_part.Conditions[index_min+1].Id
            # file_due.write(str(closest_element[i]))
            # file_due.write('\n')
            i = i + 1
    # file_due.close()
    print("--> %s seconds for FindClosestSkinElement" % (time.time() - start_time))
    return closest_element



def Find_projections(model_part,skin_model_part,tot_sur_nodes,closest_element) :
    start_time = time.time()
    # Find the PROJECTION for each of the surrogate nodes on the closest skin element
    # 1.0.1 Initialize a matrix with the coordinates of the projections
    projection_surr_nodes = [[0 for _ in range(2)] for _ in range(tot_sur_nodes)]
    file_tre = open("projection_surr_nodes.txt", "w")
    # 1.1 Run over each surrogate node take the closest skin element
    i = 0
    for node in model_part.Nodes :
        if node.Is(BOUNDARY) :
            # 1.2 take the closest skin element --> closest_element[i]
            # 1.3 Get the two nodes
            node1 = skin_model_part.Conditions[closest_element[i]].GetNodes()[0]
            node2 = skin_model_part.Conditions[closest_element[i]].GetNodes()[1]
            # 1.4 Compute m and q 
            if (node1.X - node2.X) != 0 :
                m = (node1.Y - node2.Y) / (node1.X - node2.X)
                q = node1.Y - m * node1.X
                # 1.5 Compute Q
                if m != 0 :
                    Q = node.Y + 1/m * node.X
                    projection_surr_nodes[i][0] = (Q-q) / (m+1/m)
                else : 
                    # La retta perpendicolare è del tipo x = node.X
                    projection_surr_nodes[i][0] = node.X
                projection_surr_nodes[i][1] = m * projection_surr_nodes[i][0] + q
            else :
                # La retta per i nodi 1 e 2 è verticale -> Trovo subito le proiezioni
                # print('Attenzione, skin element VERTICALE \n\n\n\n\n')
                projection_surr_nodes[i][0] = node1.X
                projection_surr_nodes[i][1] = node.Y
            # 1.7 Need to check if the point actually lies on the closest elements: compute the distance from
                # each skin node of the closest element and check that is less than the length of the element.
            check1 = (projection_surr_nodes[i][0] - node1.X)**2 + (projection_surr_nodes[i][1] - node1.Y)**2
            check2 = (projection_surr_nodes[i][0] - node2.X)**2 + (projection_surr_nodes[i][1] - node2.Y)**2
            element_length = (node1.X - node2.X)**2 + (node1.Y - node2.Y)**2
            if check1 > element_length or check2 > element_length :
                # print('-->\nNeed projection correction for node: ', node.Id)
                # Need to find the real projection
                # Take the closest node between node1 and node2
                if check1 < check2 :
                    candidate = node1
                else :
                    candidate = node2
                #-----------------------NEED TO IMPROVE--------------------------------------- 
                # search for the element with node = candidate (so I find the other unique possible second-closest element)
                for cond in skin_model_part.Conditions:
                    for nod in cond.GetGeometry():
                        if nod.Id == candidate.Id :
                            break
                    if nod.Id == candidate.Id :
                        break
                # Now we know on which element the projection lies --> We take the two nodes: node1 & node2
                node1 = skin_model_part.Conditions[cond.Id].GetNodes()[0]
                node2 = skin_model_part.Conditions[cond.Id].GetNodes()[1]
                # 1.4 Compute m and q 
                if (node1.X - node2.X) != 0 :
                    m = (node1.Y - node2.Y) / (node1.X - node2.X)
                    q = node1.Y - m * node1.X
                    # 1.5 Compute Q
                    if m != 0 :
                        Q = node.Y + 1/m * node.X
                        projection_surr_nodes[i][0] = (Q-q) / (m+1/m)
                    else : 
                        # La retta perpendicolare è del tipo x = node.X
                        projection_surr_nodes[i][0] = node.X
                    projection_surr_nodes[i][1] = m * projection_surr_nodes[i][0] + q
                else :
                    # La retta per i nodi 1 e 2 è verticale -> Trovo subito le proiezioni
                    projection_surr_nodes[i][0] = node1.X
                    projection_surr_nodes[i][1] = node.Y
                # SECOND CHECK (the projection, actually lies on an element?)
                check1 = (projection_surr_nodes[i][0] - node1.X)**2 + (projection_surr_nodes[i][1] - node1.Y)**2
                check2 = (projection_surr_nodes[i][0] - node2.X)**2 + (projection_surr_nodes[i][1] - node2.Y)**2
                element_length = (node1.X - node2.X)**2 + (node1.Y - node2.Y)**2
                if check1 > element_length or check2 > element_length :
                    # print('Need a second projection correction for node: ', node.Id)
                    # No, the projection just found does not lie on a skin element
                    # --> Take the closest node as the projection
                    projection_surr_nodes[i][0] = candidate.X
                    projection_surr_nodes[i][1] = candidate.Y
                    node.Set(INTERFACE, True)
                else :
                    # We have found the correct projection
                    print('Trovata proiezione nel secondo elemento più vicino del nodo : ', node.Id)
            file_tre.write(str(projection_surr_nodes[i][0]))
            file_tre.write('  ')
            file_tre.write(str(projection_surr_nodes[i][1]))
            file_tre.write('\n')
            i = i + 1
    file_tre.close()
    print("--> %s seconds for Find_projections" % (time.time() - start_time))
    return projection_surr_nodes



def Dirichlet_BC (model_part,skin_model_part,tot_sur_nodes,closest_element,projection_surr_nodes) :
    surr_BC = [0] * tot_sur_nodes
    i = 0
    for node in model_part.Nodes :
        if node.Is(BOUNDARY) :
            # 1.1 Initialize a Kratos vector
            d_vector = KratosMultiphysics.Array3()
            # 1.2 Compute the vector x - x_tilde
            d_vector[0] = projection_surr_nodes[i][0] - node.X
            d_vector[1] = projection_surr_nodes[i][1] - node.Y
            d_vector[2] = 0.0
            # 1.3 Take the gradient of the generic scalar vector field
            grad = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)
            # 1.4 Compute the scalar product
            correction = 0
            for j in range(2) :
                correction = correction +  d_vector[j] * grad[j]
            # 1.5 Obtain the value of the scalar field at the projection with an INTERPOLATION
            node1 = skin_model_part.Conditions[closest_element[i]].GetNodes()[0]
            node2 = skin_model_part.Conditions[closest_element[i]].GetNodes()[1]
            d_star = math.sqrt((node1.X-projection_surr_nodes[i][0])**2 + (node1.Y-projection_surr_nodes[i][1])**2)
            d = math.sqrt((node1.X-node2.X)**2 + (node1.Y-node2.Y)**2)
            # 1.5.0.1 Compute the dirichlet value at the projection
            dir_value_projection = node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - d_star / d * (node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-node2.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
            # 1.5.1 Additional check for the node with "second correction"
            if node.Is(INTERFACE) :
                if projection_surr_nodes[i][0] == node1.X and projection_surr_nodes[i][1] == node1.Y:
                    dir_value_projection = node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                    # print('node 1')
                else :
                    # print(projection_surr_nodes[i][0], node2.X)
                    dir_value_projection = node2.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                    # print('node 2')
            # 1.6 Compute the approximated Dirichlet BC at the surrogate boundary
            surr_BC[i] = dir_value_projection - correction
            # print('errore = ', abs(velocity - correction - node.X-node.Y))
            # print('\n')
            i = i+1
    return surr_BC


def Interpolation(skin_model_part,closest_element,projection_surr_nodes, i, node):
    node1 = skin_model_part.Conditions[closest_element[i]].GetNodes()[0]
    node2 = skin_model_part.Conditions[closest_element[i]].GetNodes()[1]
    d_star = math.sqrt((node1.X-projection_surr_nodes[i][0])**2 + (node1.Y-projection_surr_nodes[i][1])**2)
    d = math.sqrt((node1.X-node2.X)**2 + (node1.Y-node2.Y)**2)
    dir_value_projection = node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - d_star / d * (node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-node2.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
    # 1.5.1 Additional check for the node with "second correction"
    if node.Is(INTERFACE) :
        if projection_surr_nodes[i][0] == node1.X and projection_surr_nodes[i][1] == node1.Y:
            dir_value_projection = node1.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
        else :
            dir_value_projection = node2.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
    return dir_value_projection



def Dirichlet_BC_CFD (model_part,skin_model_part,tot_sur_nodes,closest_element,projection_surr_nodes) :
    # Compute the value of the velocity_x on the surrogate boundary using the grandient and the distance
    surr_BC_x = [0] * tot_sur_nodes
    surr_BC_y = [0] * tot_sur_nodes
    i = 0
    for node in model_part.Nodes :
        if node.Is(BOUNDARY) :
            # Need the distance vector : x - x_tilde
            # 1.1 Initialize a Kratos vector
            d_vector = KratosMultiphysics.Array3()
            # 1.2 Compute the vector x - x_tilde
            d_vector[0] = projection_surr_nodes[i][0] - node.X
            d_vector[1] = projection_surr_nodes[i][1] - node.Y
            d_vector[2] = 0.0
            # 1.3 Take the gradient of the generic scalar vector field
            grad_x = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X_GRADIENT)
            grad_y = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y_GRADIENT)
            # 1.4 Compute the scalar product
            scalar_product_x = 0
            scalar_product_y = 0
            for j in range(2) :
                scalar_product_x = scalar_product_x +  d_vector[j] * grad_x[j]
                scalar_product_y = scalar_product_y +  d_vector[j] * grad_y[j]
            correction_x = scalar_product_x
            correction_y = scalar_product_y
            # 1.5 Obtain the value of the scalar field at the projection with an INTERPOLATION
            node1 = skin_model_part.Conditions[closest_element[i]].GetNodes()[0]
            node2 = skin_model_part.Conditions[closest_element[i]].GetNodes()[1]
            d_star = math.sqrt((node1.X-projection_surr_nodes[i][0])**2 + (node1.Y-projection_surr_nodes[i][1])**2)
            d = math.sqrt((node1.X-node2.X)**2 + (node1.Y-node2.Y)**2)
            velocity_x = node1.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X) - d_star / d * (node1.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)-node2.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X))
            velocity_y = node1.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y) - d_star / d * (node1.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)-node2.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y))
            # 1.5.1 Additional check for the node with "second correction"
            if node.Is(INTERFACE) :
                if projection_surr_nodes[i][0] == node1.X and projection_surr_nodes[i][1] == node1.Y:
                    velocity_x = node1.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                    velocity_y = node1.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    # print('node 1')
                else :
                    # print(projection_surr_nodes[i][0], node2.X)
                    velocity_x = node2.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                    velocity_y = node2.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                    # print('node 2')
            # 1.6 Compute the approximated Dirichlet BC at the surrogate boundary
            surr_BC_x[i] = velocity_x - correction_x
            surr_BC_y[i] = velocity_y - correction_y
            # print('errore = ', abs(velocity - correction - node.X-node.Y))
            # print('\n')
            i = i+1
    return surr_BC_x, surr_BC_y




def Create_sub_model_part_fluid(main_model_part,iter) :
    name_sub_model_part_fluid = "sub_model_part_fluid" + "_" + str(iter)
    sub_model_part_fluid = main_model_part.CreateSubModelPart(name_sub_model_part_fluid)
    for node in main_model_part.Nodes :
        if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0 :
            sub_model_part_fluid.AddNode(node,0)
    for elem in main_model_part.Elements :
        if elem.Is(MARKER):
            sub_model_part_fluid.AddElement(elem,0)
    return sub_model_part_fluid



def Compute_error(main_model_part, sub_model_part_fluid) :
        file_due = open("error.txt", "w")
        L2_err = 0
        L2_err_area = 0
        L2_grad_err = 0
        L2_grad_err_area = 0
        KratosMultiphysics.ComputeNodalGradientProcess(
        sub_model_part_fluid,
        KratosMultiphysics.TEMPERATURE,
        KratosMultiphysics.TEMPERATURE_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()
        
        # for node in main_model_part.Nodes : 
        #     if node.Id == 198 :
        #         print(node.GetValue(NODAL_AREA))
        #     if node.Is(BOUNDARY) :
        #         node.SetValue(NODAL_AREA, 0)
        #     if node.Id == 198 :
        #         print(node.GetValue(NODAL_AREA))
        # for elem in main_model_part.Elements :
        #     if elem.Is(ACTIVE) :
        #         for node in elem.GetGeometry() :
        #             node.SetValue(NODAL_AREA, node.GetValue(NODAL_AREA) + elem.GetGeometry().Area() / 3)
        total_number_fluid_nodes = 0
        total_area = 0
        max_err = 0
        for node in main_model_part.Nodes :
            exact_grad = KratosMultiphysics.Array3()
            # exact = node.X + node.Y    # --> Tutto lineare
            # exact_grad[0] = 1.0        # --> Tutto lineare
            # exact_grad[1] = 1.0        # --> Tutto lineare
            # exact = 0.25*(9 - ((node.X)**2 + (node.Y)**2) )  # --> Paraboloide
            # exact_grad[0] = 0.25*(-2*node.X)                 # --> Paraboloide
            # exact_grad[1] = 0.25*(-2*node.Y)                 # --> Paraboloide

            # log classico        
            # exact = 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log((node.X)**2+(node.Y)**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y)
            # exact_grad[0] = 0.25 * (-2*node.X + 2*node.X / (node.X**2 + node.Y**2))  +  0.25 * math.cos(node.X) * math.sinh(node.Y)
            # exact_grad[1] = 0.25 * (-2*node.Y + 2*node.Y / (node.X**2 + node.Y**2))  +  0.25 * math.sin(node.X) * math.cosh(node.Y)
            
            # sin(x)*cos(y)
            # exact = math.sin(node.X) * math.cos(node.Y)
            # exact_grad[0] = math.cos(node.X) * math.cos(node.Y)
            # exact_grad[1] = -math.sin(node.X) * math.sin(node.Y)

            # sin(x)*cos(y)+log(1+x**2+y**2)
            exact = math.sin(node.X)*math.cos(node.Y)+math.log(1+node.X**2+node.Y**2)
            exact_grad[0] =0
            exact_grad[1] = 0

            ## Senza Logaritmo
            # exact = 0.25*(9-node.X**2-node.Y**2-2*math.log(3)) + 0.25 *math.sin(node.X) * math.sinh(node.Y) 
            # exact_grad[0] = 0.25 * (-2*node.X)  +  0.25 * math.cos(node.X) * math.sinh(node.Y)
            # exact_grad[1] = 0.25 * (-2*node.Y)  +  0.25 * math.sin(node.X) * math.cosh(node.Y)

            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0:
                if abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact) > max_err :
                    max_err = abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)
                total_number_fluid_nodes = total_number_fluid_nodes + 1
                nodal_area = node.GetValue(NODAL_AREA)
                L2_err = L2_err + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)**2
                L2_err_area = L2_err_area + nodal_area * (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)**2
                L2_grad_err = L2_grad_err + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]-exact_grad[0])**2 + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]-exact_grad[1])**2 
                L2_grad_err_area = L2_grad_err_area + nodal_area * ((node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]-exact_grad[0])**2 + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]-exact_grad[1])**2 )
                total_area = total_area + nodal_area
                file_due.write(str(abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)))
                file_due.write('       ')
                file_due.write(str(abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]-exact_grad[0])))
                file_due.write('       ')
                file_due.write(str(abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]-exact_grad[1])))
                file_due.write('       ')
                file_due.write(str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]))
                file_due.write('       ')
                file_due.write(str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]))
            else :
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
            file_due.write('\n')
            
        if total_number_fluid_nodes == 0 :
            total_number_fluid_nodes = len(sub_model_part_fluid.Nodes)
        L2_err = math.sqrt(L2_err  / total_number_fluid_nodes )
        L2_err_area = math.sqrt(L2_err_area / total_area) 
        H1_err = L2_err +  math.sqrt( L2_grad_err  / total_number_fluid_nodes )
        H1_err_area = L2_err_area + math.sqrt(L2_grad_err_area / total_area ) 
        # file_due.write(str(L2_err))
        print('Errore in norma L2 (equal areas): ', L2_err)
        print('Errore in norma L2 : ', L2_err_area)
        # print('Errore in norma H1 (equal areas): ', H1_err)
        print('Errore in norma H1 : ', H1_err_area)
        file_due.close
        return L2_err_area, H1_err_area, max_err



def ComputeGradientCoefficients (sub_model_part_fluid, model,surrogate_sub_model_part) :
    Parameters = KratosMultiphysics.Parameters("""
    {
        "model_part_name" : "ThermalModelPart",
        "boundary_sub_model_part_name" : "surrogate_sub_model_part",
        "sbm_interface_condition_name" : "LineCondition2D2N",
        "conforming_basis" : true,
        "extension_operator_type" : "MLS",
        "mls_extension_operator_order" : 1,
        "levelset_variable_name" : "DISTANCE"
    }
    """)

    # CALCOLO DEI COEFFICIENTI PER IL CALCOLO DEL GRADIENTE
    a = KratosMultiphysics.ShiftedBoundaryMeshlessInterfaceUtilityCopy2(model,Parameters)
    # a = KratosMultiphysics.ShiftedBoundaryMeshlessInterfaceUtility(model,Parameters)
    result = a.SetSurrogateBoundaryNodalGradientWeights()

    # CHECK IF THE ARE VERY_PROBLEMATIC ELEMENTS
    # number_very_problematic = 0
    # for elem in sub_model_part_fluid.Elements :
    #     if elem.Is(BOUNDARY) :
    #         count = 0
    #         for node in elem.GetGeometry() :
    #             if node.Is(BOUNDARY) :
    #                 count = count + 1
            # if count == 3 :
            #     number_very_problematic = number_very_problematic + 1
                    

    # print('Number of surr nodes : ', len(surrogate_sub_model_part.Nodes))
    # print('Number of result2 : ', len(result2)- number_very_problematic)
    # print('Number of VERY_PROBLEMATIC nodes : ', number_very_problematic)
    # return result, result2
    return result

def Compute_T_matrix (result, node, projection_surr_nodes, i) :
    my_result = result[node.Id]
    T_tilde = [[0 for _ in range(len(my_result))] for _ in range(2)]
    j = 0
    for key, value in my_result.items():
        T_tilde[0][j] = value[0]
        T_tilde[1][j] = value[1]     
        j = j + 1
    # Get the distance vector using the projections
    d_vector = [[0 for _ in range(1)] for _ in range(2)]
    d_vector[0] = projection_surr_nodes[i][0] - node.X
    d_vector[1] = projection_surr_nodes[i][1] - node.Y
    # MATRIX MULTIPLICATION
    T = np.matmul( np.transpose(T_tilde) , d_vector )
    return T

def Impose_MPC_Globally (main_model_part, result, skin_model_part, closest_element, projection_surr_nodes, T, node, j) :
    # Interpolate the value at the projection 
    dirichlet_projection = Interpolation(skin_model_part, closest_element, projection_surr_nodes, j-1, node)
    # dirichlet_projection = 0.25*(9 - ((projection_surr_nodes[j-1][0])**2 + (projection_surr_nodes[j-1][1])**2) ) 
    # dirichlet_projection =   0.25*(9-(projection_surr_nodes[j-1][0])**2-(projection_surr_nodes[j-1][1])**2-2*math.log(3) + math.log((projection_surr_nodes[j-1][0])**2+(projection_surr_nodes[j-1][1])**2)) + 0.25 *math.sin(projection_surr_nodes[j-1][0]) * math.sinh(projection_surr_nodes[j-1][1])
    # dirichlet_projection = math.sin(node.X)*math.cos(node.Y)+math.log(1+node.X**2+node.Y**2)

    my_result = result[node.Id]
    # print('\n', node.Id)
    # print('\n')
    if len(my_result) != 0 :
        DofMasterVector = []
        if node.IsNot(SLAVE) :
            CoeffVector = KratosMultiphysics.Vector(len(my_result)-1)
            ConstantVector = 0.0*KratosMultiphysics.Vector(len(my_result)-1)
        else :
            CoeffVector = KratosMultiphysics.Vector(len(my_result))
            ConstantVector = 0.0*KratosMultiphysics.Vector(len(my_result))
        i = 0
        k = 0
        Coeff_Slave = 0
        for key, value in my_result.items() :
            # print(key)
            node_master = main_model_part.GetNode(key)
            # print(key)
            if node.IsNot(SLAVE) :
                # Need to find the "node" term and bring it to the left-hand-side
                if node.Id != key :
                    DofMasterVector.append(node_master.GetDof(KratosMultiphysics.TEMPERATURE))
                    ConstantVector[i] = dirichlet_projection
                    CoeffVector[i] = - T[k]
                    i = i + 1
                    k = k + 1
                else :
                    Coeff_Slave = - T[k]
                    k = k + 1
            else :
                DofMasterVector.append(node_master.GetDof(KratosMultiphysics.TEMPERATURE))
                ConstantVector[i] = dirichlet_projection
                CoeffVector[i] = - T[i]
                i = i + 1
        CoeffVector = CoeffVector / (1-Coeff_Slave)
        ConstantVector = ConstantVector / (1-Coeff_Slave)
        CoeffMatrix = KratosMultiphysics.Matrix(np.array(CoeffVector).reshape(1,-1))
        DofSlaveVector = [node.GetDof(KratosMultiphysics.TEMPERATURE)]
        # Create the constraint
        main_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", j, DofMasterVector, DofSlaveVector, CoeffMatrix, ConstantVector)
    else :
        # Create the constraint where simply grad == 0
        node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, dirichlet_projection)
        node.Fix(KratosMultiphysics.TEMPERATURE)
        print('grad == 0, OOOOOOOOOOOOOO, sabelo')
        # exit()
    return 0 