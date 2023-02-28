import KratosMultiphysics
import json
import numpy as np

def GetStringIndexesAsInt(list_of_string_keys):
    res = [eval(i) for i in list_of_string_keys]
    res = np.array(res)
    return res


def main(OriginalNumberOfElements):

    try_flag_for_petrov_galerkin = True  #We start considering it might have a left basis
    try_flag_for_hrom = True  #We start considering it might have hrom info

    # Get the ROM settings from the RomParameters.json input file
    with open('RomParameters.json') as rom_parameters:
        rom_parameters = KratosMultiphysics.Parameters(rom_parameters.read())

    NodeIndexes = GetStringIndexesAsInt(rom_parameters["nodal_modes"].keys())



    number_of_nodes = np.max(NodeIndexes)
    right_nodal_basis = []
    left_nodal_basis = []
    for i in range(number_of_nodes):
        right_nodal_basis.append(np.array(rom_parameters["nodal_modes"][f"{i+1}"].GetMatrix(), copy=False))
        if try_flag_for_petrov_galerkin:
            try:
                left_nodal_basis.append(np.array(rom_parameters["petrov_galerkin_nodal_modes"][f"{i+1}"].GetMatrix(), copy=False))
            except:
                try_flag_for_petrov_galerkin=False

    right_nodal_basis = np.concatenate(right_nodal_basis, axis=0)


    #Node Ids
    Dimensions  = len(rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray())
    np.save('NodeIds.npy',  np.arange(1,((right_nodal_basis.shape[0]+1)/Dimensions), 1, dtype=int)   ) #this fixes the +1 issue !!!!!!!!!!

    #Right Basis
    print('the number of nodes is:', number_of_nodes, 'with',Dimensions,'dofs per node, a total number of dofs of:', number_of_nodes*Dimensions)
    print('the size of the modes matrix is: ', right_nodal_basis.shape)
    np.save('RightBasisMatrix.npy', right_nodal_basis)


    #Left Basis
    if try_flag_for_petrov_galerkin:
        left_nodal_basis = np.concatenate(left_nodal_basis, axis=0)
        print('the number of nodes is:', number_of_nodes, 'with',Dimensions,'dofs per node, a total number of dofs of:', number_of_nodes*Dimensions)
        print('the size of the modes matrix is: ', left_nodal_basis.shape)
        np.save('LeftBasisMatrix.npy', left_nodal_basis)





    try:
        ElementsVector = GetStringIndexesAsInt(rom_parameters["elements_and_weights"]["Elements"].keys())
    except:
        try_flag_for_hrom = False

    if try_flag_for_hrom:
        WeightsMatrix = []
        for ElementIndex in ElementsVector:
            WeightsMatrix.append(rom_parameters["elements_and_weights"]["Elements"][str(ElementIndex)].GetDouble())


        ConditionsVector = GetStringIndexesAsInt(rom_parameters["elements_and_weights"]["Conditions"].keys())
        for ConditionIndex in ConditionsVector:
            WeightsMatrix.append(rom_parameters["elements_and_weights"]["Conditions"][str(ConditionIndex)].GetDouble())

        np.save('WeightsMatrix.npy',np.block([WeightsMatrix]) )
        np.save('Elementsvector.npy', np.r_[ElementsVector,ConditionsVector+OriginalNumberOfElements])




if __name__=='__main__':

    OriginalNumberOfElements = 116

    main(OriginalNumberOfElements)


