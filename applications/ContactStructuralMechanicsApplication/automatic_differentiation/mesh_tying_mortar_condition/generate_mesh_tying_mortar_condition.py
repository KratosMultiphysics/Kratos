
# Import KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.ContactStructuralMechanicsApplication

# Import sympy utils
import sympy
from KratosMultiphysics.ContactStructuralMechanicsApplication import custom_sympy_fe_utilities

do_simplifications = False
mode = "c" #to output to a c++ file

impose_partion_of_unity = False

dim_combinations = [2,2,2,2,3,3,3,3,3,3,3,3]
nnodeselement_combinations = [2,2,3,3,4,4,3,3,4,4]
nnodeselement_master_combinations = [2,2,3,3,4,4,4,4,3,3]
tensor_combinations = [1,2,1,2,1,3,1,3,1,3,1,3]

lhs_string = ""
lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate< >\ntemplate< >\nvoid MeshTyingMortarCondition<TDim,TNumNodesElem,TNumNodesElemMaster>::CalculateLocalLHS<MeshTyingMortarCondition<TDim,TNumNodesElem,TNumNodesElemMaster>::TTensor>(\n    Matrix& rLocalLHS,\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    const DofData<TTensor>& rDofData\n    )\n{\n    // We get the mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"

lhs_template_end_string = "\n}\n"

rhs_string = ""

rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\ntemplate<>\nvoid MeshTyingMortarCondition<TDim,TNumNodesElem, TNumNodesElemMaster>::CalculateLocalRHS<MeshTyingMortarCondition<TDim,TNumNodesElem,TNumNodesElemMaster>::TTensor>(\n    Vector& rLocalRHS,\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    const DofData<TTensor>& rDofData\n    )\n{\n    // Initialize values\n    const BoundedMatrix<double, TNumNodes, TTensor>& u1 = rDofData.u1;\n    const BoundedMatrix<double, TNumNodesMaster, TTensor>& u2 = rDofData.u2;\n\n    const BoundedMatrix<double, TNumNodes, TTensor>& lm = rDofData.LagrangeMultipliers; \n\n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n"

rhs_template_end_string = "\n}\n"

for dim, nnodeselement, nnodeselement_master, tensor in zip(dim_combinations, nnodeselement_combinations, nnodeselement_master_combinations, tensor_combinations):

    if ((nnodeselement == 3) or (dim == 2 and nnodeselement == 4)):
        nnodes =  2
    else:
        if (nnodeselement == 4):
            nnodes = 3
        else:
            nnodes = 4

    if ((nnodeselement_master == 3) or (dim == 2 and nnodeselement_master == 4)):
        nnodes_master =  2
    else:
        if (nnodeselement_master == 4):
            nnodes_master = 3
        else:
            nnodes_master = 4

    number_dof = tensor * (nnodes_master + 2 * nnodes)

    #Defining the unknowns
    u1 = custom_sympy_fe_utilities.DefineMatrix('u1',nnodes,tensor) #u1(i,j) is displacement of node i component j at domain 1
    u2 = custom_sympy_fe_utilities.DefineMatrix('u2',nnodes_master,tensor) #u2(i,j) is displacement of node i component j at domain 2
    lm = custom_sympy_fe_utilities.DefineMatrix('lm',nnodes,tensor)

    # Define test functions
    w1 = custom_sympy_fe_utilities.DefineMatrix('w1',nnodes,tensor)
    w2 = custom_sympy_fe_utilities.DefineMatrix('w2',nnodes_master,tensor)
    wlm = custom_sympy_fe_utilities.DefineMatrix('wlm',nnodes, tensor)

    DOperator = custom_sympy_fe_utilities.DefineMatrix('DOperator',nnodes,nnodes)
    MOperator = custom_sympy_fe_utilities.DefineMatrix('MOperator',nnodes,nnodes_master)

    # Define dofs & test function vector
    dofs = sympy.Matrix( sympy.zeros(number_dof, 1) )
    testfunc = sympy.Matrix( sympy.zeros(number_dof, 1) )
    count = 0

    for i in range(0,nnodes_master):
        for k in range(0,tensor):
            dofs[count] = u2[i,k]
            testfunc[count] = w2[i,k]
            count+=1
    for i in range(0,nnodes):
        for k in range(0,tensor):
            dofs[count] = u1[i,k]
            testfunc[count] = w1[i,k]
            count+=1
    for i in range(0,nnodes):
        for k in range(0,tensor):
            dofs[count] = lm[i,k]
            testfunc[count] = wlm[i,k]
            count+=1

    print("dofs = ",dofs)
    print("testfunc = ",testfunc)

    #############################################################################
    #############################################################################
    ########################## FUNCTIONAL DEFINITION ############################
    #############################################################################
    #############################################################################

    # Defining the residual
    Du1Mu2 = DOperator * u1 - MOperator * u2
    Dw1Mw2 = DOperator * w1 - MOperator * w2

    # Compute galerkin functional
    rv_galerkin = 0
    # Defining the functional
    for node in range(nnodes):
        rv_galerkin -= (lm.row(node) * (Dw1Mw2.row(node)).transpose())[0,0]
        rv_galerkin -= (wlm.row(node) * (Du1Mu2.row(node)).transpose())[0,0]

    if(do_simplifications):
        rv_galerkin = sympy.simplify(rv_galerkin)

    #############################################################################
    # Complete functional
    rv = sympy.Matrix( sympy.zeros(1, 1) )
    rv[0,0] = rv_galerkin

    rhs,lhs = custom_sympy_fe_utilities.Compute_RHS_and_LHS(rv.copy(), testfunc, dofs, False)
    print("LHS= ",lhs.shape)
    print("RHS= ",rhs.shape)
    print("LHS and RHS have been created!")

    lhs_out = custom_sympy_fe_utilities.OutputMatrix_CollectingFactors(lhs,"lhs", mode, 1, number_dof)
    rhs_out = custom_sympy_fe_utilities.OutputVector_CollectingFactors(rhs,"rhs", mode, 1, number_dof)
    print("Substitution strings are ready....")

    lhs_string += lhs_template_begin_string
    lhs_string += lhs_out
    lhs_string += lhs_template_end_string

    rhs_string += rhs_template_begin_string
    rhs_string += rhs_out
    rhs_string += rhs_template_end_string

    lhs_string = lhs_string.replace("TDim", str(dim))
    lhs_string = lhs_string.replace("TNumNodesElemMaster", str(nnodeselement_master))
    lhs_string = lhs_string.replace("TNumNodesElem", str(nnodeselement))
    if (tensor == 1):
        lhs_string = lhs_string.replace("TTensor", "ScalarValue")
    elif (tensor == 2):
        lhs_string = lhs_string.replace("TTensor", "Vector2DValue")
    elif (tensor == 3):
        lhs_string = lhs_string.replace("TTensor", "Vector3DValue")
    lhs_string = lhs_string.replace("TNumNodesMaster", str(nnodes_master))
    lhs_string = lhs_string.replace("TNumNodes", str(nnodes))
    lhs_string = lhs_string.replace("MatrixSize", str(lhs.shape[0]))

    rhs_string = rhs_string.replace("TDim", str(dim))
    rhs_string = rhs_string.replace("TNumNodesElemMaster", str(nnodeselement_master))
    rhs_string = rhs_string.replace("TNumNodesElem", str(nnodeselement))
    if (tensor == 1):
        rhs_string = rhs_string.replace("TTensor", "ScalarValue")
    elif (tensor == 2):
        rhs_string = rhs_string.replace("TTensor", "Vector2DValue")
    elif (tensor == 3):
        rhs_string = rhs_string.replace("TTensor", "Vector3DValue")
    rhs_string = rhs_string.replace("TNumNodesMaster", str(nnodes_master))
    rhs_string = rhs_string.replace("TNumNodes", str(nnodes))
    rhs_string = rhs_string.replace("MatrixSize", str(rhs.shape[0]))
    lhs_string = lhs_string.replace("lhs(", "rLocalLHS(")
    rhs_string = rhs_string.replace("rhs[", "rLocalRHS[")

#############################################################################
################################# FINAL SAVING ##############################
#############################################################################

input = open("mesh_tying_mortar_condition_template.cpp",'r').read()
outputstring = input.replace("// replace_lhs", lhs_string)
outputstring = outputstring.replace("// replace_rhs", rhs_string)
output = open("mesh_tying_mortar_condition.cpp",'w')
output.write(outputstring)
output.close()

print("Strings have been replaced...")

print("Process Finished..................")