from sympy import *
from custom_sympy_fe_utilities import *
import operator 

do_simplifications = False
mode = "c" #to output to a c++ file

separate_derivatives = True
impose_partion_of_unity = False

dim_combinations = [2,2,2,2,3,3,3,3]
nnodeselement_combinations = [3,3,4,4,4,4,8,8]
tensor_combinations = [1,2,1,2,1,3,1,3]

lhs_string = ""
lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate< >\ntemplate< >\nboost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize> MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLocalLHS<MatrixSize>(\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    DofData& rDofData,\n    const unsigned int& rMasterElementIndex,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    boost::numeric::ublas::bounded_matrix<double, MatrixSize, MatrixSize> lhs;\n\n    // We get the mortar operators\n    const boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> DOperator = rMortarConditionMatrices.DOperator;\n\n"

lhs_template_end_string = "\n\n    return lhs;\n}\n"

rhs_string = ""
rhs_template_begin_scalar_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\ntemplate<>\narray_1d<double, MatrixSize> MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLocalRHS<MatrixSize>(\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    DofData& rDofData,\n    const unsigned int& rMasterElementIndex,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    array_1d<double,MatrixSize> rhs;\n\n    // Master segment info\n    GeometryType& CurrentMasterElement = mThisMasterConditions[rMasterElementIndex]->GetGeometry();\n\n    // Initialize values\n    const array_1d<double, NumNodes> u1 = ContactUtilities::GetVariableVector<NumNodes>(this->GetGeometry(), TEMPERATURE, 0);\n    const array_1d<double, NumNodes> u2 = ContactUtilities::GetVariableVector<NumNodes>(CurrentMasterElement, TEMPERATURE, 0);\n\n    const array_1d<double, NumNodes> lm = ContactUtilities::GetVariableVector<NumNodes>(this->GetGeometry(), SCALAR_LAGRANGE_MULTIPLIER, 0); \n\n    // Mortar operators\n    const bounded_matrix<double, NumNodes, NumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, NumNodes, NumNodes> DOperator = rMortarConditionMatrices.DOperator;\n\n"

rhs_template_begin_vector_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\ntemplate<>\narray_1d<double, MatrixSize> MeshTyingMortarCondition<TDim,TNumNodesElem,TTensor>::CalculateLocalRHS<MatrixSize>(\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    DofData& rDofData,\n    const unsigned int& rMasterElementIndex,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    array_1d<double,MatrixSize> rhs;\n\n    // Master segment info\n    GeometryType& CurrentMasterElement = mThisMasterConditions[rMasterElementIndex]->GetGeometry();\n\n    // Initialize values\n    const bounded_matrix<double, NumNodes, TDim> u1 = ContactUtilities::GetVariableMatrix<TDim,NumNodes>(this->GetGeometry(), DISPLACEMENT, 0);\n    const bounded_matrix<double, NumNodes, TDim> u2 = ContactUtilities::GetVariableMatrix<TDim,NumNodes>(CurrentMasterElement, DISPLACEMENT, 0);\n\n    const bounded_matrix<double, NumNodes, TDim> lm = ContactUtilities::GetVariableMatrix<TDim,NumNodes>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); \n\n    // Mortar operators\n    const bounded_matrix<double, NumNodes, NumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, NumNodes, NumNodes> DOperator = rMortarConditionMatrices.DOperator;\n\n"

rhs_template_end_string = "\n\n    return rhs;\n}\n"

for dim, nnodeselement, tensor in zip(dim_combinations, nnodeselement_combinations, tensor_combinations):
    
    if ((nnodeselement == 3) or (dim == 2 and nnodeselement == 4)):
        nnodes =  2
    else:
        if (nnodeselement == 4):
            nnodes = 3
        else:
            nnodes = 4
            
    number_dof = tensor * (3 * nnodes)

    if tensor == 1:
        #Defining the unknowns
        u1 = DefineVector('u1',nnodes) #u1(i,j) is displacement of node i component j at domain 1
        u2 = DefineVector('u2',nnodes) #u2(i,j) is displacement of node i component j at domain 2
        lm = DefineVector('lm',nnodes) 
        
        # Define test functions
        w1 = DefineVector('w1',nnodes)
        w2 = DefineVector('w2',nnodes)
        wlm = DefineVector('wlm',nnodes)
    else:
        #Defining the unknowns
        u1 = DefineMatrix('u1',nnodes,dim) #u1(i,j) is displacement of node i component j at domain 1
        u2 = DefineMatrix('u2',nnodes,dim) #u2(i,j) is displacement of node i component j at domain 2
        lm = DefineMatrix('lm',nnodes,dim) 
        
        # Define test functions
        w1 = DefineMatrix('w1',nnodes,dim)
        w2 = DefineMatrix('w2',nnodes,dim)
        wlm = DefineMatrix('wlm',nnodes, dim)
            
    DOperator = DefineMatrix('DOperator',nnodes,nnodes) 
    MOperator = DefineMatrix('MOperator',nnodes,nnodes) 

    # Define dofs & test function vector
    dofs = Matrix( zeros(number_dof, 1) )
    testfunc = Matrix( zeros(number_dof, 1) )
    count = 0
    if (tensor == 1):
        for i in range(0,nnodes):
                dofs[count] = u2[i]
                testfunc[count] = w2[i]
                count+=1
        for i in range(0,nnodes):
                dofs[count] = u1[i]
                testfunc[count] = w1[i]
                count+=1
        for i in range(0,nnodes):
            dofs[count] = lm[i]
            testfunc[count] = wlm[i]
            count+=1
    else:
        for i in range(0,nnodes):
            for k in range(0,dim):
                dofs[count] = u2[i,k]
                testfunc[count] = w2[i,k]
                count+=1
        for i in range(0,nnodes):
            for k in range(0,dim):
                dofs[count] = u1[i,k]
                testfunc[count] = w1[i,k]
                count+=1
        for i in range(0,nnodes):
            for k in range(0,dim):
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

    # Compute galerkin functional 
    rv_galerkin = 0
    for node in range(nnodes):
        if (tensor == 1):
            # Defining the normal gap
            Du1Mu2 = DOperator * u1 - MOperator * u2
            Dw1Mw2 = DOperator * w1 - MOperator * w2
            # Defining the functional
            rv_galerkin -= Dw1Mw2[node] * lm[node]
            rv_galerkin -= Du1Mu2[node] * wlm[node]
        else:
            for dvalue in range(dim):
                # Defining the normal gap
                Du1Mu2 = DOperator * u1.col(dvalue) - MOperator * u2.col(dvalue)
                Dw1Mw2 = DOperator * w1.col(dvalue) - MOperator * w2.col(dvalue)
                # Defining the functional
                rv_galerkin -= Dw1Mw2[node] * lm[node, dvalue]
                rv_galerkin -= Du1Mu2[node] * wlm[node, dvalue]

    if(do_simplifications):
        rv_galerkin = simplify(rv_galerkin)

    #############################################################################
    # Complete functional
    rv = Matrix( zeros(1, 1) )
    rv[0,0] = rv_galerkin 

    rhs,lhs = Compute_RHS_and_LHS(rv.copy(), testfunc, dofs, False)
    print("LHS= ",lhs.shape)
    print("RHS= ",rhs.shape)
    print("LHS and RHS have been created!")

    lhs_out = OutputMatrix_CollectingFactors(lhs,"lhs", mode, 1, number_dof)
    rhs_out = OutputVector_CollectingFactors(rhs,"rhs", mode, 1, number_dof)
    print("Substitution strings are ready....")
    
    lhs_string += lhs_template_begin_string
    lhs_string += lhs_out
    lhs_string += lhs_template_end_string
    
    if (tensor == 1):
        rhs_string += rhs_template_begin_scalar_string
    else:
        rhs_string += rhs_template_begin_vector_string
    rhs_string += rhs_out
    rhs_string += rhs_template_end_string

    lhs_string = lhs_string.replace("TDim", str(dim))
    lhs_string = lhs_string.replace("TNumNodesElem", str(nnodeselement))
    if (tensor == 1):
        lhs_string = lhs_string.replace("TTensor", "ScalarValue")
    elif (tensor == 2):
        lhs_string = lhs_string.replace("TTensor", "Vector2DValue")
    elif (tensor == 3):
        lhs_string = lhs_string.replace("TTensor", "Vector3DValue")
    lhs_string = lhs_string.replace("NumNodes", str(nnodes))
    lhs_string = lhs_string.replace("MatrixSize", str(lhs.shape[0]))
    
    rhs_string = rhs_string.replace("TDim", str(dim))
    rhs_string = rhs_string.replace("TNumNodesElem", str(nnodeselement))
    if (tensor == 1):
        rhs_string = rhs_string.replace("TTensor", "ScalarValue")
    elif (tensor == 2):
        rhs_string = rhs_string.replace("TTensor", "Vector2DValue")
    elif (tensor == 3):
        rhs_string = rhs_string.replace("TTensor", "Vector3DValue")
    rhs_string = rhs_string.replace("NumNodes", str(nnodes))
    rhs_string = rhs_string.replace("MatrixSize", str(rhs.shape[0]))

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
