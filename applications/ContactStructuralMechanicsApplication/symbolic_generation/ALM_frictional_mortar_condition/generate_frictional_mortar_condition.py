from sympy import *
from custom_sympy_fe_utilities import *
import operator 

do_simplifications = False
mode = "c" #to output to a c++ file

separate_derivatives = True
impose_partion_of_unity = False

# For debug    
#dim_combinations = [2] 
#nnodes_combinations = [2]

def convert_chain_int_int(list_slip_stick):
    value = 0
    count = 0
    for slip_stick in list_slip_stick:
        value += slip_stick * 3**count
        count += 1
    return value
        
## Debug
dim_combinations = [2]
nnodes_combinations = [2]

#dim_combinations = [2,3,3]
#nnodes_combinations = [2,3,4]

def ternary (n, size):
    nums = [0] * size
    count = 0
    while n:
        n, r = divmod(n, 3)
        nums[count] = r
        count += 1
    return nums

def itri(size):
    sol = []
    
    comb = 3**size
    
    for i in range(comb):
        sol.append(ternary(i, size))
    
    return sol

def real_norm(input):

    output = 0
    
    for i in range(input.shape[1]):
        output += input[i]**2
        
    output = real_root(output, 2)
        
    return output

lhs_string = ""
lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nbounded_matrix<double, MatrixSize, MatrixSize> AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::CalculateLocalLHS(\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const unsigned int& rMasterElementIndex,\n        const unsigned int& rActiveInactive\n        )\n{\n    bounded_matrix<double,MatrixSize,MatrixSize> lhs;\n    \n    // Master segment info\n    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();\n\n    // Initialize values\n    const bounded_matrix<double, TNumNodes, TDim> u1 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> u1old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 1);\n    const bounded_matrix<double, TNumNodes, TDim> u2 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(CurrentMasterElement, DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> u2old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(CurrentMasterElement, DISPLACEMENT, 1);\n    const bounded_matrix<double, TNumNodes, TDim> X1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(this->GetGeometry(), false);\n    const bounded_matrix<double, TNumNodes, TDim> X2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(CurrentMasterElement, false);\n    \n    const bounded_matrix<double, TNumNodes, TDim> lm = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); \n    \n    const bounded_matrix<double, TNumNodes, TDim> normalslave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  NORMAL);\n    const bounded_matrix<double, TNumNodes, TDim> tangentxislave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  TANGENT_XI);\n    const bounded_matrix<double, TNumNodes, TDim> tangentetaslave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  TANGENT_ETA);\n    \n    // Augmentation parameters\n    double scale_factor = 1.0;\n    double penalty_parameter = 0.0;\n    if (GetProperties().Has(SCALE_FACTOR) == true)\n    {\n        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);\n    }\n    if (GetProperties().Has(PENALTY_FACTOR) == true)\n    {\n        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);\n    }\n    \n    // Mortar operators\n    const bounded_matrix<double, TNumNodes, TNumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, TNumNodes, TNumNodes> DOperator = rMortarConditionMatrices.DOperator;\n    // Mortar operators derivatives\n    const array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2> DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n    const array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2> DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n\n    // We get the friction coefficient\n    const array_1d<double, TNumNodes> mu = GetFrictionCoefficient();\n\n"

lhs_template_end_string = "\n\n    return lhs;\n}\n"

rhs_string = ""
rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\narray_1d<double, MatrixSize> AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes>::CalculateLocalRHS(\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const unsigned int& rMasterElementIndex,\n        const unsigned int& rActiveInactive\n        )\n{\n    array_1d<double,MatrixSize> rhs;\n\n    // Master segment info\n    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();\n\n    // Initialize values\n    const bounded_matrix<double, TNumNodes, TDim> u1 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> u1old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 1);\n    const bounded_matrix<double, TNumNodes, TDim> u2 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(CurrentMasterElement, DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> u2old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(CurrentMasterElement, DISPLACEMENT, 1);\n    const bounded_matrix<double, TNumNodes, TDim> X1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(this->GetGeometry(), false);\n    const bounded_matrix<double, TNumNodes, TDim> X2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(CurrentMasterElement, false);\n    \n    const bounded_matrix<double, TNumNodes, TDim> lm = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); \n    \n    const bounded_matrix<double, TNumNodes, TDim> normalslave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  NORMAL);\n    const bounded_matrix<double, TNumNodes, TDim> tangentxislave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  TANGENT_XI);\n    const bounded_matrix<double, TNumNodes, TDim> tangentetaslave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  TANGENT_ETA);\n    \n    // Augmentation parameters\n    double scale_factor = 1.0;\n    double penalty_parameter = 0.0;\n    if (GetProperties().Has(SCALE_FACTOR) == true)\n    {\n        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);\n    }\n    if (GetProperties().Has(PENALTY_FACTOR) == true)\n    {\n        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);\n    }\n    \n    // Mortar operators\n    const bounded_matrix<double, TNumNodes, TNumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, TNumNodes, TNumNodes> DOperator = rMortarConditionMatrices.DOperator;\n    // We get the friction coefficient\n\n    const array_1d<double, TNumNodes> mu = GetFrictionCoefficient();\n\n"

rhs_template_end_string = "\n\n    return rhs;\n}\n"

for dim, nnodes in zip(dim_combinations, nnodes_combinations):

    slip_stick_combinations = itri(nnodes)
    
    number_dof = dim * (3 * nnodes)
        
    slip_stick_comb = 0

    for slip_stick in slip_stick_combinations: # Change output in combination of this
        
        slip_stick_comb += 1

        #Defining the unknowns
        u1 = DefineMatrix('u1',nnodes,dim) #u1(i,j) is the current displacement of node i component j at domain 1
        u2 = DefineMatrix('u2',nnodes,dim) #u2(i,j) is the current displacement of node i component j at domain 2
        u1old = DefineMatrix('u1old',nnodes,dim) #u1(i,j) is the previous displacement of node i component j at domain 1
        u2old = DefineMatrix('u2old',nnodes,dim) #u2(i,j) is the previous displacement of node i component j at domain 2
        lm = DefineMatrix('lm',nnodes,dim)
        # Normal and tangets of the slave
        normalslave = DefineMatrix('normalslave',nnodes,dim)
        tangentxislave = DefineMatrix('tangentxislave',nnodes,dim)
        tangentetaslave = DefineMatrix('tangentetaslave',nnodes,dim)
        
        # The resultant tangent
        tangentslave = DefineMatrix('tangentslave',nnodes,dim)
        
        # Define test functions
        w1 = DefineMatrix('w1',nnodes,dim)
        w2 = DefineMatrix('w2',nnodes,dim)
        wlm = DefineMatrix('wlm',nnodes,dim)
        
        # Defining normal and tangent components 
        lmnormal = DefineVector('lmnormal',nnodes)
        wlmnormal = DefineVector('wlmnormal',nnodes)
        lmtangentxi = DefineVector('lmtangentxi',nnodes)
        wlmtangentxi = DefineVector('wlmtangentxi',nnodes)
        lmtangenteta = DefineVector('lmtangenteta',nnodes)
        wlmtangenteta = DefineVector('wlmtangenteta',nnodes)
        
        # The resultant tangent LM
        lmtangent = DefineMatrix('lmtangent',nnodes,dim)
        wlmtangent = DefineMatrix('wlmtangent',nnodes,dim)
        
        for node in range(nnodes):
            lmnormal[node] = lm.row(node) * normalslave.row(node).transpose()
            wlmnormal[node] = wlm.row(node) * normalslave.row(node).transpose()
            lmtangentxi[node] = lm.row(node) * tangentxislave.row(node).transpose()
            wlmtangentxi[node] = wlm.row(node) * tangentxislave.row(node).transpose()
            lmtangenteta[node] = lm.row(node) * tangentetaslave.row(node).transpose()
            wlmtangenteta[node] = wlm.row(node) * tangentetaslave.row(node).transpose()
            
            # We calculate the LM tangent resultant
            for idim in range(dim):
                lmtangent[node,idim] = tangentxislave[node,idim] * lmtangentxi[node] + tangentetaslave[node,idim] * lmtangenteta[node]
                wlmtangent[node,idim] = tangentxislave[node,idim] * wlmtangentxi[node] + tangentetaslave[node,idim] * wlmtangenteta[node]
        
        # Now we can compute the resultant tangent
        tangentslave = lmtangent.copy()
        for node in range(nnodes):
            tangentslave /= real_norm(lmtangent.row(node))
            
        # Defining additional variables
        gap = DefineVector('gap',nnodes) 
        slipxi = DefineVector('slipxi',nnodes) 
        slipeta = DefineVector('slipeta',nnodes) 
        DOperator = DefineMatrix('DOperator',nnodes,nnodes) 
        MOperator = DefineMatrix('MOperator',nnodes,nnodes) 
        #DOperatorold = DefineMatrix('DOperatorold',nnodes,nnodes) 
        #MOperatorold = DefineMatrix('MOperatorold',nnodes,nnodes) 

        # Define other parameters
        X1 = DefineMatrix('X1',nnodes,dim)
        X2 = DefineMatrix('X2',nnodes,dim)
        x1 = X1 + u1 
        x2 = X2 + u2 
        x1old = X1 + u1old 
        x2old = X2 + u2old 

        #Define other symbols
        mu  = DefineVector('mu',nnodes) 
        penalty_parameter  = Symbol('penalty_parameter',  positive=True)
        scale_factor = Symbol('scale_factor', positive=True)

        # Define variables list for later derivation
        u1_var = []
        u2_var = []
        lm_var = []
        CreateVariableMatrixList(u1_var, u1)
        CreateVariableMatrixList(u2_var, u2)
        u12_var=u1_var.copy()
        u1_lm_var=u1_var.copy()
        CreateVariableMatrixList(u12_var, u2)
        CreateVariableMatrixList(lm_var, lm)
        CreateVariableMatrixList(u1_lm_var, lm)
        all_var=u12_var.copy()
        CreateVariableMatrixList(all_var, lm)

        # Force the variables to be dependendant of the DOF
        #normalslave = DefineDofDependencyMatrix(normalslave, u1_var)
        #tangentxislave = DefineDofDependencyMatrix(tangentxislave, u1_var)
        #tangentetaslave = DefineDofDependencyMatrix(tangentetaslave, u1_var)
        DOperator = DefineDofDependencyMatrix(DOperator, u12_var) # If you consider Gitterle you need to keep the old operators
        MOperator = DefineDofDependencyMatrix(MOperator, u12_var)

        # Defining the normal gap
        Dx1Mx2 = DOperator * x1 - MOperator * x2
        Dx1oldMx2old = DOperator * x1old - MOperator * x2old
        Dw1Mw2 = DOperator * w1 - MOperator * w2
        for node in range(nnodes): 
            gap[node]  = Dx1Mx2.row(node) * - normalslave.row(node).transpose()
            slipxi[node]  = Dx1oldMx2old.row(node) * tangentxislave.row(node).transpose()
            slipeta[node]  = Dx1oldMx2old.row(node) * tangentetaslave.row(node).transpose()

        # Define dofs & test function vector
        dofs = Matrix( zeros(number_dof, 1) )
        testfunc = Matrix( zeros(number_dof, 1) )
        count = 0
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

        # Compute galerkin functional # NOTE: Maybe you can define a different penalty and scale factor in the tangent direction
        rv_galerkin = 0
        for node in range(nnodes):
            slip = slip_stick[node] 
            if (slip == 0):  
                rv_galerkin -= 0.5/penalty_parameter * scale_factor**2.0 * lmnormal[node] * wlmnormal[node]
                rv_galerkin -= 0.5/penalty_parameter * scale_factor**2.0 * lmtangentxi[node] * wlmtangentxi[node]
                rv_galerkin -= 0.5/penalty_parameter * scale_factor**2.0 * lmtangenteta[node] * wlmtangenteta[node]
            else:
                rv_galerkin += ((((scale_factor * lmnormal[node] + penalty_parameter * gap[node]) * normalslave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                rv_galerkin +=  scale_factor * gap[node] * wlmnormal[node]
                
                if (slip == 1): # Slip 
                    rv_galerkin -= (((mu[node] * (scale_factor * lmnormal[node] + penalty_parameter * gap[node]) * tangentslave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                    rv_galerkin -=  (0.5/penalty_parameter * scale_factor * (scale_factor * lmtangent.row(node) + mu[node] * (scale_factor * lmnormal[node] + penalty_parameter * gap[node]) * tangentslave.row(node)) * wlmtangent.row(node).transpose())[0,0] 
                else: # Stick 
                    rv_galerkin += ((((scale_factor * lmtangentxi[node] + penalty_parameter * slipxi[node]) * tangentxislave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                    rv_galerkin +=  scale_factor * slipxi[node] * wlmtangentxi[node]
                    rv_galerkin += ((((scale_factor * lmtangenteta[node] + penalty_parameter * slipeta[node]) * tangentetaslave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                    rv_galerkin +=  scale_factor * slipeta[node] * wlmtangenteta[node]

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

        if (slip_stick_comb == 1):
            lhs_string += lhs_template_begin_string
            lhs_string += "    if (rActiveInactive == " + str(convert_chain_int_int(slip_stick)) + " )\n    {\n    "
        else:
            lhs_string += "    else if (rActiveInactive == "+str(convert_chain_int_int(slip_stick)) + " )\n    {\n    "
        lhs_string += lhs_out.replace("\n","\n    ")
        lhs_string += "}\n"    
        
        if (slip_stick_comb == len(slip_stick_combinations)):
            lhs_string += lhs_template_end_string
    
        if (slip_stick_comb == 1):
            rhs_string += rhs_template_begin_string
            rhs_string += "    if (rActiveInactive == "+str(convert_chain_int_int(slip_stick)) + " )\n    {\n    "
        else:
            rhs_string += "    else if (rActiveInactive == "+str(convert_chain_int_int(slip_stick)) + " )\n    {\n    "
        rhs_string += rhs_out.replace("\n","\n    ")
        rhs_string += "}\n"
        
        if (slip_stick_comb == len(slip_stick_combinations)):
            rhs_string += rhs_template_end_string
            
    lhs_string = lhs_string.replace("TDim", str(dim))
    lhs_string = lhs_string.replace("TNumNodes", str(nnodes))
    lhs_string = lhs_string.replace("MatrixSize", str(lhs.shape[0]))
    lhs_string = lhs_string.replace("SIZEDERIVATIVES2", str(2 * (nnodes * dim)))

    rhs_string = rhs_string.replace("TDim", str(dim))
    rhs_string = rhs_string.replace("TNumNodes", str(nnodes))
    rhs_string = rhs_string.replace("MatrixSize", str(lhs.shape[0]))

    ##############################################################################
    ##############################################################################
    ##################### DEFINE VARIABLES AND DERIVATIVES #######################
    ##############################################################################
    ##############################################################################

    var_strings = []
    var_strings_subs = []
    var_strings_aux_subs = []
    der_var_strings = []
    der_var_list = []
    der_var_used_index = []

    #var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(normalslave, "normalslave", "normalslave", u1_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
    #var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(tangentxislave, "tangentxislave", "tangentxislave", u1_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
    #var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(tangentetaslave, "tangentetaslave", "tangentetaslave", u1_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
    var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(DOperator, "DOperator", "DOperator", u12_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
    var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(MOperator, "MOperator", "MOperator", u12_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")

    #############################################################################
    ############################### SUBSTITUTION ################################
    #############################################################################

    # Replace all
    lhs_string = lhs_string.replace("//subsvar_", "")
    rhs_string = rhs_string.replace("//subsvar_", "")

    for index in range(len(der_var_strings)):
        if (separate_derivatives == True):
            lhs_string = lhs_string.replace(der_var_strings[index], der_var_list[index])
            rhs_string = rhs_string.replace(der_var_strings[index], der_var_list[index])
        else:
            aux_out = ccode(der_var_strings_subs[index])
            if ("// Not supported in C:") in aux_out:
                print("WARNING")
                print(der_var_strings[index])
                print(der_var_strings_subs[index])
            lhs_string = lhs_string.replace(der_var_strings[index], aux_out)
            rhs_string = rhs_string.replace(der_var_strings[index], aux_out)
            
    for index in range(len(var_strings)):
        lhs_string = lhs_string.replace(var_strings[index], var_strings_subs[index])
        rhs_string = rhs_string.replace(var_strings[index], var_strings_subs[index])

    lhs_string = SubstituteIndex(lhs_string, mode, number_dof)
    rhs_string = SubstituteIndex(rhs_string, mode, number_dof)
    lhs_string = lhs_string.replace("array[1]d", "array_1d") # Repair the definition
    rhs_string = rhs_string.replace("array[1]d", "array_1d") # Repair the definition

#############################################################################
############################### SIMPLIFICATION ##############################
#############################################################################

lhs_string = lhs_string.replace("pow(", "std::pow(")
lhs_string = lhs_string.replace("sqrt(", "std::sqrt(")
rhs_string = rhs_string.replace("pow(", "std::pow(")
rhs_string = rhs_string.replace("sqrt(", "std::sqrt(")

#for dof in reversed(range(len(u1_var))):
    #lhs_string = lhs_string.replace("Deltanormalslave"+str(dof), "Deltanormalslave["+str(dof)+"]")
for dof in reversed(range(len(u12_var))):
    lhs_string = lhs_string.replace("DeltaDOperator"+str(dof), "DeltaDOperator["+str(dof)+"]")
for dof in reversed(range(len(u12_var))):
    lhs_string = lhs_string.replace("DeltaMOperator"+str(dof), "DeltaMOperator["+str(dof)+"]")
    
#############################################################################
################################# FINAL SAVING ##############################
#############################################################################

input = open("ALM_frictional_mortar_contact_condition_template.cpp",'r').read()
outputstring = input.replace("// replace_lhs", lhs_string)
outputstring = outputstring.replace("// replace_rhs", rhs_string)
output = open("ALM_frictional_mortar_contact_condition.cpp",'w')
output.write(outputstring)
output.close()

print("Strings have been replaced...")

print("Process Finished..................")
