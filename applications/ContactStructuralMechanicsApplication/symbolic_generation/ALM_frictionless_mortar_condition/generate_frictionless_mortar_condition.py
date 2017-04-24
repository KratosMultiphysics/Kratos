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

def convert_active_inactive_int(list_active):
    value = 0
    count = 0
    for active_inactive in list_active:
        value += active_inactive * 2**count
        count += 1
    return value
        
## Debug
#dim_combinations = [2]
#nnodes_combinations = [2]

dim_combinations = [2,3,3]
nnodes_combinations = [2,3,4]

lhs_string = ""
lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\ntemplate<>\nbounded_matrix<double, TMatrixSize, TMatrixSize> AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateLocalLHS<TMatrixSize>(\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const unsigned int& rMasterElementIndex,\n        const unsigned int& rActiveInactive\n        )\n{\n    bounded_matrix<double,TMatrixSize,TMatrixSize> lhs;\n    \n    // Master segment info\n    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();\n\n    // Initialize values\n    const bounded_matrix<double, TNumNodes, TDim> u1 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> u2 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(CurrentMasterElement, DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> X1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(this->GetGeometry(), false);\n    const bounded_matrix<double, TNumNodes, TDim> X2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(CurrentMasterElement, false);\n    \n    const array_1d<double, TNumNodes> lmnormal = ContactUtilities::GetVariableVector<TNumNodes>(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); \n    \n    const bounded_matrix<double, TNumNodes, TDim> normalslave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  NORMAL); \n    \n    // Augmentation parameters\n    double scale_factor = 1.0;\n    double penalty_parameter = 0.0;\n    if (GetProperties().Has(SCALE_FACTOR) == true)\n    {\n        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);\n    }\n    if (GetProperties().Has(PENALTY_FACTOR) == true)\n    {\n        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);\n    }\n    \n    // Mortar operators\n    const bounded_matrix<double, TNumNodes, TNumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, TNumNodes, TNumNodes> DOperator = rMortarConditionMatrices.DOperator;\n    // Mortar operators derivatives\n    const array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2> DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n    const array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2> DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n\n"

lhs_template_end_string = "\n\n    return lhs;\n}\n"

rhs_string = ""
rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\ntemplate<>\narray_1d<double, TMatrixSize> AugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes>::CalculateLocalRHS<TMatrixSize>(\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const unsigned int& rMasterElementIndex,\n        const unsigned int& rActiveInactive\n        )\n{\n    array_1d<double,TMatrixSize> rhs;\n\n    // Master segment info\n    GeometryType& CurrentMasterElement = mThisMasterElements[rMasterElementIndex]->GetGeometry();\n\n    // Initialize values\n    const bounded_matrix<double, TNumNodes, TDim> u1 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> u2 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(CurrentMasterElement, DISPLACEMENT, 0);\n    const bounded_matrix<double, TNumNodes, TDim> X1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(this->GetGeometry(), false);\n    const bounded_matrix<double, TNumNodes, TDim> X2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(CurrentMasterElement, false);\n    \n    const array_1d<double, TNumNodes> lmnormal = ContactUtilities::GetVariableVector<TNumNodes>(this->GetGeometry(), NORMAL_CONTACT_STRESS, 0); \n    \n    const bounded_matrix<double, TNumNodes, TDim> normalslave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(),  NORMAL); \n    \n    // Augmentation parameters\n    double scale_factor = 1.0;\n    double penalty_parameter = 0.0;\n    if (GetProperties().Has(SCALE_FACTOR) == true)\n    {\n        scale_factor  = GetProperties().GetValue(SCALE_FACTOR);\n    }\n    if (GetProperties().Has(PENALTY_FACTOR) == true)\n    {\n        penalty_parameter = GetProperties().GetValue(PENALTY_FACTOR);\n    }\n    \n    // Mortar operators\n    const bounded_matrix<double, TNumNodes, TNumNodes> MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, TNumNodes, TNumNodes> DOperator = rMortarConditionMatrices.DOperator;\n\n"

rhs_template_end_string = "\n\n    return rhs;\n}\n"

for dim, nnodes in zip(dim_combinations, nnodes_combinations):

    from sympy.utilities.iterables import ibin
    active_inactive_combinations = list(ibin(nnodes, 'all'))
    
    active_inactive_comb = 0
    
    for active_inactive in active_inactive_combinations: # Change output in combination of this
        
        active_inactive_comb += 1
        number_dof = dim * (2 * nnodes) + nnodes

        #Defining the unknowns
        u1 = DefineMatrix('u1',nnodes,dim) #u1(i,j) is displacement of node i component j at domain 1
        u2 = DefineMatrix('u2',nnodes,dim) #u2(i,j) is displacement of node i component j at domain 2
        lmnormal = DefineVector('lmnormal',nnodes) 
        gap = DefineVector('gap',nnodes) 
        DOperator = DefineMatrix('DOperator',nnodes,nnodes) 
        MOperator = DefineMatrix('MOperator',nnodes,nnodes) 

        # Define other parameters
        # Normal and tangets of the slave
        normalslave = DefineMatrix('normalslave',nnodes,dim)
            
        X1 = DefineMatrix('X1',nnodes,dim)
        X2 = DefineMatrix('X2',nnodes,dim)
        x1 = X1 + u1 
        x2 = X2 + u2 

        #Define other symbols
        penalty_parameter  = Symbol('penalty_parameter',  positive=True)
        scale_factor = Symbol('scale_factor', positive=True)

        # Define test functions
        w1 = DefineMatrix('w1',nnodes,dim)
        w2 = DefineMatrix('w2',nnodes,dim)
        wlmnormal = DefineVector('wlmnormal',nnodes)

        # Define variables list for later derivation
        u1_var = []
        u2_var = []
        lm_var = []
        CreateVariableMatrixList(u1_var, u1)
        CreateVariableMatrixList(u2_var, u2)
        u12_var=u1_var.copy()
        u1_lm_var=u1_var.copy()
        CreateVariableMatrixList(u12_var, u2)
        CreateVariableVectorList(lm_var, lmnormal)
        CreateVariableVectorList(u1_lm_var, lmnormal)
        all_var=u12_var.copy()
        CreateVariableVectorList(all_var, lmnormal)

        # Force the variables to be dependendant of the DOF
        #normalslave = DefineDofDependencyMatrix(normalslave, u1_var)
        DOperator = DefineDofDependencyMatrix(DOperator, u12_var)
        MOperator = DefineDofDependencyMatrix(MOperator, u12_var)

        # Defining the normal gap
        Dx1Mx2 = DOperator * x1 - MOperator * x2
        Dw1Mw2 = DOperator * w1 - MOperator * w2
        for node in range(nnodes): 
            gap[node]  = Dx1Mx2.row(node) * - normalslave.row(node).transpose() # NOTE: SIGN?¿

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
            dofs[count] = lmnormal[i]
            testfunc[count] = wlmnormal[i]
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
            active = active_inactive[node] 
            if (active == 1):                
                rv_galerkin += ((((scale_factor * lmnormal[node] + penalty_parameter * gap[node]) * normalslave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                rv_galerkin +=  scale_factor * gap[node] * wlmnormal[node]
            else:
                rv_galerkin += - 0.5/penalty_parameter * scale_factor**2.0 * lmnormal[node] * wlmnormal[node]

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

        if (active_inactive_comb == 1):
            lhs_string += lhs_template_begin_string
            lhs_string += "    if (rActiveInactive == " + str(convert_active_inactive_int(active_inactive)) + " )\n    {\n    "
        else:
            lhs_string += "    else if (rActiveInactive == "+str(convert_active_inactive_int(active_inactive)) + " )\n    {\n    "
        lhs_string += lhs_out.replace("\n","\n    ")
        lhs_string += "}\n"    
        
        if (active_inactive_comb == len(active_inactive_combinations)):
            lhs_string += lhs_template_end_string
    
        if (active_inactive_comb == 1):
            rhs_string += rhs_template_begin_string
            rhs_string += "    if (rActiveInactive == "+str(convert_active_inactive_int(active_inactive)) + " )\n    {\n    "
        else:
            rhs_string += "    else if (rActiveInactive == "+str(convert_active_inactive_int(active_inactive)) + " )\n    {\n    "
        rhs_string += rhs_out.replace("\n","\n    ")
        rhs_string += "}\n"
        
        if (active_inactive_comb == len(active_inactive_combinations)):
            rhs_string += rhs_template_end_string
        
        lhs_string = lhs_string.replace("TDim", str(dim))
        lhs_string = lhs_string.replace("TNumNodes", str(nnodes))
        lhs_string = lhs_string.replace("TMatrixSize", str(lhs.shape[0]))
        lhs_string = lhs_string.replace("SIZEDERIVATIVES2", str(2 * (nnodes * dim)))

        rhs_string = rhs_string.replace("TDim", str(dim))
        rhs_string = rhs_string.replace("TNumNodes", str(nnodes))
        rhs_string = rhs_string.replace("TMatrixSize", str(lhs.shape[0]))
    
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

out = open("lhs.txt",'w')
out.write(lhs_string)
out.close()

out = open("rhs.txt",'w')
out.write(rhs_string)
out.close()

print("Strings have been replaced...")

print("Process Finished..................")
