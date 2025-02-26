

# Import KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.ContactStructuralMechanicsApplication

# Import sympy utils
import sympy
from KratosMultiphysics.ContactStructuralMechanicsApplication import custom_sympy_fe_utilities

do_simplifications = False
mode = "c" #to output to a c++ file

impose_partion_of_unity = False

def convert_active_inactive_int(list_active):
    value = 0
    count = 0
    for active_inactive in list_active:
        value += active_inactive * 2**count
        count += 1
    return value

# Debug
debug = False
#debug = True # NOTE: COMMENT FOR NOT DEBUG
debug_counter = 0

dim_combinations = [2,3,3,3,3]
nnodes_combinations = [2,3,4,3,4]
nnodes_master_combinations = [2,3,4,4,3]
normal_combs = 2

lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nvoid AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::CalculateLocalLHS(\n    Matrix& rLocalLHS,\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    const DerivativeDataType& rDerivativeData,\n    const IndexType rActiveInactive,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    // Initialize\n    for (std::size_t i = 0; i < MatrixSize; ++i)\n        for (std::size_t j = 0; j < MatrixSize; ++j)\n            rLocalLHS(i, j) = 0.0;\n\n    // The geometry of the condition\n    const GeometryType& r_geometry = this->GetParentGeometry();\n\n    // Initialize values\n    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n    \n    const BoundedMatrix<double, TNumNodes, TDim> LM = MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(this->GetParentGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0);\n    \n    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n\n    // The ALM parameters\n    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(this->GetParentGeometry(), DYNAMIC_FACTOR);\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    \n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n    // Mortar operators derivatives\n    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodesMaster>, SIZEDERIVATIVES2>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n\n"

lhs_template_end_string = "\n}\n"

rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nvoid AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::CalculateLocalRHS(\n    Vector& rLocalRHS,\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    const DerivativeDataType& rDerivativeData,\n    const IndexType rActiveInactive,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    // Initialize\n    for (std::size_t i = 0; i < MatrixSize; ++i)\n        rLocalRHS[i] = 0.0;\n\n    // The geometry of the condition\n    const GeometryType& r_geometry = this->GetParentGeometry();\n\n    // Initialize values\n    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n    \n    const BoundedMatrix<double, TNumNodes, TDim> LM = MortarUtilities::GetVariableMatrix<TDim, TNumNodes>(this->GetParentGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0);\n    \n    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n\n    // The ALM parameters\n    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(this->GetParentGeometry(), DYNAMIC_FACTOR);\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    \n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n\n"

rhs_template_end_string = "\n}\n"

# We track the output each time
output_count = 0
total_combs = normal_combs * len(nnodes_combinations)

for normalvar in range(normal_combs):

    if normalvar == 0:
        normalvarstring = "false"
    else:
        normalvarstring = "true"

    if normalvar == 1:
        lhs_template_begin_string += "    const array_1d<BoundedMatrix<double, TNumNodes, TDim>, SIZEDERIVATIVES1>& DeltaNormalSlave = rDerivativeData.DeltaNormalSlave;\n\n"

    for dim, nnodes, nnodes_master in zip(dim_combinations, nnodes_combinations, nnodes_master_combinations):

        # Update counter and clear strings
        lhs_string = ""
        rhs_string = ""
        output_count += 1
        number_dof = dim * (nnodes_master + 2 * nnodes)

        #Defining the unknowns
        u1 = custom_sympy_fe_utilities.DefineMatrix('u1',nnodes,dim)        #u1(i,j) is displacement of node i component j at domain 1
        u2 = custom_sympy_fe_utilities.DefineMatrix('u2',nnodes_master,dim) #u2(i,j) is displacement of node i component j at domain 2
        LM = custom_sympy_fe_utilities.DefineMatrix('LM',nnodes,dim)
        NormalGap = custom_sympy_fe_utilities.DefineVector('NormalGap',nnodes)
        DOperator = custom_sympy_fe_utilities.DefineMatrix('DOperator',nnodes,nnodes)
        MOperator = custom_sympy_fe_utilities.DefineMatrix('MOperator',nnodes,nnodes_master)

        # Define other parameters
        # Normal and tangents of the slave
        NormalSlave = custom_sympy_fe_utilities.DefineMatrix('NormalSlave',nnodes,dim)

        X1 = custom_sympy_fe_utilities.DefineMatrix('X1',nnodes,dim)
        X2 = custom_sympy_fe_utilities.DefineMatrix('X2',nnodes_master,dim)
        x1 = X1 + u1
        x2 = X2 + u2

        #Define other symbols
        DynamicFactor  = custom_sympy_fe_utilities.DefineVector('DynamicFactor',nnodes)
        PenaltyParameter  = custom_sympy_fe_utilities.DefineVector('PenaltyParameter',nnodes)
        ScaleFactor = sympy.Symbol('ScaleFactor', positive=True)

        # Define test functions
        w1 = custom_sympy_fe_utilities.DefineMatrix('w1',nnodes,dim)
        w2 = custom_sympy_fe_utilities.DefineMatrix('w2',nnodes_master,dim)
        wLM = custom_sympy_fe_utilities.DefineMatrix('wLMN',nnodes,dim)

        # Defining normal and tangent components
        LMNormal = custom_sympy_fe_utilities.DefineVector('LMNormal',nnodes)
        wLMNormal = custom_sympy_fe_utilities.DefineVector('wLMNormal',nnodes)

        # The resultant tangent LM
        LMTangent = custom_sympy_fe_utilities.DefineMatrix('LMTangent',nnodes,dim)
        wLMTangent = custom_sympy_fe_utilities.DefineMatrix('wLMTangent',nnodes,dim)

        for node in range(nnodes):
            LMNormal[node] = LM.row(node).dot(NormalSlave.row(node))
            wLMNormal[node] = wLM.row(node).dot(NormalSlave.row(node))

            # We calculate the LM tangent resultant
            for idim in range(dim):
                LMTangent[node,idim] = LM[node,idim] - LMNormal[node] * NormalSlave[node,idim]
                wLMTangent[node,idim] = wLM[node,idim] - wLMNormal[node] * NormalSlave[node,idim]

        # Define variables list for later derivation
        u1_var = []
        u2_var = []
        lm_var = []
        custom_sympy_fe_utilities.CreateVariableMatrixList(u1_var, u1)
        custom_sympy_fe_utilities.CreateVariableMatrixList(u2_var, u2)
        u12_var=u1_var.copy()
        u1_lm_var=u1_var.copy()
        custom_sympy_fe_utilities.CreateVariableMatrixList(u12_var, u2)
        custom_sympy_fe_utilities.CreateVariableMatrixList(lm_var, LM)
        custom_sympy_fe_utilities.CreateVariableMatrixList(u1_lm_var, LM)
        all_var=u12_var.copy()
        custom_sympy_fe_utilities.CreateVariableMatrixList(all_var, LM)

        # Force the variables to be dependendant of the DOF
        if normalvar == 1:
            NormalSlave = custom_sympy_fe_utilities.DefineDofDependencyMatrix(NormalSlave, u1_var)
        DOperator = custom_sympy_fe_utilities.DefineDofDependencyMatrix(DOperator, u12_var)
        MOperator = custom_sympy_fe_utilities.DefineDofDependencyMatrix(MOperator, u12_var)

        # Defining the normal NormalGap
        Dx1Mx2 = DOperator * x1 - MOperator * x2
        Dw1Mw2 = DOperator * w1 - MOperator * w2

        for node in range(nnodes):
            NormalGap[node] = - Dx1Mx2.row(node).dot(NormalSlave.row(node))

        # Define dofs & test function vector
        dofs = sympy.Matrix( sympy.zeros(number_dof, 1) )
        testfunc = sympy.Matrix( sympy.zeros(number_dof, 1) )
        count = 0
        for i in range(0,nnodes_master):
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
                dofs[count] = LM[i,k]
                testfunc[count] = wLM[i,k]
                count+=1
        print("dofs = ",dofs)
        print("testfunc = ",testfunc)

        #############################################################################
        #############################################################################
        ########################## FUNCTIONAL DEFINITION ############################
        #############################################################################
        #############################################################################

        # Compute galerkin functional
        lhs_string += lhs_template_begin_string
        rhs_string += rhs_template_begin_string
        if debug_counter == 0:
            for node in range(nnodes):
                for active in range(2):
                    rv_galerkin = 0
                    if active == 1: # Active
                        augmented_lm = (ScaleFactor * LM.row(node) + PenaltyParameter[node] * NormalGap[node] * NormalSlave.row(node))
                        rv_galerkin += DynamicFactor[node] * (augmented_lm).dot(Dw1Mw2.row(node))
                        rv_galerkin += ScaleFactor * NormalGap[node] * wLMNormal[node]
                        rv_galerkin -= ScaleFactor**2/PenaltyParameter[node] * (wLMTangent.row(node).dot(LMTangent.row(node)))
                    else: # Inactive
                        rv_galerkin -= ScaleFactor**2/PenaltyParameter[node] * (wLM.row(node).dot(LM.row(node)))

                    if do_simplifications:
                        rv_galerkin = sympy.simplify(rv_galerkin)

                    #############################################################################
                    # Complete functional
                    rv = sympy.Matrix( sympy.zeros(1, 1) )
                    rv[0,0] = rv_galerkin

                    rhs,lhs = custom_sympy_fe_utilities.Compute_RHS_and_LHS(rv.copy(), testfunc, dofs, False)
                    print("LHS= ",lhs.shape)
                    print("RHS= ",rhs.shape)
                    print("LHS and RHS have been created!")

                    lhs_out = custom_sympy_fe_utilities.OutputMatrix_CollectingFactorsNonZero(lhs,"lhs", mode, 1, number_dof)
                    rhs_out = custom_sympy_fe_utilities.OutputVector_CollectingFactorsNonZero(rhs,"rhs", mode, 1, number_dof)
                    print("Substitution strings are ready....")

                    if active == 0:
                        lhs_string += "    \n    // NODE " + str(node) + "\n"
                        lhs_string += "    if (r_geometry[" + str(node) + "].IsNot(ACTIVE)) { // INACTIVE\n    "
                    else:
                        lhs_string += "} else { // ACTIVE\n    "
                    lhs_string += lhs_out.replace("\n","\n    ")
                    if active == 1:
                        lhs_string += "}"

                    if active == 0:
                        rhs_string += "    \n    // NODE " + str(node) + "\n"
                        rhs_string += "    if (r_geometry[" + str(node) + "].IsNot(ACTIVE)) { // INACTIVE\n    "
                    else:
                        rhs_string += "} else { // ACTIVE\n    "
                    rhs_string += rhs_out.replace("\n","\n    ")
                    if active == 1:
                        rhs_string += "}"

        lhs_string += lhs_template_end_string
        rhs_string += rhs_template_end_string
        if debug:
            debug_counter += 1

        lhs_string = lhs_string.replace("TDim", str(dim))
        lhs_string = lhs_string.replace("TNumNodesMaster", str(nnodes_master))
        lhs_string = lhs_string.replace("TNumNodes", str(nnodes))
        lhs_string = lhs_string.replace("MatrixSize", str(lhs.shape[0]))
        lhs_string = lhs_string.replace("TNormalVariation", normalvarstring)
        lhs_string = lhs_string.replace("SIZEDERIVATIVES1", str(((nnodes) * dim)))
        lhs_string = lhs_string.replace("SIZEDERIVATIVES2", str(((nnodes + nnodes_master) * dim)))

        rhs_string = rhs_string.replace("TDim", str(dim))
        rhs_string = rhs_string.replace("TNumNodesMaster", str(nnodes_master))
        rhs_string = rhs_string.replace("TNumNodes", str(nnodes))
        rhs_string = rhs_string.replace("TNormalVariation", normalvarstring)
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

        if normalvar == 1:
            var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = custom_sympy_fe_utilities.DefineVariableLists(NormalSlave, "NormalSlave", "NormalSlave", u1_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
        var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = custom_sympy_fe_utilities.DefineVariableLists(DOperator, "DOperator", "DOperator", u12_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
        var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = custom_sympy_fe_utilities.DefineVariableLists(MOperator, "MOperator", "MOperator", u12_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")

        #############################################################################
        ############################### SUBSTITUTION ################################
        #############################################################################

        # Replace all
        lhs_string = lhs_string.replace("//subsvar_", "")
        rhs_string = rhs_string.replace("//subsvar_", "")

        for index in range(len(der_var_strings)):
            lhs_string = lhs_string.replace(der_var_strings[index], der_var_list[index])
            rhs_string = rhs_string.replace(der_var_strings[index], der_var_list[index])

        for index in range(len(var_strings)):
            lhs_string = lhs_string.replace(var_strings[index], var_strings_subs[index])
            rhs_string = rhs_string.replace(var_strings[index], var_strings_subs[index])

        lhs_string = custom_sympy_fe_utilities.SubstituteIndex(lhs_string, mode, number_dof)
        rhs_string = custom_sympy_fe_utilities.SubstituteIndex(rhs_string, mode, number_dof)
        lhs_string = lhs_string.replace("array[1]d", "array_1d") # Repair the definition
        rhs_string = rhs_string.replace("array[1]d", "array_1d") # Repair the definition

        #############################################################################
        ############################### SIMPLIFICATION ##############################
        #############################################################################

        lhs_string = lhs_string.replace("pow(", "std::pow(")
        lhs_string = lhs_string.replace("sqrt(", "std::sqrt(")
        rhs_string = rhs_string.replace("pow(", "std::pow(")
        rhs_string = rhs_string.replace("sqrt(", "std::sqrt(")
        lhs_string = lhs_string.replace("lhs(", "rLocalLHS(")
        rhs_string = rhs_string.replace("rhs[", "rLocalRHS[")

        for dof in reversed(range(len(u1_var))):
            lhs_string = lhs_string.replace("DeltaNormalSlave"+str(dof), "DeltaNormalSlave["+str(dof)+"]")
        for dof in reversed(range(len(u12_var))):
            lhs_string = lhs_string.replace("DeltaDOperator"+str(dof), "DeltaDOperator["+str(dof)+"]")
        for dof in reversed(range(len(u12_var))):
            lhs_string = lhs_string.replace("DeltaMOperator"+str(dof), "DeltaMOperator["+str(dof)+"]")

        #############################################################################
        ################################# FINAL SAVING ##############################
        #############################################################################

        if output_count == 1:
            first_input = open("ALM_frictionless_components_mortar_contact_condition_template.cpp", 'r').read()
            lhs_string += "// replace_lhs"
            outputstring = first_input.replace("// replace_lhs", lhs_string)
        else:
            input = open("ALM_frictionless_components_mortar_contact_condition.cpp", 'r').read()
            if output_count < total_combs:
                lhs_string += "// replace_lhs"
            outputstring = input.replace("// replace_lhs", lhs_string)

        if output_count < total_combs:
            rhs_string += "// replace_rhs"

        outputstring = outputstring.replace("// replace_rhs", rhs_string)
        output = open("ALM_frictionless_components_mortar_contact_condition.cpp",'w')
        output.write(outputstring)
        output.close()

print("Strings have been replaced...")

print("Process Finished..................")
