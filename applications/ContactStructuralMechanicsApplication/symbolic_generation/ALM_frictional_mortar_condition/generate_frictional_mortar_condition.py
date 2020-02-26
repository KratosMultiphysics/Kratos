from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.ContactStructuralMechanicsApplication

# Import sympy utils
import sympy
from KratosMultiphysics.ContactStructuralMechanicsApplication import custom_sympy_fe_utilities

do_simplifications = False
mode = "c" #to output to a c++ file

impose_partion_of_unity = False

def convert_chain_int_int(list_slip_stick):
    value = 0
    count = 0
    for slip_stick in list_slip_stick:
        value += slip_stick * 3**count
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

def real_norm(input):

    output = 0

    for i in range(input.shape[1]):
        output += input[i]**2

    output = sympy.real_root(output, 2)

    return output

lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nvoid AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::CalculateLocalLHS(\n    Matrix& rLocalLHS,\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    const DerivativeDataType& rDerivativeData,\n    const IndexType rActiveInactive,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    // Initialize\n    for (std::size_t i = 0; i < MatrixSize; ++i)\n        for (std::size_t j = 0; j < MatrixSize; ++j)\n            rLocalLHS(i, j) = 0.0;\n\n    // The geometry of the condition\n    const GeometryType& r_geometry = this->GetParentGeometry();\n\n    // Initialize values\n    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    const BoundedMatrix<double, TNumNodes, TDim>& u1old = rDerivativeData.u1old;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2old = rDerivativeData.u2old;\n    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n    \n    const BoundedMatrix<double, TNumNodes, TDim> LM = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_geometry, VECTOR_LAGRANGE_MULTIPLIER, 0);\n    \n    // The normal and tangent vectors\n    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n    const BoundedMatrix<double, TNumNodes, TDim> TangentSlave = MortarUtilities::ComputeTangentMatrix<TNumNodes,TDim>(r_geometry);\n\n    // The ALM parameters\n    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    const double TangentFactor = rDerivativeData.TangentFactor;\n    \n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& MOperatorold = mPreviousMortarOperators.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperatorold = mPreviousMortarOperators.DOperator;\n\n    // Mortar operators derivatives\n    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n\n    // We get the friction coefficient\n    const array_1d<double, TNumNodes> mu = GetFrictionCoefficient();\n\n//    // The delta time\n//    const double delta_time = rCurrentProcessInfo[DELTA_TIME];\n\n    const double OperatorThreshold = rCurrentProcessInfo[OPERATOR_THRESHOLD];\n    const double norm_delta_M = norm_frobenius(MOperator - MOperatorold);\n    const double norm_delta_D = norm_frobenius(DOperator - DOperatorold);\n    const bool is_objetive = (norm_delta_D > OperatorThreshold && norm_delta_M > OperatorThreshold) ? true : false;\n    this->Set(MODIFIED, !is_objetive);\n"

lhs_template_end_string = "}\n"

rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nvoid AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes, TNormalVariation, TNumNodesMaster>::CalculateLocalRHS(\n    Vector& rLocalRHS,\n    const MortarConditionMatrices& rMortarConditionMatrices,\n    const DerivativeDataType& rDerivativeData,\n    const IndexType rActiveInactive,\n    const ProcessInfo& rCurrentProcessInfo\n    )\n{\n    // Initialize\n    for (std::size_t i = 0; i < MatrixSize; ++i)\n        rLocalRHS[i] = 0.0;\n\n    // The geometry of the condition\n    const GeometryType& r_geometry = this->GetParentGeometry();\n\n    // Initialize values\n    const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    const BoundedMatrix<double, TNumNodes, TDim>& u1old = rDerivativeData.u1old;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2 = rDerivativeData.u2;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& u2old = rDerivativeData.u2old;\n    const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    const BoundedMatrix<double, TNumNodesMaster, TDim>& X2 = rDerivativeData.X2;\n    \n    const BoundedMatrix<double, TNumNodes, TDim> LM = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(r_geometry, VECTOR_LAGRANGE_MULTIPLIER, 0);\n    \n    // The normal and tangent vectors\n    const BoundedMatrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n    const BoundedMatrix<double, TNumNodes, TDim> TangentSlave = MortarUtilities::ComputeTangentMatrix<TNumNodes,TDim>(r_geometry);\n\n    // The ALM parameters\n    const array_1d<double, TNumNodes> DynamicFactor = MortarUtilities::GetVariableVector<TNumNodes>(r_geometry, DYNAMIC_FACTOR);\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    const double TangentFactor = rDerivativeData.TangentFactor;\n    \n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& MOperatorold = mPreviousMortarOperators.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperatorold = mPreviousMortarOperators.DOperator;\n\n    // We get the friction coefficient\n    const array_1d<double, TNumNodes> mu = GetFrictionCoefficient();\n\n//    // The delta time\n//    const double delta_time = rCurrentProcessInfo[DELTA_TIME];\n\n    const double OperatorThreshold = rCurrentProcessInfo[OPERATOR_THRESHOLD];\n    const double norm_delta_M = norm_frobenius(MOperator - MOperatorold);\n    const double norm_delta_D = norm_frobenius(DOperator - DOperatorold);\n    const bool is_objetive = (norm_delta_D > OperatorThreshold && norm_delta_M > OperatorThreshold) ? true : false;\n    this->Set(MODIFIED, !is_objetive);\n"

rhs_template_end_string = "}\n"

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

        # Update counter and calculate dof
        lhs_string = ""
        rhs_string = ""
        output_count += 1
        number_dof = dim * (2 * nnodes + nnodes_master)

        # Defining the unknowns
        u1 = custom_sympy_fe_utilities.DefineMatrix('u1', nnodes, dim, "Symbol")              # u1(i,j) is the current displacement of node i component j at domain 1
        u2 = custom_sympy_fe_utilities.DefineMatrix('u2', nnodes_master, dim, "Symbol")       # u2(i,j) is the current displacement of node i component j at domain 2
        u1old = custom_sympy_fe_utilities.DefineMatrix('u1old', nnodes, dim, "Symbol")        # u1(i,j) is the previous displacement of node i component j at domain 1
        u2old = custom_sympy_fe_utilities.DefineMatrix('u2old', nnodes_master, dim, "Symbol") # u2(i,j) is the previous displacement of node i component j at domain 2
        LM = custom_sympy_fe_utilities.DefineMatrix('LM', nnodes, dim, "Symbol")              # Lm are the Lagrange multipliers

        # Normal and tangets of the slave
        if normalvar == 1:
            NormalSlave = custom_sympy_fe_utilities.DefineMatrix('NormalSlave', nnodes, dim)
        else:
            NormalSlave = custom_sympy_fe_utilities.DefineMatrix('NormalSlave', nnodes, dim, "Symbol")

        # The resultant tangent
        if normalvar == 1:
            TangentSlave = custom_sympy_fe_utilities.DefineMatrix('TangentSlave', nnodes, dim)
        else:
            TangentSlave = custom_sympy_fe_utilities.DefineMatrix('TangentSlave', nnodes, dim, "Symbol")

        # Define test functions
        w1 = custom_sympy_fe_utilities.DefineMatrix('w1',nnodes,dim, "Symbol")
        w2 = custom_sympy_fe_utilities.DefineMatrix('w2',nnodes_master,dim, "Symbol")
        wLM = custom_sympy_fe_utilities.DefineMatrix('wLM',nnodes,dim, "Symbol")

        # Defining normal and tangent components
        LMNormal = custom_sympy_fe_utilities.DefineVector('LMNormal', nnodes)
        wLMNormal = custom_sympy_fe_utilities.DefineVector('wLMNormal', nnodes)

        # The resultant tangent LM
        LMTangent = custom_sympy_fe_utilities.DefineMatrix('LMTangent', nnodes, dim)
        wLMTangent = custom_sympy_fe_utilities.DefineMatrix('wLMTangent', nnodes, dim)

        # Defining additional variables
        NormalGap = custom_sympy_fe_utilities.DefineVector('NormalGap', nnodes)
        NormalwGap = custom_sympy_fe_utilities.DefineVector('NormalwGap', nnodes)
        TangentSlipNonObjective = custom_sympy_fe_utilities.DefineMatrix('TangentSlipNonObjective', nnodes, dim)
        TangentwSlipNonObjective = custom_sympy_fe_utilities.DefineMatrix('TangentwSlipNonObjective', nnodes, dim)
        TangentSlipObjective = custom_sympy_fe_utilities.DefineMatrix('TangentSlipObjective', nnodes, dim)
        TangentwSlipObjective = custom_sympy_fe_utilities.DefineMatrix('TangentwSlipObjective', nnodes, dim)
        DOperator = custom_sympy_fe_utilities.DefineMatrix('DOperator', nnodes, nnodes)
        MOperator = custom_sympy_fe_utilities.DefineMatrix('MOperator', nnodes, nnodes_master)
        DOperatorold = custom_sympy_fe_utilities.DefineMatrix('DOperatorold',nnodes,nnodes, "Symbol")
        MOperatorold = custom_sympy_fe_utilities.DefineMatrix('MOperatorold',nnodes,nnodes_master, "Symbol")

        # Define other parameters
        X1 = custom_sympy_fe_utilities.DefineMatrix('X1',nnodes,dim)
        X2 = custom_sympy_fe_utilities.DefineMatrix('X2',nnodes_master,dim)
        x1 = X1 + u1
        x2 = X2 + u2
        x1old = X1 + u1old
        x2old = X2 + u2old

        # Define other symbols
        mu = custom_sympy_fe_utilities.DefineVector('mu',nnodes, "Symbol")
        DynamicFactor = custom_sympy_fe_utilities.DefineVector('DynamicFactor',nnodes, "Symbol")
        PenaltyParameter = custom_sympy_fe_utilities.DefineVector('PenaltyParameter',nnodes, "Symbol")
        delta_time = sympy.Symbol('delta_time', positive=True)
        ScaleFactor = sympy.Symbol('ScaleFactor', positive=True)
        TangentFactor = sympy.Symbol('TangentFactor', positive=True)

        # Define variables list for later derivation
        u1_var = []
        u2_var = []
        LM_var = []
        custom_sympy_fe_utilities.CreateVariableMatrixList(u1_var, u1)
        custom_sympy_fe_utilities.CreateVariableMatrixList(u2_var, u2)
        u12_var=u1_var.copy()
        u1_LM_var=u1_var.copy()
        custom_sympy_fe_utilities.CreateVariableMatrixList(u12_var, u2)
        custom_sympy_fe_utilities.CreateVariableMatrixList(LM_var, LM)
        custom_sympy_fe_utilities.CreateVariableMatrixList(u1_LM_var, LM)
        all_var=u12_var.copy()
        custom_sympy_fe_utilities.CreateVariableMatrixList(all_var, LM)

        # Force the variables to be dependendant of the DOF
        if normalvar == 1:
            NormalSlave = custom_sympy_fe_utilities.DefineDofDependencyMatrix(NormalSlave, u1_var)
        DOperator = custom_sympy_fe_utilities.DefineDofDependencyMatrix(DOperator, u12_var) # If you consider Gitterle you need to keep the old operators
        MOperator = custom_sympy_fe_utilities.DefineDofDependencyMatrix(MOperator, u12_var)

        # Definitions of components
        for node in range(nnodes):
            LMNormal[node] = LM.row(node).dot(NormalSlave.row(node))
            wLMNormal[node] = wLM.row(node).dot(NormalSlave.row(node))

            # We calculate the LM tangent resultant
            #aux = LM.row(node).dot(TangentSlave.row(node)) * TangentSlave.row(node)
            #aux_w = wLM.row(node).dot(TangentSlave.row(node)) * TangentSlave.row(node)
            for idim in range(dim):
                #LMTangent[node,idim] = aux[idim]
                #wLMTangent[node,idim] = aux_w[idim]
                LMTangent[node,idim] = LM[node,idim] - LMNormal[node] * NormalSlave[node, idim]
                wLMTangent[node,idim] = wLM[node,idim] - wLMNormal[node] * NormalSlave[node, idim]

        ## Explicit definition tangent
        #for node in range(nnodes):
            #aux = LMTangent.row(node)/real_norm(LMTangent.row(node))
            #for idim in range(dim):
                #TangentSlave[node,idim] = aux[idim]

        # Defining the normal NormalGap and tangent slip
        Dx1Mx2 = DOperator * x1 - MOperator * x2
        #Dx1oldMx2old = DOperator * x1old - MOperator * x2old
        Dw1Mw2 = DOperator * w1 - MOperator * w2

        # Delta operations
        DDeltax1MDeltax2 = DOperator * (x1 - x1old) - MOperator * (x2 - x2old)
        DeltaDx1DeltaMx2 = (DOperator - DOperatorold) * x1 - (MOperator - MOperatorold) * x2
        DeltaDw1DeltaMw2 = (DOperator - DOperatorold) * w1 - (MOperator - MOperatorold) * w2
        for node in range(nnodes):
            NormalGap[node] = - Dx1Mx2.row(node).dot(NormalSlave.row(node))
            NormalwGap[node] = Dw1Mw2.row(node).dot(NormalSlave.row(node))
            #gap_time_derivative = - Dx1Mx2.row(node)/delta_time
            gap_time_derivative_non_objective = - DDeltax1MDeltax2.row(node)/delta_time
            gap_time_derivative_non_objective_w = Dw1Mw2.row(node)/delta_time
            gap_time_derivative_objective = DeltaDx1DeltaMx2.row(node)/delta_time
            gap_time_derivative_objective_w = - DeltaDw1DeltaMw2.row(node)/delta_time

            # Direct computation
            auxTangentSlipNonObjective = delta_time * (gap_time_derivative_non_objective - gap_time_derivative_non_objective.dot(NormalSlave.row(node)) * NormalSlave.row(node))
            auxTangentwSlipNonObjective = delta_time * (gap_time_derivative_non_objective_w - gap_time_derivative_non_objective_w.dot(NormalSlave.row(node)) * NormalSlave.row(node))
            auxTangentSlipObjective = delta_time * (gap_time_derivative_objective - gap_time_derivative_objective.dot(NormalSlave.row(node)) * NormalSlave.row(node))
            auxTangentwSlipObjective = delta_time * (gap_time_derivative_objective_w - gap_time_derivative_objective_w.dot(NormalSlave.row(node)) * NormalSlave.row(node))

            ## Enforced
            #auxTangentSlipNonObjective = delta_time * gap_time_derivative_non_objective.dot(TangentSlave.row(node)) * TangentSlave.row(node)
            #auxTangentwSlipNonObjective = delta_time * gap_time_derivative_non_objective_w.dot(TangentSlave.row(node)) * TangentSlave.row(node)
            #auxTangentSlipObjective = delta_time * gap_time_derivative_objective.dot(TangentSlave.row(node)) * TangentSlave.row(node)
            #auxTangentwSlipObjective = delta_time * gap_time_derivative_objective_w.dot(TangentSlave.row(node)) * TangentSlave.row(node)
            for idim in range(dim):
                TangentSlipNonObjective[node, idim] = auxTangentSlipNonObjective[idim]
                TangentwSlipNonObjective[node, idim] = auxTangentwSlipNonObjective[idim]
                TangentSlipObjective[node, idim] = auxTangentSlipObjective[idim]
                TangentwSlipObjective[node, idim] = auxTangentwSlipObjective[idim]

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

        # Compute galerkin functional # NOTE: This is for Galerkin functional
        lhs_string += lhs_template_begin_string
        rhs_string += rhs_template_begin_string
        if debug_counter == 0:
            for node in range(nnodes):
                for slip in range(5):
                    rv_galerkin = 0
                    if slip == 0: # Inactive
                        rv_galerkin -= ScaleFactor**2 / PenaltyParameter[node] * LMNormal[node] * wLMNormal[node]
                        rv_galerkin -= ScaleFactor**2 / (PenaltyParameter[node] * TangentFactor) * (LMTangent.row(node)).dot(wLMTangent.row(node))
                        #rv_galerkin -= ScaleFactor**2/PenaltyParameter[node] * (wLM.row(node).dot(LM.row(node)))
                    else:
                        rv_galerkin += ScaleFactor * NormalGap[node] * wLMNormal[node]

                        augmented_normal_contact_pressure = (ScaleFactor * LMNormal[node] + PenaltyParameter[node] * NormalGap[node])
                        normal_augmented_contact_pressure = augmented_normal_contact_pressure * NormalSlave.row(node)
                        #rv_galerkin += DynamicFactor[node] * augmented_normal_contact_pressure * NormalwGap[node]

                        augmented_lm = (ScaleFactor * LM.row(node) + PenaltyParameter[node] * NormalGap[node] * NormalSlave.row(node))

                        if slip == 1 or slip == 2: # Slip
                            augmented_tangent_contact_pressure = - mu[node] * augmented_normal_contact_pressure * TangentSlave.row(node)
                            modified_augmented_tangent_lm = ScaleFactor * LMTangent.row(node) - augmented_tangent_contact_pressure

                            #rv_galerkin -= (ScaleFactor / (PenaltyParameter[node] * TangentFactor)) * modified_augmented_tangent_lm.dot(wLMTangent.row(node))

                            #rv_galerkin -= ScaleFactor**2/PenaltyParameter[node] * (wLMTangent.row(node).dot(LMTangent.row(node))) # Frictionless
                            rv_galerkin -= ScaleFactor**2/PenaltyParameter[node] * (wLMTangent.row(node).dot(LMTangent.row(node) - augmented_tangent_contact_pressure/ScaleFactor))

                            #if slip == 1: # Objective
                                #rv_galerkin += DynamicFactor[node] * augmented_tangent_contact_pressure.dot(TangentwSlipObjective.row(node))
                            #else:
                                #rv_galerkin += DynamicFactor[node] * augmented_tangent_contact_pressure.dot(TangentwSlipNonObjective.row(node))
                        else: # Stick
                            if slip == 3: # Objective
                                augmented_tangent_contact_pressure = ScaleFactor * LMTangent.row(node) + TangentFactor * PenaltyParameter[node] * TangentSlipObjective.row(node)

                                augmented_lm += TangentFactor * PenaltyParameter[node] * TangentSlipObjective.row(node)

                                rv_galerkin += ScaleFactor * (TangentSlipObjective.row(node)).dot(wLMTangent.row(node))
                                #rv_galerkin += DynamicFactor[node] * augmented_tangent_contact_pressure.dot(TangentwSlipObjective.row(node))
                            else: # Non-Objective
                                augmented_tangent_contact_pressure = ScaleFactor * LMTangent.row(node) + TangentFactor * PenaltyParameter[node] * TangentSlipNonObjective.row(node)

                                augmented_lm += TangentFactor * PenaltyParameter[node] * TangentSlipNonObjective.row(node)

                                rv_galerkin += ScaleFactor * (TangentSlipNonObjective.row(node)).dot(wLMTangent.row(node))
                                #rv_galerkin += DynamicFactor[node] * augmented_tangent_contact_pressure.dot(TangentwSlipNonObjective.row(node))

                        rv_galerkin += DynamicFactor[node] * (augmented_lm).dot(Dw1Mw2.row(node))
                        #rv_galerkin += DynamicFactor[node] * (normal_augmented_contact_pressure + augmented_tangent_contact_pressure).dot(Dw1Mw2.row(node))

                    if do_simplifications:
                        rv_galerkin = sympy.simplify(rv_galerkin)

                    #############################################################################
                    # Complete functional
                    rv = sympy.Matrix(sympy.zeros(1, 1))
                    rv[0,0] = rv_galerkin

                    rhs,lhs = custom_sympy_fe_utilities.Compute_RHS_and_LHS(rv.copy(), testfunc, dofs, False)
                    print("LHS= ", lhs.shape)
                    print("RHS= ", rhs.shape)
                    print("LHS and RHS have been created!")

                    lhs_out = custom_sympy_fe_utilities.OutputMatrix_CollectingFactorsNonZero(lhs, "lhs", mode, 1, number_dof)
                    rhs_out = custom_sympy_fe_utilities.OutputVector_CollectingFactorsNonZero(rhs, "rhs", mode, 1, number_dof)
                    print("Substitution strings are ready....")

                    if slip == 0:
                        lhs_string += "    \n    // NODE " + str(node) + "\n"
                        lhs_string += "    if (r_geometry["+str(node)+"].IsNot(ACTIVE)) { // INACTIVE\n    "
                        lhs_string += lhs_out.replace("\n","\n    ")
                    elif slip == 1 or slip == 2:
                        if slip == 1:
                            lhs_string += "} else if (r_geometry["+str(node)+"].Is(SLIP)) { // ACTIVE-SLIP\n        if (is_objetive) { // OBJECTIVE-SLIP\n        "
                            lhs_string += lhs_out.replace("\n","\n        ")
                        else:
                            lhs_string += "} else { // NONOBJECTIVE-SLIP\n        "
                            lhs_string += lhs_out.replace("\n","\n        ")
                            lhs_string += "}\n    "
                    else:
                        if slip == 3:
                            lhs_string += "} else { // ACTIVE-STICK\n        if (is_objetive) { // OBJECTIVE-STICK\n        "
                            lhs_string += lhs_out.replace("\n","\n        ")
                        else:
                            lhs_string += "} else { // NONOBJECTIVE-STICK\n        "
                            lhs_string += lhs_out.replace("\n","\n        ")
                    if slip == 4:
                        lhs_string += "}\n    }\n"

                    if slip == 0:
                        rhs_string += "    \n    // NODE " + str(node) + "\n"
                        rhs_string += "    if (r_geometry["+str(node)+"].IsNot(ACTIVE)) { // INACTIVE\n    "
                        rhs_string += rhs_out.replace("\n","\n    ")
                    elif slip == 1 or slip == 2:
                        if slip == 1:
                            rhs_string += "} else if (r_geometry["+str(node)+"].Is(SLIP)) { // ACTIVE-SLIP\n        if (is_objetive) { // OBJECTIVE-SLIP\n        "
                            rhs_string += rhs_out.replace("\n","\n        ")
                        else:
                            rhs_string += "} else { // NONOBJECTIVE-SLIP\n        "
                            rhs_string += rhs_out.replace("\n","\n        ")
                            rhs_string += "}\n    "
                    else:
                        if slip == 3:
                            rhs_string += "} else { // ACTIVE-STICK\n        if (is_objetive) { // OBJECTIVE-STICK\n        "
                            rhs_string += rhs_out.replace("\n","\n        ")
                        else:
                            rhs_string += "} else { // NONOBJECTIVE-STICK\n        "
                            rhs_string += rhs_out.replace("\n","\n        ")
                    if slip == 4:
                        rhs_string += "}\n    }\n"

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
            first_input = open("ALM_frictional_mortar_contact_condition_template.cpp",'r').read()
            lhs_string += "// replace_lhs"
            outputstring = first_input.replace("// replace_lhs", lhs_string)
        else:
            input = open("ALM_frictional_mortar_contact_condition.cpp",'r').read()
            if output_count < total_combs:
                lhs_string += "// replace_lhs"
            outputstring = input.replace("// replace_lhs", lhs_string)

        if (output_count < total_combs):
            rhs_string += "// replace_rhs"

        outputstring = outputstring.replace("// replace_rhs", rhs_string)
        output = open("ALM_frictional_mortar_contact_condition.cpp",'w')
        output.write(outputstring)
        output.close()

print("Strings have been replaced...")

print("Process Finished..................")
