from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from sympy import *
from custom_sympy_fe_utilities import *
import operator 

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
        
## Debug
#dim_combinations = [2]
#nnodes_combinations = [2]

dim_combinations = [2,3,3]
nnodes_combinations = [2,3,4]

def real_norm(input):

    output = 0
    
    for i in range(input.shape[1]):
        output += input[i]**2
        
    output = real_root(output, 2)
        
    return output

lhs_string = ""
lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nbounded_matrix<double, MatrixSize, MatrixSize> AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes, TNormalVariation>::CalculateLocalLHS(\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const DerivativeDataType& rDerivativeData,\n        const unsigned int rActiveInactive\n        )\n{\n    bounded_matrix<double,MatrixSize,MatrixSize> lhs = ZeroMatrix(MatrixSize, MatrixSize);\n\n    // Initialize values\n    const bounded_matrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    const bounded_matrix<double, TNumNodes, TDim>& u1old = rDerivativeData.u1old;\n    const bounded_matrix<double, TNumNodes, TDim>& u2 = rDerivativeData.u2;\n    const bounded_matrix<double, TNumNodes, TDim>& u2old = rDerivativeData.u2old;\n    const bounded_matrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    const bounded_matrix<double, TNumNodes, TDim>& X2 = rDerivativeData.X2;\n    \n    const bounded_matrix<double, TNumNodes, TDim>& LM = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); \n    \n    const bounded_matrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n// The ALM parameters\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    const double TangentFactor = rDerivativeData.TangentFactor;\n    \n    // Mortar operators\n    const bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n    // Mortar operators derivatives\n    const array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n    const array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n\n    // We get the friction coefficient\n    const array_1d<double, TNumNodes>& mu = GetFrictionCoefficient();\n\n"

lhs_template_end_string = "\n\n    return lhs;\n}\n"

rhs_string = ""
rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\narray_1d<double, MatrixSize> AugmentedLagrangianMethodFrictionalMortarContactCondition<TDim,TNumNodes, TNormalVariation>::CalculateLocalRHS(\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const DerivativeDataType& rDerivativeData,\n        const unsigned int rActiveInactive\n        )\n{\n    array_1d<double,MatrixSize> rhs(MatrixSize, 0.0);\n\n    // Initialize values\n    const bounded_matrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    const bounded_matrix<double, TNumNodes, TDim>& u1old = rDerivativeData.u1old;\n    const bounded_matrix<double, TNumNodes, TDim>& u2 = rDerivativeData.u2;\n    const bounded_matrix<double, TNumNodes, TDim>& u2old = rDerivativeData.u2old;\n    const bounded_matrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    const bounded_matrix<double, TNumNodes, TDim>& X2 = rDerivativeData.X2;\n    \n    const bounded_matrix<double, TNumNodes, TDim>& LM = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), VECTOR_LAGRANGE_MULTIPLIER, 0); \n    \n    const bounded_matrix<double, TNumNodes, TDim>& NormalSlave = rDerivativeData.NormalSlave;\n// The ALM parameters\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    const double TangentFactor = rDerivativeData.TangentFactor;\n    \n    // Mortar operators\n    const bounded_matrix<double, TNumNodes, TNumNodes>& MOperator = rMortarConditionMatrices.MOperator;\n    const bounded_matrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n    // We get the friction coefficient\n\n    const array_1d<double, TNumNodes>& mu = GetFrictionCoefficient();\n\n"

rhs_template_end_string = "\n\n    return rhs;\n}\n"

for normalvar in range(2):
    
    if normalvar == 0:
        normalvarstring = "false"
    else:
        normalvarstring = "true"
        
    if normalvar == 1:
        lhs_template_begin_string += "   const array_1d<bounded_matrix<double, TNumNodes, TDim>,  (TNumNodes * TDim)> DeltaNormalSlave = rDerivativeData.DeltaNormalSlave;\n\n"

    for dim, nnodes in zip(dim_combinations, nnodes_combinations):
        
        number_dof = dim * (3 * nnodes)

        #Defining the unknowns
        u1 = DefineMatrix('u1',nnodes,dim) #u1(i,j) is the current displacement of node i component j at domain 1
        u2 = DefineMatrix('u2',nnodes,dim) #u2(i,j) is the current displacement of node i component j at domain 2
        u1old = DefineMatrix('u1old',nnodes,dim) #u1(i,j) is the previous displacement of node i component j at domain 1
        u2old = DefineMatrix('u2old',nnodes,dim) #u2(i,j) is the previous displacement of node i component j at domain 2
        LM = DefineMatrix('LM',nnodes,dim)
        # Normal and tangets of the slave
        NormalSlave = DefineMatrix('NormalSlave',nnodes,dim)
        
        # The resultant tangent
        TangentSlave = DefineMatrix('TangentSlave',nnodes,dim)
        
        # Define test functions
        w1 = DefineMatrix('w1',nnodes,dim)
        w2 = DefineMatrix('w2',nnodes,dim)
        wLM = DefineMatrix('wLM',nnodes,dim)
        
        # Defining normal and tangent components 
        LMNormal = DefineVector('LMNormal',nnodes)
        wLMNormal = DefineVector('wLMNormal',nnodes)
        
        # The resultant tangent LM
        LMTangent = DefineMatrix('LMTangent',nnodes,dim)
        wLMTangent = DefineMatrix('wLMTangent',nnodes,dim)
        
        for node in range(nnodes):
            LMNormal[node] = LM.row(node) * NormalSlave.row(node).transpose()
            wLMNormal[node] = wLM.row(node) * NormalSlave.row(node).transpose()
            
            # We calculate the LM tangent resultant
            for idim in range(dim):
                LMTangent[node,idim] = LM[node,idim] - LMNormal[node] * NormalSlave[node,idim]
                wLMTangent[node,idim] = wLM[node,idim] - wLMNormal[node] * NormalSlave[node,idim]
            
            tangentnorm = real_norm(LMTangent.row(node))
            for idim in range(dim):
                TangentSlave[node,idim] = LMTangent[node,idim]/tangentnorm
            
        # Defining additional variables
        NormalGap = DefineVector('NormalGap',nnodes) 
        TangentSlip = DefineMatrix('TangentSlip',nnodes, dim) 
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
        PenaltyParameter  = DefineVector('PenaltyParameter',nnodes) 
        ScaleFactor = Symbol('ScaleFactor', positive=True)
        TangentFactor = Symbol('TangentFactor', positive=True)

        # Define variables list for later derivation
        u1_var = []
        u2_var = []
        LM_var = []
        CreateVariableMatrixList(u1_var, u1)
        CreateVariableMatrixList(u2_var, u2)
        u12_var=u1_var.copy()
        u1_LM_var=u1_var.copy()
        CreateVariableMatrixList(u12_var, u2)
        CreateVariableMatrixList(LM_var, LM)
        CreateVariableMatrixList(u1_LM_var, LM)
        all_var=u12_var.copy()
        CreateVariableMatrixList(all_var, LM)

        # Force the variables to be dependendant of the DOF
        if normalvar == 1:
            NormalSlave = DefineDofDependencyMatrix(NormalSlave, u1_var)
        DOperator = DefineDofDependencyMatrix(DOperator, u12_var) # If you consider Gitterle you need to keep the old operators
        MOperator = DefineDofDependencyMatrix(MOperator, u12_var)

        # Defining the normal NormalGap and tangent slip
        Dx1Mx2 = DOperator * x1 - MOperator * x2
        Dx1oldMx2old = DOperator * x1old - MOperator * x2old
        Dw1Mw2 = DOperator * w1 - MOperator * w2
        for node in range(nnodes): 
            NormalGap[node]  = Dx1Mx2.row(node) * - NormalSlave.row(node).transpose()
            auxTangentSlip = (Dx1oldMx2old.row(node) * TangentSlave.row(node).transpose())
            for idim in range(dim):
                TangentSlip[node,idim]  = auxTangentSlip * TangentSlave[node,idim]

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

        # We compute the augmented tangent LM
        hatLMTangent = DefineMatrix('hatLMTangent',nnodes,dim)
        HatTangentSlave = DefineMatrix('HatTangentSlave',nnodes,dim)
        for node in range(nnodes):
            for idim in range(dim):
                hatLMTangent[node,idim] = ScaleFactor * LMTangent[node,idim] + TangentSlip[node,idim] * PenaltyParameter[node] * TangentFactor
        
        # Now we can compute the resultant tangent
        HatTangentSlave = hatLMTangent.copy()
        for node in range(nnodes):
            AuxHatTangentNorm = real_norm(HatTangentSlave.row(node))
            for idim in range(dim):
                HatTangentSlave[node,idim] /= AuxHatTangentNorm

        # Compute galerkin functional # NOTE: Maybe you can define a different penalty and scale factor in the tangent direction
        lhs_string += lhs_template_begin_string
        rhs_string += rhs_template_begin_string
        for node in range(nnodes):
            for slip in range(3):
                rv_galerkin = 0
                if (slip == 0):  
                    rv_galerkin -=  ScaleFactor**2.0/PenaltyParameter[node] * LMNormal[node] * wLMNormal[node]
                    rv_galerkin -= (ScaleFactor**2.0/(PenaltyParameter[node] * TangentFactor) * LMTangent.row(node) * wLMTangent.row(node).transpose())[0,0]
                else:
                    rv_galerkin += ((((ScaleFactor * LMNormal[node] + PenaltyParameter[node] * NormalGap[node]) * NormalSlave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                    rv_galerkin +=  ScaleFactor * NormalGap[node] * wLMNormal[node]
                    
                    if (slip == 1): # Slip 
                        rv_galerkin -= (((mu[node] * (ScaleFactor * LMNormal[node] + PenaltyParameter[node] * NormalGap[node]) * HatTangentSlave.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                        rv_galerkin -=  (ScaleFactor**2.0/(PenaltyParameter[node] * TangentFactor) * (ScaleFactor * LMTangent.row(node) + mu[node] * (ScaleFactor * LMNormal[node] + PenaltyParameter[node] * NormalGap[node]) * HatTangentSlave.row(node)) * wLMTangent.row(node).transpose())[0,0] 
                    else: # Stick 
                        rv_galerkin += (((ScaleFactor * LMTangent.row(node) + PenaltyParameter[node] * TangentSlip.row(node))) * Dw1Mw2.row(node).transpose())[0,0]
                        rv_galerkin +=  (ScaleFactor * TangentSlip.row(node) * wLMTangent.row(node).transpose())[0,0]

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

                lhs_out = OutputMatrix_CollectingFactorsNonZero(lhs,"lhs", mode, 1, number_dof)
                rhs_out = OutputVector_CollectingFactorsNonZero(rhs,"rhs", mode, 1, number_dof)
                print("Substitution strings are ready....")

                if (slip == 0):
                    lhs_string += "    \n    // NODE " + str(node) + "\n"
                    lhs_string += "    if (this->GetGeometry()["+str(node)+"].Is(ACTIVE) == false ) // INACTIVE\n    {\n    "
                elif (slip == 1):
                    lhs_string += "    else if (this->GetGeometry()["+str(node)+"].Is(SLIP) == true ) // ACTIVE-SLIP\n    {\n    "
                else:
                    lhs_string += "    else // ACTIVE-STICK\n    {\n    "
                    
                lhs_string += lhs_out.replace("\n","\n    ")
                lhs_string += "}\n"    
            
                if (slip == 0):
                    rhs_string += "    \n    // NODE " + str(node) + "\n"
                    rhs_string += "    if (this->GetGeometry()["+str(node)+"].Is(ACTIVE) == false ) // INACTIVE\n    {\n    "
                elif (slip == 1):
                    rhs_string += "    else if (this->GetGeometry()["+str(node)+"].Is(SLIP) == true ) // ACTIVE-SLIP\n    {\n    "
                else:
                    rhs_string += "    else // ACTIVE-STICK\n    {\n    "
                rhs_string += rhs_out.replace("\n","\n    ")
                rhs_string += "}\n"
                
        lhs_string += lhs_template_end_string
        rhs_string += rhs_template_end_string
                
        lhs_string = lhs_string.replace("TDim", str(dim))
        lhs_string = lhs_string.replace("TNumNodes", str(nnodes))
        lhs_string = lhs_string.replace("MatrixSize", str(lhs.shape[0]))
        lhs_string = lhs_string.replace("TNormalVariation", normalvarstring)
        lhs_string = lhs_string.replace("SIZEDERIVATIVES2", str(2 * (nnodes * dim)))

        rhs_string = rhs_string.replace("TDim", str(dim))
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
            var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(NormalSlave, "NormalSlave", "NormalSlave", u1_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
        var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(DOperator, "DOperator", "DOperator", u12_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")
        var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list = DefineVariableLists(MOperator, "MOperator", "MOperator", u12_var, var_strings, var_strings_subs, var_strings_aux_subs, der_var_strings, der_var_list, "matrix")

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

for dof in reversed(range(len(u1_var))):
    lhs_string = lhs_string.replace("DeltaNormalSlave"+str(dof), "DeltaNormalSlave["+str(dof)+"]")
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
