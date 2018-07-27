from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ContactStructuralMechanicsApplication  import *

from sympy import *
from sympy.physics.vector import *
from custom_sympy_fe_utilities import *
import operator

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

### Debug
#dim_combinations = [2]
#nnodes_combinations = [2]

dim_combinations = [2,3,3]
nnodes_combinations = [2,3,4]

lhs_string = ""
lhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nvoid DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation>::CalculateLocalLHS(\n        Matrix& rLocalLHS,\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const DerivativeDataType& rDerivativeData,\n        const IndexType rActiveInactive,\n    const ProcessInfo& rCurrentProcessInfo\n        )\n{\n    // Initialize values\n    // const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    // const BoundedMatrix<double, TNumNodes, TDim>& u2 = rDerivativeData.u2;\n    // const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    // const BoundedMatrix<double, TNumNodes, TDim>& X2 = rDerivativeData.X2;\n    \n    const array_1d<double, TNumNodes> LMSlaveNormal = MortarUtilities::GetVariableVector<TNumNodes>(this->GetGeometry(), LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, 0);\n    const array_1d<double, TNumNodes> LMMasterNormal = MortarUtilities::GetVariableVector<TNumNodes>(this->GetPairedGeometry(), LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, 0);\n\n    // The ALM parameters\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    // const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    \n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n    // Mortar operators derivatives\n    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaMOperator = rMortarConditionMatrices.DeltaMOperator;\n    const array_1d<BoundedMatrix<double, TNumNodes, TNumNodes>, SIZEDERIVATIVES2>& DeltaDOperator = rMortarConditionMatrices.DeltaDOperator;\n\n"

lhs_template_end_string = "\n\n}\n"

rhs_string = ""
rhs_template_begin_string = "\n/***********************************************************************************/\n/***********************************************************************************/\n\ntemplate<>\nvoid DoubleAugmentedLagrangianMethodFrictionlessMortarContactCondition<TDim,TNumNodes, TNormalVariation>::CalculateLocalRHS(\n        Vector& rLocalRHS,\n        const MortarConditionMatrices& rMortarConditionMatrices,\n        const DerivativeDataType& rDerivativeData,\n        const IndexType rActiveInactive,\n    const ProcessInfo& rCurrentProcessInfo\n        )\n{\n    // Initialize values\n    // const BoundedMatrix<double, TNumNodes, TDim>& u1 = rDerivativeData.u1;\n    // const BoundedMatrix<double, TNumNodes, TDim>& u2 = rDerivativeData.u2;\n    // const BoundedMatrix<double, TNumNodes, TDim>& X1 = rDerivativeData.X1;\n    // const BoundedMatrix<double, TNumNodes, TDim>& X2 = rDerivativeData.X2;\n    \n    const array_1d<double, TNumNodes> LMSlaveNormal = MortarUtilities::GetVariableVector<TNumNodes>(this->GetGeometry(), LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, 0);\n    const array_1d<double, TNumNodes> LMMasterNormal = MortarUtilities::GetVariableVector<TNumNodes>(this->GetPairedGeometry(), LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, 0);\n\n    // The ALM parameters\n    const double ScaleFactor = rDerivativeData.ScaleFactor;\n    // const array_1d<double, TNumNodes>& PenaltyParameter = rDerivativeData.PenaltyParameter;\n    \n    // Mortar operators\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& MOperator = rMortarConditionMatrices.MOperator;\n    const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rMortarConditionMatrices.DOperator;\n\n"

rhs_template_end_string = "\n\n}\n"

for normalvar in range(2):

    if normalvar == 0:
        normalvarstring = "false"
    else:
        normalvarstring = "true"

    for dim, nnodes in zip(dim_combinations, nnodes_combinations):

        from sympy.utilities.iterables import ibin
        active_inactive_combinations = list(ibin(nnodes, 'all'))

        active_inactive_comb = 0

        for active_inactive in active_inactive_combinations: # Change output in combination of this

            active_inactive_comb += 1
            number_dof = dim * (2 * nnodes) + 2 * nnodes

            #Defining the unknowns
            u1 = DefineMatrix('u1', nnodes, dim) #u1(i,j) is displacement of node i component j at domain 1
            u2 = DefineMatrix('u2', nnodes, dim) #u2(i,j) is displacement of node i component j at domain 2
            LMSlaveNormal = DefineVector('LMSlaveNormal', nnodes)
            LMMasterNormal = DefineVector('LMMasterNormal', nnodes)
            DOperator = DefineMatrix('DOperator', nnodes, nnodes)
            MOperator = DefineMatrix('MOperator', nnodes, nnodes)

            X1 = DefineMatrix('X1', nnodes, dim)
            X2 = DefineMatrix('X2', nnodes, dim)
            x1 = X1 + u1
            x2 = X2 + u2

            #Define other symbols
            PenaltyParameter  = DefineVector('PenaltyParameter', nnodes)
            ScaleFactor = Symbol('ScaleFactor', positive=True)

            # Define test functions
            w1 = DefineMatrix('w1', nnodes, dim)
            w2 = DefineMatrix('w2', nnodes, dim)
            wLMSlaveNormal = DefineVector('wLMSlaveNormal', nnodes)
            wLMMasterNormal = DefineVector('wLMMasterNormal', nnodes)

            # Define variables list for later derivation
            u1_var = []
            u2_var = []
            CreateVariableMatrixList(u1_var, u1)
            CreateVariableMatrixList(u2_var, u2)
            u12_var=u1_var.copy()
            CreateVariableMatrixList(u12_var, u2)

            # Force the variables to be dependendant of the DOF
            DOperator = DefineDofDependencyMatrix(DOperator, u12_var)
            MOperator = DefineDofDependencyMatrix(MOperator, u12_var)

            # Defining the normal NormalGap
            DLM1MLM2 = DOperator * LMSlaveNormal - MOperator * LMMasterNormal

            # Define dofs & test function vector
            dofs = Matrix( zeros(number_dof, 1) )
            testfunc = Matrix( zeros(number_dof, 1) )
            count = 0
            for i in range(0, nnodes):
                for k in range(0, dim):
                    dofs[count] = u2[i, k]
                    testfunc[count] = w2[i, k]
                    count += 1
            for i in range(0, nnodes):
                for k in range(0, dim):
                    dofs[count] = u1[i, k]
                    testfunc[count] = w1[i, k]
                    count += 1
            for i in range(0, nnodes):
                dofs[count] = LMSlaveNormal[i]
                testfunc[count] = wLMSlaveNormal[i]
                count += 1
            for i in range(0, nnodes):
                dofs[count] = LMMasterNormal[i]
                testfunc[count] = wLMMasterNormal[i]
                count += 1
            print("dofs = ", dofs)
            print("testfunc = ", testfunc)

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
                    rv_galerkin += ScaleFactor * DLM1MLM2[node] * wLMSlaveNormal[node] # TODO: Think about adding the penalty x gap contribution

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

            lhs_out = OutputMatrix_CollectingFactorsNonZero(lhs, "lhs", mode, 1, number_dof)
            rhs_out = OutputVector_CollectingFactorsNonZero(rhs, "rhs", mode, 1, number_dof)
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
lhs_string = lhs_string.replace("lhs(", "rLocalLHS(")
rhs_string = rhs_string.replace("rhs[", "rLocalRHS[")

for dof in reversed(range(len(u12_var))):
    lhs_string = lhs_string.replace("DeltaDOperator" + str(dof), "DeltaDOperator[" + str(dof) + "]")
for dof in reversed(range(len(u12_var))):
    lhs_string = lhs_string.replace("DeltaMOperator" + str(dof), "DeltaMOperator[" + str(dof) + "]")

#############################################################################
################################# FINAL SAVING ##############################
#############################################################################

input = open("DALM_frictionless_mortar_contact_condition_template.cpp",'r').read()
outputstring = input.replace("// replace_lhs", lhs_string)
outputstring = outputstring.replace("// replace_rhs", rhs_string)
output = open("DALM_frictionless_mortar_contact_condition.cpp",'w')
output.write(outputstring)
output.close()

print("Strings have been replaced...")

print("Process Finished..................")
