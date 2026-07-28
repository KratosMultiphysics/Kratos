//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes


// External includes


// Project includes
#include "elements/embedded_nodal_variable_calculation_element_simplex.h"


namespace Kratos {

template <class TVarType>
void EmbeddedNodalVariableCalculationElementSimplex<TVarType>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <class TVarType>
int EmbeddedNodalVariableCalculationElementSimplex<TVarType>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if(ErrorCode != 0) {
        return ErrorCode;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(auto &r_node : this->GetGeometry()) {
        if constexpr (std::is_same_v<TVarType, double>) {
            KRATOS_ERROR_IF(!r_node.SolutionStepsDataHas(NODAL_MAUX)) << "Missing NODAL_MAUX variable on solution step data for node " << r_node.Id() << std::endl;;
        } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
            KRATOS_ERROR_IF(!r_node.SolutionStepsDataHas(NODAL_VAUX)) << "Missing NODAL_VAUX variable on solution step data for node " << r_node.Id() << std::endl;;
        } else {
            KRATOS_ERROR << "Unsupported variable type" << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH("");
}

template <class TVarType>
void EmbeddedNodalVariableCalculationElementSimplex<TVarType>::CalculateLeftHandSide(
    MatrixType &rLeftHandSideMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    std::size_t expected_matrix_size = 0;

    if constexpr (std::is_same_v<TVarType, double>) { 
        expected_matrix_size = 2;
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        expected_matrix_size = 6;
    } else {
        KRATOS_ERROR << "Unsupported variable type" << std::endl;
    }

    // Check size
    if (rLeftHandSideMatrix.size1() != expected_matrix_size || rLeftHandSideMatrix.size2() != expected_matrix_size) {
        rLeftHandSideMatrix.resize(expected_matrix_size, expected_matrix_size, false);
    }

    // Initialize LHS. This is required since not all the entries of the matrix are iterated
    rLeftHandSideMatrix = ZeroMatrix(expected_matrix_size,expected_matrix_size);

    // Get the element shape function values from the normalized distance to node 0
    const auto &rN = this->GetDistanceBasedShapeFunctionValues();

    // Compute the Gramm matrix and gradient penalty term
    const double penalty = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT];
    std::array<double, 2> aux_penalty{{penalty, -penalty}};
    if constexpr (std::is_same_v<TVarType, double>) { 
        for (unsigned int i = 0; i < 2; ++i) {
            for (unsigned int j = 0; j < 2; ++j) {
                rLeftHandSideMatrix(i, j) = rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j];
            }
        }
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        for (unsigned int i = 0; i < 2; ++i) {
            for (unsigned int j = 0; j < 2; ++j) {
                for (unsigned int k = 0; k < 3; ++k) {
                    rLeftHandSideMatrix(i * 3 + k, j * 3 + k) = rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j];
                }
            }
        }
    }
}

template <class TVarType>
void EmbeddedNodalVariableCalculationElementSimplex<TVarType>::CalculateRightHandSide(
    VectorType &rRigthHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    std::size_t expected_matrix_size = 0;

    if constexpr (std::is_same_v<TVarType, double>) { 
        expected_matrix_size = 2;
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        expected_matrix_size = 6;
    } else {
        KRATOS_ERROR << "Unsupported variable type" << std::endl;
    }

    // Check size
    if (rRigthHandSideVector.size() != expected_matrix_size) {
        rRigthHandSideVector.resize(expected_matrix_size, false);
    }

    // Get the element shape function values from the normalized distance to node 0
    const auto &r_geom = this->GetGeometry();
    const auto &rN = this->GetDistanceBasedShapeFunctionValues();

    // Compute the data and penalty Right Hand Side contributions
    const double penalty = rCurrentProcessInfo[GRADIENT_PENALTY_COEFFICIENT];
    std::array<double, 2> aux_penalty{{penalty, -penalty}};

    if constexpr (std::is_same_v<TVarType, double>) { 
        const double &rData = this->GetValue(NODAL_MAUX);
        for (unsigned int i = 0; i < 2; ++i) {
            rRigthHandSideVector(i) = rN[i] * rData;
            for (unsigned int j = 0; j < 2; ++j) {
                rRigthHandSideVector(i) -=  (rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j]) * r_geom[j].FastGetSolutionStepValue(NODAL_MAUX);
            }
        }
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        const array_1d<double,3> &rData = this->GetValue(NODAL_VAUX);
        for (unsigned int i = 0; i < 2; ++i) {
            const auto &r_aux = r_geom[i].FastGetSolutionStepValue(NODAL_VAUX);
            for (unsigned int k = 0; k < 3; ++k) {
                rRigthHandSideVector(i * 3 + k) = rN[i] * rData(k);
                for (unsigned int j = 0; j < 2; ++j) {
                    rRigthHandSideVector(i * 3 + k) -= (rN[i] * rN[j] + aux_penalty[i] * aux_penalty[j]) * r_aux[k];
                }
            }
        }
    }
}

template <class TVarType>
void EmbeddedNodalVariableCalculationElementSimplex<TVarType>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    std::size_t expected_matrix_size = 0;

    if constexpr (std::is_same_v<TVarType, double>) { 
        expected_matrix_size = 2;
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        expected_matrix_size = 6;
    } else {
        KRATOS_ERROR << "Unsupported variable type" << std::endl;
    }

    if (rResult.size() != expected_matrix_size) {
        rResult.resize(expected_matrix_size, false);
    }

    if constexpr (std::is_same_v<TVarType, double>) { 
        const unsigned int pos = (this->GetGeometry())[0].GetDofPosition(NODAL_MAUX);
        for (unsigned int i = 0; i < 2; i++) {
            rResult[i] = (this->GetGeometry())[i].GetDof(NODAL_MAUX, pos).EquationId();
        }
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        const unsigned int x_pos = (this->GetGeometry())[0].GetDofPosition(NODAL_VAUX_X);
        for (unsigned int i = 0; i < 2; i++) {
            rResult[i * 3] = GetGeometry()[i].GetDof(NODAL_VAUX_X, x_pos).EquationId();
            rResult[i * 3 + 1] = GetGeometry()[i].GetDof(NODAL_VAUX_Y, x_pos + 1).EquationId();
            rResult[i * 3 + 2] = GetGeometry()[i].GetDof(NODAL_VAUX_Z, x_pos + 2).EquationId();
        }
    }
}

template <class TVarType>
void EmbeddedNodalVariableCalculationElementSimplex<TVarType>::GetDofList(
    DofsVectorType &rElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    std::size_t expected_matrix_size = 0;

    if constexpr (std::is_same_v<TVarType, double>) { 
        expected_matrix_size = 2;
    } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
        expected_matrix_size = 6;
    } else {
        KRATOS_ERROR << "Unsupported variable type" << std::endl;
    }

    if (rElementalDofList.size() != expected_matrix_size) {
        rElementalDofList.resize(expected_matrix_size);
    }

    for (unsigned int i = 0; i < 2; i++) {
        if constexpr (std::is_same_v<TVarType, double>) { 
            rElementalDofList[i] = (this->GetGeometry())[i].pGetDof(NODAL_MAUX);
        } else if constexpr (std::is_same_v<TVarType, array_1d<double, 3>>) {
            rElementalDofList[i * 3] = (this->GetGeometry())[i].pGetDof(NODAL_VAUX_X);
            rElementalDofList[i * 3 + 1] = (this->GetGeometry())[i].pGetDof(NODAL_VAUX_Y);
            rElementalDofList[i * 3 + 2] = (this->GetGeometry())[i].pGetDof(NODAL_VAUX_Z);
        }
    }
}

template <class TVarType>
const array_1d<double, 2> EmbeddedNodalVariableCalculationElementSimplex<TVarType>::GetDistanceBasedShapeFunctionValues()
{
    const double d = this->GetValue(DISTANCE);
    array_1d<double, 2> N;
    N[0] = 1.0 - d;
    N[1] = d;
    return N;
}

template class KRATOS_API(KRATOS_CORE) EmbeddedNodalVariableCalculationElementSimplex<double>;
template class KRATOS_API(KRATOS_CORE) EmbeddedNodalVariableCalculationElementSimplex<Kratos::array_1d<double, 3ul>>;

} // namespace Kratos.
