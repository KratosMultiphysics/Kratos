/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//
//  Authors: Tobias Teschemacher
*/

// System includes
#include "includes/define.h"
#include "utilities/math_utils.h"

// External includes

// Project includes
#include "custom_elements/base_discrete_element.h"

#include "iga_application.h"
#include "iga_application_variables.h"

#include "geometries/geometry.h"

namespace Kratos
{
    void BaseDiscreteElement::Initialize()
    {
        KRATOS_TRY

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != 1)
            mConstitutiveLawVector.resize(1);

        InitializeMaterial();

        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::InitializeMaterial()
    {
        KRATOS_TRY

            if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
                mConstitutiveLawVector[0] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[0]->InitializeMaterial(GetProperties(),
                    GetGeometry(),
                    GetValue(SHAPE_FUNCTION_VALUES)
                );
            }
            else
                KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
                mConstitutiveLawVector[0]->ResetMaterial(GetProperties(),
                    GetGeometry(),
                    GetValue(SHAPE_FUNCTION_VALUES)
                );
        }

        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const unsigned int number_of_control_points = GetGeometry().size();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 3;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    };

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;

        const unsigned int number_of_control_points = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    };

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_ERROR << "You have called to the CalculateMassMatrix() from the base class BaseDiscreteElement" << std::endl;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY;
        unsigned int number_of_control_points = GetGeometry().size();

        // Resizing as needed the LHS
        unsigned int mat_size = number_of_control_points * 3;

        if (rDampingMatrix.size1() != mat_size)
            rDampingMatrix.resize(mat_size, mat_size, false);

        noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);

        // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        double alpha = 0.0;
        if (GetProperties().Has(RAYLEIGH_ALPHA))
            alpha = GetProperties()[RAYLEIGH_ALPHA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_ALPHA))
            alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

        double beta = 0.0;
        if (GetProperties().Has(RAYLEIGH_BETA))
            beta = GetProperties()[RAYLEIGH_BETA];
        else if (rCurrentProcessInfo.Has(RAYLEIGH_BETA))
            beta = rCurrentProcessInfo[RAYLEIGH_BETA];

        // Rayleigh Damping Matrix: alpha*M + beta*K

        // 2.-Calculate StiffnessMatrix:
        if (beta > 0.0)
        {
            MatrixType StiffnessMatrix = Matrix();
            VectorType ResidualVector = Vector();
            this->CalculateAll(StiffnessMatrix, ResidualVector, rCurrentProcessInfo, true, false);
            noalias(rDampingMatrix) += beta * StiffnessMatrix;
        }

        // 3.-Calculate MassMatrix:
        if (alpha > 0.0)
        {
            MatrixType MassMatrix = Matrix();
            this->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);
            noalias(rDampingMatrix) += alpha * MassMatrix;
        }

        KRATOS_CATCH("")
    }


    /***********************************************************************************/
    /// CalculateOnIntegrationPoints
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rOutput.size() != 1)
            rOutput.resize(1);

        bool flag = false;
        mConstitutiveLawVector[0]->GetValue(rVariable, flag);
        rOutput[0] = flag;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rOutput.size() != 1)
            rOutput.resize(1);

        mConstitutiveLawVector[0]->GetValue(rVariable, rOutput[0]);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rOutput.size() != 1)
            rOutput.resize(1);

        mConstitutiveLawVector[0]->GetValue(rVariable, rOutput[0]);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rOutput.size() != 1)
            rOutput.resize(1);

        mConstitutiveLawVector[0]->GetValue(rVariable, rOutput[0]);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::CalculateOnIntegrationPoints(
        const Variable<Matrix >& rVariable,
        std::vector< Matrix >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rOutput.size() != 1)
            rOutput.resize(1);

        mConstitutiveLawVector[0]->GetValue(rVariable, rOutput[0]);
    }

    /***********************************************************************************/
    /// SetValueOnIntegrationPoints
    /***********************************************************************************/
    void BaseDiscreteElement::SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        mConstitutiveLawVector[0]->SetValue(rVariable,
            rValues[0],
            rCurrentProcessInfo
        );
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);

        for (SizeType i = 0; i < number_of_control_points; i++)
        {
            NodeType & iNode = GetGeometry()[i];

            Vector external_variable = N[i] * rValues[0] + iNode.GetValue(rVariable);
            iNode.SetValue(rVariable, external_variable);
        }

        mConstitutiveLawVector[0]->SetValue(rVariable,
            rValues[0],
            rCurrentProcessInfo
        );
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        mConstitutiveLawVector[0]->SetValue(rVariable,
            rValues[0],
            rCurrentProcessInfo
        );
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::SetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rVariable == CONSTITUTIVE_LAW)
            mConstitutiveLawVector[0] = rValues[0];
    }

    /***********************************************************************************/
    /// GetValueOnIntegrationPoints
    /***********************************************************************************/
    void BaseDiscreteElement::GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == VELOCITY) {
            const int& number_of_control_points = GetGeometry().size();
            Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);

            array_1d<double, 3> velocity = ZeroVector(3);
            for (SizeType i = 0; i < number_of_control_points; i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, 0);

                velocity[0] += N[i] * vel[0];
                velocity[1] += N[i] * vel[1];
                velocity[2] += N[i] * vel[2];
            }
            rValues[0] = velocity;
        }
        else {
            CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        const unsigned int size = GetGeometry().IntegrationPoints().size();

        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == COORDINATES) {
            const int& number_of_control_points = GetGeometry().size();
            Vector N = this->GetValue(SHAPE_FUNCTION_VALUES);
            Vector condition_coords = ZeroVector(3);
            for (SizeType i = 0; i < number_of_control_points; i++)
            {
                const NodeType & iNode = GetGeometry()[i];
                const array_1d<double, 3>& coords = iNode.Coordinates();

                condition_coords[0] += N[i] * coords[0];
                condition_coords[1] += N[i] * coords[1];
                condition_coords[2] += N[i] * coords[2];
            }
            rValues[0] = condition_coords;
        }
        else {
            CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetValuesVector(
        Vector& rValues,
        int Step
    )
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i)
        {
            const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const unsigned int index = i * 3;

            rValues[i]     = displacement[0];
            rValues[i + 1] = displacement[1];
            rValues[i + 2] = displacement[2];
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetFirstDerivativesVector(
        Vector& rValues,
        int Step
    )
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            const unsigned int index = i * 3;

            rValues[index]     = velocity[0];
            rValues[index + 1] = velocity[1];
            rValues[index + 2] = velocity[2];
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void BaseDiscreteElement::GetSecondDerivativesVector(
        Vector& rValues,
        int Step
    )
    {
        const unsigned int number_of_control_points = GetGeometry().size();
        const unsigned int mat_size = number_of_control_points * 3;

        if (rValues.size() != mat_size)
            rValues.resize(mat_size, false);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const unsigned int index = i * 3;

            rValues[index]     = acceleration[0];
            rValues[index + 1] = acceleration[1];
            rValues[index + 2] = acceleration[2];
        }
    }

    //************************************************************************************/
    //************************************************************************************/
    void BaseDiscreteElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_ERROR << "You have called to the CalculateAll() from the base class BaseDiscreteElement" << std::endl;
    }

    //***********************************************************************************/
    //***********************************************************************************/
    void BaseDiscreteElement::Jacobian(const Matrix& DN_De,
        Matrix& Jacobian,
        const int& rWorkingSpaceDimension,
        const int& rLocalSpaceDimension)
    {
        const unsigned int number_of_control_points = GetGeometry().size();

        Jacobian.resize(rWorkingSpaceDimension, rLocalSpaceDimension);

        Jacobian.clear();
        for (unsigned int i = 0; i < number_of_control_points; i++)
        {
            for (unsigned int k = 0; k<rWorkingSpaceDimension; k++)
            {
                for (unsigned int m = 0; m<rLocalSpaceDimension; m++)
                {
                    Jacobian(k, m) += (GetGeometry()[i]).Coordinates()[k] * DN_De(i, m);
                }
            }
        }
    }

/***********************************************************************************/
/***********************************************************************************/
} // Namespace Kratos


