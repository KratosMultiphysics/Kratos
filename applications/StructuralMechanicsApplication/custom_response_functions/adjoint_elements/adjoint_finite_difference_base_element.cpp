// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

// System includes

// External includes

// Project includes
#include "adjoint_finite_difference_base_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/response_data.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferencingBaseElement::AdjointFiniteDifferencingBaseElement(IndexType NewId,
                        GeometryType::Pointer pGeometry)
                        : Element(NewId, pGeometry) {}

AdjointFiniteDifferencingBaseElement::AdjointFiniteDifferencingBaseElement(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties, Element::Pointer pPrimalElement)
                    : Element(NewId, pGeometry, pProperties) {

    mpPrimalElement = pPrimalElement;
    }

AdjointFiniteDifferencingBaseElement::~AdjointFiniteDifferencingBaseElement() {}

Element::Pointer AdjointFiniteDifferencingBaseElement::Create(IndexType NewId,
                NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties,
                Element::Pointer pPrimalElement) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_shared<AdjointFiniteDifferencingBaseElement>(
        NewId, rGeom.Create(rThisNodes), pProperties, pPrimalElement);
}

void AdjointFiniteDifferencingBaseElement::EquationIdVector(EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node
    if(rResult.size() != num_dofs)
        rResult.resize(num_dofs, false);

    GeometryType & geom = this->GetGeometry();

    for(SizeType i = 0; i < geom.size(); i++)
    {
        int index = i * 6;
        NodeType & iNode = geom[i];

        rResult[index]     = iNode.GetDof(ADJOINT_DISPLACEMENT_X).EquationId();
        rResult[index + 1] = iNode.GetDof(ADJOINT_DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = iNode.GetDof(ADJOINT_DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = iNode.GetDof(ADJOINT_ROTATION_X).EquationId();
        rResult[index + 4] = iNode.GetDof(ADJOINT_ROTATION_Y).EquationId();
        rResult[index + 5] = iNode.GetDof(ADJOINT_ROTATION_Z).EquationId();
    }

}

void AdjointFiniteDifferencingBaseElement::GetDofList(DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo) {

    const SizeType num_dofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node
    rElementalDofList.resize(0);
    rElementalDofList.reserve(num_dofs);

    GeometryType & geom = this->GetGeometry();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        NodeType & iNode = geom[i];

        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_X));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Y));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_DISPLACEMENT_Z));

        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_X));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Y));
        rElementalDofList.push_back(iNode.pGetDof(ADJOINT_ROTATION_Z));
    }
}

void AdjointFiniteDifferencingBaseElement::GetValuesVector(Vector& rValues, int Step) {

    const SizeType num_dofs = 6 * GetGeometry().PointsNumber(); // 6 dofs per node
    if(rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    const GeometryType & geom = this->GetGeometry();
    const int dimension = geom.WorkingSpaceDimension();

    for (SizeType i = 0; i < geom.size(); i++)
    {
        const NodeType & iNode = geom[i];
        const array_1d<double,3>& disp = iNode.FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const array_1d<double,3>& rot = iNode.FastGetSolutionStepValue(ADJOINT_ROTATION, Step);

        const SizeType index = i * dimension * 2;
        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
}

void AdjointFiniteDifferencingBaseElement::Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rVariable == STRESS_DISP_DERIV_ON_GP)
    {
       this->CalculateStressDisplacementDerivative(STRESS_ON_GP, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
    {
        std::string design_variable_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_variable_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_variable_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_variable_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
    }
    else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                          std::vector<double>& rOutput,
                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")

}

void AdjointFiniteDifferencingBaseElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                    std::vector<double>& rValues,
                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    KRATOS_CATCH("")
}

int AdjointFiniteDifferencingBaseElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geom = GetGeometry();

    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    // TODO generic way of doing these checks without checking the dofs..

    // Check dofs
    for (unsigned int i = 0; i < r_geom.size(); i++)
    {
        auto& r_node = r_geom[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);
    }

    return 0;

    KRATOS_CATCH("")
}

// Sensitivity functions

void AdjointFiniteDifferencingBaseElement::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")

}

void AdjointFiniteDifferencingBaseElement::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")

}
void AdjointFiniteDifferencingBaseElement::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable,
                                                const Variable<Vector>& rStressVariable, Matrix& rOutput,
                                                const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                            Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;

    KRATOS_CATCH("")
}

// private

double AdjointFiniteDifferencingBaseElement::GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rDesignVariable)
{
    KRATOS_TRY;

    if ( mpPrimalElement->GetProperties().Has(rDesignVariable) )
    {
        const double variable_value = mpPrimalElement->GetProperties()[rDesignVariable];
        return variable_value;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

double AdjointFiniteDifferencingBaseElement::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        KRATOS_ERROR << "NOT_IMPLEMENTED" << std::endl;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingBaseElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer,  Element );
}

void AdjointFiniteDifferencingBaseElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );

}

} // namespace Kratos

