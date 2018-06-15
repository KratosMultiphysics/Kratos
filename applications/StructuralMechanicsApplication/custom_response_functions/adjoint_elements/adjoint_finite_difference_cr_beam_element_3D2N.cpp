// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#include "adjoint_finite_difference_cr_beam_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/response_data.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferenceCrBeamElement::AdjointFiniteDifferenceCrBeamElement(IndexType NewId,
    GeometryType::Pointer pGeometry)
    : AdjointFiniteDifferencingBaseElement(NewId, pGeometry)
{
}

AdjointFiniteDifferenceCrBeamElement::AdjointFiniteDifferenceCrBeamElement(IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties, Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(NewId, pGeometry, pProperties, pPrimalElement)
{
    mpPrimalBeamElement = dynamic_pointer_cast<CrBeamElement3D2N>(pPrimalElement);
}

AdjointFiniteDifferenceCrBeamElement::~AdjointFiniteDifferenceCrBeamElement()
{
}

Element::Pointer AdjointFiniteDifferenceCrBeamElement::Create(IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties, Element::Pointer pPrimalElement) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_shared<AdjointFiniteDifferenceCrBeamElement>(
        NewId, rGeom.Create(rThisNodes), pProperties, pPrimalElement);
}

void AdjointFiniteDifferenceCrBeamElement::EquationIdVector(EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo) {

    const int number_of_nodes = this->GetGeometry().PointsNumber();
    const int dimension = this->GetGeometry().WorkingSpaceDimension();
    const unsigned int local_size = number_of_nodes * dimension * 2;

    if (rResult.size() != local_size) rResult.resize(local_size);

    GeometryType & geom = this->GetGeometry();

    for (int i = 0; i < number_of_nodes; ++i)
    {
        int index = i * number_of_nodes * dimension;
        rResult[index]     = geom[i].GetDof(ADJOINT_DISPLACEMENT_X).EquationId();
        rResult[index + 1] = geom[i].GetDof(ADJOINT_DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = geom[i].GetDof(ADJOINT_DISPLACEMENT_Z).EquationId();

        rResult[index + 3] = geom[i].GetDof(ADJOINT_ROTATION_X).EquationId();
        rResult[index + 4] = geom[i].GetDof(ADJOINT_ROTATION_Y).EquationId();
        rResult[index + 5] = geom[i].GetDof(ADJOINT_ROTATION_Z).EquationId();
    }

}

void AdjointFiniteDifferenceCrBeamElement::GetDofList(DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo) {

    GeometryType & geom = this->GetGeometry();

    const int number_of_nodes = geom.PointsNumber();
    const int dimension = geom.WorkingSpaceDimension();
    const unsigned int local_size = number_of_nodes * dimension * 2;

    if (rElementalDofList.size() != local_size) {
        rElementalDofList.resize(local_size);
    }

    for (int i = 0; i < number_of_nodes; ++i)
    {
        int index = i * number_of_nodes * dimension;
        rElementalDofList[index]     = geom[i].pGetDof(ADJOINT_DISPLACEMENT_X);
        rElementalDofList[index + 1] = geom[i].pGetDof(ADJOINT_DISPLACEMENT_Y);
        rElementalDofList[index + 2] = geom[i].pGetDof(ADJOINT_DISPLACEMENT_Z);

        rElementalDofList[index + 3] = geom[i].pGetDof(ADJOINT_ROTATION_X);
        rElementalDofList[index + 4] = geom[i].pGetDof(ADJOINT_ROTATION_Y);
        rElementalDofList[index + 5] = geom[i].pGetDof(ADJOINT_ROTATION_Z);
    }
}

void AdjointFiniteDifferenceCrBeamElement::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY
    const GeometryType & geom = this->GetGeometry();
    const int number_of_nodes = geom.PointsNumber();
    const int dimension = geom.WorkingSpaceDimension();
    const unsigned int num_dofs = number_of_nodes * dimension * 2;

    if (rValues.size() != num_dofs)
        rValues.resize(num_dofs, false);

    for (int i = 0; i < number_of_nodes; ++i)
    {
        const SizeType index = i * dimension * 2;

        const array_1d<double,3>& disp = geom[i].FastGetSolutionStepValue(ADJOINT_DISPLACEMENT, Step);
        const array_1d<double,3>& rot = geom[i].FastGetSolutionStepValue(ADJOINT_ROTATION, Step);

        rValues[index]     = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];

        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
        rValues[index + 5] = rot[2];
    }
    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceCrBeamElement::Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_WATCH(rVariable)
    KRATOS_WATCH(mpPrimalElement)
    KRATOS_WATCH(this)

    if(rVariable == STRESS_ON_GP || rVariable == STRESS_ON_NODE)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

        KRATOS_WATCH(this->GetValue(TRACED_STRESS_TYPE))

        std::vector< array_1d<double, 3 > > stress_vector;
        int direction_1 = 0;
        bool stress_is_moment = true;

        switch (traced_stress_type)
        {
            case TracedStressType::MX:
            {
                direction_1 = 0;
                break;
            }
            case TracedStressType::MY:
            {
                direction_1 = 1;
                break;
            }
            case TracedStressType::MZ:
            {
                direction_1 = 2;
                break;
            }
            case TracedStressType::FX:
            {
                direction_1 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FY:
            {
                direction_1 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZ:
            {
                direction_1 = 2;
                stress_is_moment = false;
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }

        if(stress_is_moment)
            mpPrimalElement->GetValueOnIntegrationPoints(MOMENT, stress_vector, rCurrentProcessInfo);
        else
            mpPrimalElement->GetValueOnIntegrationPoints(FORCE, stress_vector, rCurrentProcessInfo);

        if(rVariable == STRESS_ON_GP)
        {
            const unsigned int&  GP_num = GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);

            rOutput.resize(GP_num);
            for(unsigned int i = 0; i < GP_num ; i++)
            {
                rOutput(i) = stress_vector[i][direction_1];
            }
        }
        else if(rVariable == STRESS_ON_NODE)
        {
            rOutput.resize(2);
            rOutput(0) = 2 * stress_vector[0][direction_1] - stress_vector[1][direction_1];
            rOutput(1) = 2 * stress_vector[2][direction_1] - stress_vector[1][direction_1];
        }
    }
    else
    {
        rOutput.resize(3);
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceCrBeamElement::Calculate(const Variable<Matrix >& rVariable,
                        Matrix& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rVariable == STRESS_DISP_DERIV_ON_GP)
    {
            this->CalculateStressDisplacementDerivative(STRESS_ON_GP, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DISP_DERIV_ON_NODE)
    {
        this->CalculateStressDisplacementDerivative(STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
    {
        const std::string design_varible_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_GP, rOutput, rCurrentProcessInfo);
        }
    }
    else if(rVariable == STRESS_DESIGN_DERIVATIVE_ON_NODE)
    {
        std::string design_varible_name = this->GetValue( DESIGN_VARIABLE_NAME );

        if (KratosComponents<Variable<double>>::Has(design_varible_name) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(design_varible_name) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(design_varible_name);
            this->CalculateStressDesignVariableDerivative(r_variable, STRESS_ON_NODE, rOutput, rCurrentProcessInfo);
        }
    }
        else
    {
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceCrBeamElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != 3 || GetGeometry().size() != 2)
    << "The beam element works only in 3D and with 2 noded elements" << "" << std::endl;

    mpPrimalBeamElement->CheckVariables();
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    mpPrimalBeamElement->CheckProperties();

    // check dofs
    GeometryType& r_geom = GetGeometry();
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

double AdjointFiniteDifferenceCrBeamElement::GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
        double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
        double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
        double L = std::sqrt(dx*dx + dy*dy + dz*dz);
        return L;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceCrBeamElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceCrBeamElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


