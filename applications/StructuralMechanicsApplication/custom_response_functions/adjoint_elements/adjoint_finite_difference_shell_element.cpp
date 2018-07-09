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
#include "adjoint_finite_difference_shell_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/response_data.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferencingShellElement::AdjointFiniteDifferencingShellElement(Element::Pointer pPrimalElement)
                    : AdjointFiniteDifferencingBaseElement(pPrimalElement) {}

AdjointFiniteDifferencingShellElement::~AdjointFiniteDifferencingShellElement() {}

Element::Pointer AdjointFiniteDifferencingShellElement::Create(Element::Pointer pPrimalElement) const
{
    return Kratos::make_shared<AdjointFiniteDifferencingShellElement>(pPrimalElement);
}

void AdjointFiniteDifferencingShellElement::EquationIdVector(EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    const SizeType num_dofs = 6 * GetGeometry().PointsNumber();
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

void AdjointFiniteDifferencingShellElement::GetDofList(DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo) {

    const SizeType num_dofs = 6 * GetGeometry().PointsNumber();
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

void AdjointFiniteDifferencingShellElement::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY
    const SizeType num_dofs = 6 * GetGeometry().PointsNumber();
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
    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingShellElement::Calculate(const Variable<Vector >& rVariable,
                           Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_gps = GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());

    if(rVariable == STRESS_ON_GP)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

        int direction_1 = 0;
        int direction_2 = 0;
        std::vector<Matrix> stress_vector;
        bool stress_is_moment = true;

        switch (traced_stress_type)
        {
            case TracedStressType::MXX:
            {
                direction_1 = 0;
                direction_2 = 0;
                break;
            }
            case TracedStressType::MXY:
            {
                direction_1 = 0;
                direction_2 = 1;
                break;
            }
            case TracedStressType::MXZ:
            {
                direction_1 = 0;
                direction_2 = 2;
                break;
            }
            case TracedStressType::MYX:
            {
                direction_1 = 1;
                direction_2 = 0;
                break;
            }
            case TracedStressType::MYY :
            {
                direction_1 = 1;
                direction_2 = 1;
                break;
            }
            case TracedStressType::MYZ:
            {
                direction_1 = 1;
                direction_2 = 2;
                break;
            }
            case TracedStressType::MZX:
            {
                direction_1 = 2;
                direction_2 = 0;
                break;
            }
            case TracedStressType::MZY:
            {
                direction_1 = 2;
                direction_2 = 1;
                break;
            }
            case TracedStressType::MZZ :
            {
                direction_1 = 2;
                direction_2 = 2;
                break;
            }
            case TracedStressType::FXX :
            {
                direction_1 = 0;
                direction_2 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FXY:
            {
                direction_1 = 0;
                direction_2 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FXZ:
            {
                direction_1 = 0;
                direction_2 = 2;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FYX:
            {
                direction_1 = 1;
                direction_2 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FYY:
            {
                direction_1 = 1;
                direction_2 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FYZ:
            {
                direction_1 = 1;
                direction_2 = 2;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZX:
            {
                direction_1 = 2;
                direction_2 = 0;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZY:
            {
                direction_1 = 2;
                direction_2 = 1;
                stress_is_moment = false;
                break;
            }
            case TracedStressType::FZZ:
            {
                direction_1 = 2;
                direction_2 = 2;
                stress_is_moment = false;
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }

        if(stress_is_moment)
            mpPrimalElement->GetValueOnIntegrationPoints(SHELL_MOMENT_GLOBAL, stress_vector, rCurrentProcessInfo);
        else
            mpPrimalElement->GetValueOnIntegrationPoints(SHELL_FORCE_GLOBAL, stress_vector, rCurrentProcessInfo);

        rOutput.resize(num_gps);
        for(size_t i = 0; i < num_gps; i++)
        {
            rOutput(i) = stress_vector[i](direction_1, direction_2);
        }

    }
    else
    {
        rOutput.resize(num_gps);
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

int AdjointFiniteDifferencingShellElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    GeometryType& r_geom = GetGeometry();

    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(SHELL_CROSS_SECTION);
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);

    // check properties
    KRATOS_ERROR_IF(this->pGetProperties() == nullptr) << "Properties not provided for element " << this->Id() << std::endl;

    const PropertiesType & props = this->GetProperties();

    if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
    {
        const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
        KRATOS_ERROR_IF(section == nullptr) << "SHELL_CROSS_SECTION not provided for element " << this->Id() << std::endl;

        section->Check(props, r_geom, rCurrentProcessInfo);
    }
    else // ... allow the automatic creation of a homogeneous section from a material and a thickness
    {
        KRATOS_ERROR_IF_NOT(props.Has(CONSTITUTIVE_LAW)) << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
        const ConstitutiveLaw::Pointer& claw = props[CONSTITUTIVE_LAW];
        KRATOS_ERROR_IF(claw == nullptr) << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;

        KRATOS_ERROR_IF_NOT(props.Has(THICKNESS)) <<  "THICKNESS not provided for element " <<  this->Id() << std::endl;
        KRATOS_ERROR_IF(props[THICKNESS] <= 0.0) << "wrong THICKNESS value provided for element " << this->Id() << std::endl;

        ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
        dummySection->BeginStack();
        dummySection->AddPly(0.0, 5, props);
        dummySection->EndStack();
        dummySection->SetSectionBehavior(ShellCrossSection::Thin);
        dummySection->Check(props, r_geom, rCurrentProcessInfo);
    }

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

// private

double AdjointFiniteDifferencingShellElement::GetPerturbationSizeCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        double dx, dy, dz, L = 0.0;

        const GeometryType& geometry = mpPrimalElement->GetGeometry();

        dx = geometry[1].X0() - geometry[0].X0();
        dy = geometry[1].Y0() - geometry[0].Y0();
        dz = geometry[1].Z0() - geometry[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = geometry[2].X0() - geometry[1].X0();
        dy = geometry[2].Y0() - geometry[1].Y0();
        dz = geometry[2].Z0() - geometry[1].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = geometry[2].X0() - geometry[0].X0();
        dy = geometry[2].Y0() - geometry[0].Y0();
        dz = geometry[2].Z0() - geometry[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        L /= 3.0;

        return L;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferencingShellElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement );
}

void AdjointFiniteDifferencingShellElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, AdjointFiniteDifferencingBaseElement );

}

} // namespace Kratos

