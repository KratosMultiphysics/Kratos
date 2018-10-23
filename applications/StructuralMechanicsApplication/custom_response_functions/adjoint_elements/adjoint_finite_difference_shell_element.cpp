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
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferencingShellElement::AdjointFiniteDifferencingShellElement(Element::Pointer pPrimalElement)
                    : AdjointFiniteDifferencingBaseElement(pPrimalElement, true) {}

AdjointFiniteDifferencingShellElement::~AdjointFiniteDifferencingShellElement() {}

int AdjointFiniteDifferencingShellElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int return_value = AdjointFiniteDifferencingBaseElement::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!
    CheckVariables();
    CheckDofs();
    CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(GetGeometry().Area() < std::numeric_limits<double>::epsilon()*1000)
        << "Element #" << Id() << " has an Area of zero!" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

// private

void AdjointFiniteDifferencingShellElement::CheckVariables()
{
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ROTATION);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
    KRATOS_CHECK_VARIABLE_KEY(SHELL_ORTHOTROPIC_LAYERS);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);
    KRATOS_CHECK_VARIABLE_KEY(SHELL_CROSS_SECTION);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_ROTATION);
}

void AdjointFiniteDifferencingShellElement::CheckDofs()
{
    GeometryType& r_geom = GetGeometry();
    // verify that the dofs exist
    for (IndexType i = 0; i < r_geom.size(); ++i)
    {
        auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_ROTATION_Z, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);

        KRATOS_ERROR_IF(r_node.GetBufferSize() < 2) << "This Element needs "
            << "at least a buffer size = 2" << std::endl;
    }
}

void AdjointFiniteDifferencingShellElement::CheckProperties(const ProcessInfo& rCurrentProcessInfo)
{
    // check properties
    if(pGetProperties() == nullptr)
        KRATOS_ERROR << "Properties not provided for element " << Id() << std::endl;

    const PropertiesType & props = GetProperties();

    const GeometryType& geom = GetGeometry(); // TODO check if this can be const

    if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
    {
        const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
        if(section == nullptr)
            KRATOS_ERROR << "SHELL_CROSS_SECTION not provided for element " << Id() << std::endl;

        section->Check(props, geom, rCurrentProcessInfo);
    }
    else if (props.Has(SHELL_ORTHOTROPIC_LAYERS))
    {
        CheckSpecificProperties();

        // perform detailed orthotropic check later in shell_cross_section
    }
    else // ... allow the automatic creation of a homogeneous section from a material and a thickness
    {
        CheckSpecificProperties();

        // TODO is this needed???? => it is, the dummy is needed for "Check" => unify!
        ShellCrossSection::Pointer dummySection = Kratos::make_shared<ShellCrossSection>(ShellCrossSection());
        dummySection->BeginStack();
        dummySection->AddPly(0, 5, GetProperties());
        dummySection->EndStack();
        dummySection->SetSectionBehavior(ShellCrossSection::Thick);
        dummySection->Check(props, geom, rCurrentProcessInfo);
    }

}

void AdjointFiniteDifferencingShellElement::CheckSpecificProperties()
{
    const PropertiesType & r_props = GetProperties();

    if (!r_props.Has(CONSTITUTIVE_LAW))
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << Id() << std::endl;
    const ConstitutiveLaw::Pointer& claw = r_props[CONSTITUTIVE_LAW];
    if (claw == nullptr)
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << Id() << std::endl;

    if(!r_props.Has(THICKNESS))
        KRATOS_ERROR << "THICKNESS not provided for element " << Id() << std::endl;
    if(r_props[THICKNESS] <= 0.0)
        KRATOS_ERROR << "wrong THICKNESS value provided for element " << Id() << std::endl;

    if(!r_props.Has(DENSITY))
        KRATOS_ERROR << "DENSITY not provided for element " << Id() << std::endl;
    if(r_props[DENSITY] < 0.0)
        KRATOS_ERROR << "wrong DENSITY value provided for element " << Id() << std::endl;

    // TODO for thick shell Stenberg stabilization is not checked
}


double AdjointFiniteDifferencingShellElement::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable)
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

