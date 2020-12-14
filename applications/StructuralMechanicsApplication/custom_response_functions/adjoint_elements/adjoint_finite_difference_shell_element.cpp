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
#include "custom_elements/shell_thin_element_3D3N.hpp"


namespace Kratos
{

template <class TPrimalElement>
int AdjointFiniteDifferencingShellElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int return_value = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(this->mHasRotationDofs) << "Adjoint shell element does not have rotation dofs!" << std::endl;
    KRATOS_ERROR_IF_NOT(this->mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!
    this->CheckDofs();
    this->CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(this->GetGeometry().Area() < std::numeric_limits<double>::epsilon()*1000)
        << "Element #" << this->Id() << " has an Area of zero!" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

// private

template <class TPrimalElement>
void AdjointFiniteDifferencingShellElement<TPrimalElement>::CheckDofs() const
{
    const GeometryType& r_geom = this->GetGeometry();
    // verify that the dofs exist
    for (IndexType i = 0; i < r_geom.size(); ++i)
    {
        const auto& r_node = r_geom[i];
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

template <class TPrimalElement>
void AdjointFiniteDifferencingShellElement<TPrimalElement>::CheckProperties(const ProcessInfo& rCurrentProcessInfo) const
{
    // check properties
    if(this->pGetProperties() == nullptr)
        KRATOS_ERROR << "Properties not provided for element " << this->Id() << std::endl;

    const PropertiesType & props = this->GetProperties();

    const GeometryType& geom = this->GetGeometry(); // TODO check if this can be const

    if (props.Has(SHELL_ORTHOTROPIC_LAYERS))
    {
        this->CheckSpecificProperties();

        // perform detailed orthotropic check later in shell_cross_section
    }
    else // ... allow the automatic creation of a homogeneous section from a material and a thickness
    {
        this->CheckSpecificProperties();

        // TODO is this needed???? => it is, the dummy is needed for "Check" => unify!
        ShellCrossSection::Pointer dummySection = Kratos::make_shared<ShellCrossSection>(ShellCrossSection());
        dummySection->BeginStack();
        dummySection->AddPly(0, 5, this->GetProperties());
        dummySection->EndStack();
        dummySection->SetSectionBehavior(ShellCrossSection::Thick);
        dummySection->Check(props, geom, rCurrentProcessInfo);
    }

}

template <class TPrimalElement>
void AdjointFiniteDifferencingShellElement<TPrimalElement>::CheckSpecificProperties() const
{
    const PropertiesType & r_props = this->GetProperties();

    if (!r_props.Has(CONSTITUTIVE_LAW))
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
    const ConstitutiveLaw::Pointer& claw = r_props[CONSTITUTIVE_LAW];
    if (claw == nullptr)
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;

    if(!r_props.Has(THICKNESS))
        KRATOS_ERROR << "THICKNESS not provided for element " << this->Id() << std::endl;
    if(r_props[THICKNESS] <= 0.0)
        KRATOS_ERROR << "wrong THICKNESS value provided for element " << this->Id() << std::endl;

    if(!r_props.Has(DENSITY))
        KRATOS_ERROR << "DENSITY not provided for element " << this->Id() << std::endl;
    if(r_props[DENSITY] < 0.0)
        KRATOS_ERROR << "wrong DENSITY value provided for element " << this->Id() << std::endl;

    // TODO for thick shell Stenberg stabilization is not checked
}


template <class TPrimalElement>
double AdjointFiniteDifferencingShellElement<TPrimalElement>::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE_SENSITIVITY)
    {
        double dx, dy, dz, L = 0.0;

        const GeometryType& geometry = this->mpPrimalElement->GetGeometry();

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

template <class TPrimalElement>
void AdjointFiniteDifferencingShellElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferencingShellElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);

}

template class AdjointFiniteDifferencingShellElement<ShellThinElement3D3N>;

} // namespace Kratos

