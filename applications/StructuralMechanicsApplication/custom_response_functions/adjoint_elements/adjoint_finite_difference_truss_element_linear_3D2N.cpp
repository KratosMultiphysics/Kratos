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


#include "adjoint_finite_difference_truss_element_linear_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferenceTrussElementLinear::AdjointFiniteDifferenceTrussElementLinear(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(pPrimalElement)
{
}

AdjointFiniteDifferenceTrussElementLinear::~AdjointFiniteDifferenceTrussElementLinear()
{
}

void AdjointFiniteDifferenceTrussElementLinear::Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    ProcessInfo copy_process_info = rCurrentProcessInfo;
    Vector dummy_RHS;
    // This call of RHS is necessary in order to ensure that constitutive law saves the primal stresses.
    // This is possible since within the computation of RHS CalculateMaterialResponse of the CL is called.
    mpPrimalElement->CalculateRightHandSide(dummy_RHS, copy_process_info);

    if(rVariable == STRESS_ON_GP)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

        const SizeType  GP_num = (mpPrimalElement->GetGeometry().IntegrationPoints()).size();
        if (rOutput.size() != GP_num) 
            rOutput.resize(GP_num, false);

        switch (traced_stress_type)
        {
            case TracedStressType::FX:
            {
                std::vector< array_1d<double, 3 > > force_vector;
                mpPrimalElement->GetValueOnIntegrationPoints(FORCE, force_vector, rCurrentProcessInfo);
                for(IndexType i = 0; i < GP_num ; ++i)
                    rOutput(i) = force_vector[i][0];
                break;
            }
            case TracedStressType::PK2X:
            {
                std::vector<Vector> stress_vector;
                mpPrimalElement->GetValueOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo);
                for(IndexType i = 0; i < GP_num ; ++i)
                    rOutput(i) = stress_vector[i][0];
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }
    }
    else
    {
        rOutput.resize(1);
        rOutput.clear();
    }

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceTrussElementLinear::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    KRATOS_ERROR_IF(dimension != 3 || number_of_nodes != 2)
    << "The beam element works only in 3D and with 2 noded elements" << "" << std::endl;

    // verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(CROSS_AREA);
    KRATOS_CHECK_VARIABLE_KEY(ADJOINT_DISPLACEMENT);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (IndexType i = 0; i < number_of_nodes; ++i) 
    {
        NodeType &r_node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
    }

    KRATOS_ERROR_IF(this->GetProperties().Has(CROSS_AREA) == false || this->GetProperties()[CROSS_AREA] <= numerical_limit)
    << "CROSS_AREA not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF(this->GetProperties().Has(YOUNG_MODULUS) == false || this->GetProperties()[YOUNG_MODULUS] <= numerical_limit)
    << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( this->GetProperties().Has(DENSITY) )
    << "DENSITY not provided for this element" << this->Id() << std::endl;

    // TODO: how to check the constitutive law?
    // In the case of the truss only YOUNG_MODULUS and DENSITY are checked (Key and meaningful value)
    // Note: if constitutivelaw is a nullptr is checked in intialize of truss element.
    //if(this->mpConstitutiveLaw != nullptr) 
    //  this->mpConstitutiveLaw->Check(this->GetProperties(),this->GetGeometry(),rCurrentProcessInfo);
    
    return 0;

    KRATOS_CATCH("")
}

double AdjointFiniteDifferenceTrussElementLinear::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
        const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
        const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
        const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

        KRATOS_ERROR_IF(L<=numerical_limit)
            << "Pertubation size for shape derivatives of element" << this->Id() << "~ 0" << std::endl;
        return L;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTrussElementLinear::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceTrussElementLinear::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


