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


#include "adjoint_finite_difference_truss_element_3D2N.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"


namespace Kratos
{

AdjointFiniteDifferenceTrussElement::AdjointFiniteDifferenceTrussElement(Element::Pointer pPrimalElement)
    : AdjointFiniteDifferencingBaseElement(pPrimalElement)
{
}

AdjointFiniteDifferenceTrussElement::~AdjointFiniteDifferenceTrussElement()
{
}

void AdjointFiniteDifferenceTrussElement::Calculate(const Variable<Vector >& rVariable,
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

void AdjointFiniteDifferenceTrussElement::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rStressVariable == STRESS_ON_GP)
    {
        TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs_per_node = dimension; 
        const SizeType num_dofs = num_nodes * num_dofs_per_node;
        const SizeType num_GP = (mpPrimalElement->GetGeometry().IntegrationPoints()).size();
        Vector length_derivative_vector;
        double derivative_pre_factor = 0.0;
   
        switch (traced_stress_type)
        {
            case TracedStressType::FX:
            {
                const double E = mpPrimalElement->GetProperties()[YOUNG_MODULUS];
                const double A = mpPrimalElement->GetProperties()[CROSS_AREA];
                const double L0 = CalculateReferenceLength();
                const double l = CalculateCurrentLength();
                double prestress = 0.00;
                if (mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2))
                    prestress = mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
                std::vector<Vector> GL_strain;  
                mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, GL_strain, rCurrentProcessInfo);
                const double GL_strain_X = GL_strain[0][0]; //one Gauss-Point result is enough due to constant strains.
        
                derivative_pre_factor = A / L0 * (E * GL_strain_X + prestress + E * l * l / (L0 * L0));
                break;
            }
            case TracedStressType::PK2X:
            {
                const double E = mpPrimalElement->GetProperties()[YOUNG_MODULUS];
                const double l = CalculateCurrentLength();
                const double L0 = CalculateReferenceLength();
                derivative_pre_factor = E * l / (L0 * L0);
                break;
            }
            default:
                KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
        }
        KRATOS_ERROR_IF(std::abs(derivative_pre_factor) <  numerical_limit)
            << "pre factor of stress displacement derivative is ~0!" << std::endl;

        rOutput.resize(num_dofs, num_GP);
        rOutput.clear();
        this->CalculateCurrentLengthDisplacementDerivative(length_derivative_vector);
        for(IndexType i = 0; i < num_dofs; ++i)
        {
             for(IndexType j = 0; j < num_GP; ++j)
                rOutput(i, j) = length_derivative_vector[i] * derivative_pre_factor; 
        }
    }
    else
        KRATOS_ERROR << "Stress displacement derivative only available for Gauss-points quantities!" << std::endl;

    KRATOS_CATCH("")
}

int AdjointFiniteDifferenceTrussElement::Check(const ProcessInfo& rCurrentProcessInfo)
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

double AdjointFiniteDifferenceTrussElement::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable)
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE)
        return this->CalculateReferenceLength();
    else
        return 1.0;

    KRATOS_CATCH("")
}

double AdjointFiniteDifferenceTrussElement::CalculateReferenceLength() 
{
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double L = std::sqrt(dx * dx + dy * dy + dz * dz);

    KRATOS_ERROR_IF(L<=numerical_limit)
     << "Reference Length of element" << this->Id() << "~ 0" << std::endl;
    return L;
    KRATOS_CATCH("")
}

double AdjointFiniteDifferenceTrussElement::CalculateCurrentLength() 
{
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double du =
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double dv =
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double dw =
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                               (dw + dz) * (dw + dz));

    KRATOS_ERROR_IF(l<=numerical_limit)
     << "Current Length of element" << this->Id() << "~ 0" << std::endl;
    return l;
    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTrussElement::CalculateCurrentLengthDisplacementDerivative(Vector& rDerivativeVector)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    KRATOS_ERROR_IF_NOT(num_dofs == 6) << "This implementation expects 6 dofs!" << std::endl;

    if (rDerivativeVector.size() != num_dofs)
        rDerivativeVector.resize(num_dofs, false);

    const double l = CalculateCurrentLength();
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double u1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X); 
    const double u2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double v1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y); 
    const double v2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double w1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z); 
    const double w2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z);

    rDerivativeVector[0] = -1.0 * (u2 - u1 + dx) / l;
    rDerivativeVector[1] = -1.0 * (v2 - v1 + dy) / l;
    rDerivativeVector[2] = -1.0 * (w2 - w1 + dz) / l;
    rDerivativeVector[3] = -1.0 * rDerivativeVector[0];
    rDerivativeVector[4] = -1.0 * rDerivativeVector[1];
    rDerivativeVector[5] = -1.0 * rDerivativeVector[2];
   
    KRATOS_CATCH("")
}

void AdjointFiniteDifferenceTrussElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

void AdjointFiniteDifferenceTrussElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, AdjointFiniteDifferencingBaseElement);
}

} // namespace Kratos.


