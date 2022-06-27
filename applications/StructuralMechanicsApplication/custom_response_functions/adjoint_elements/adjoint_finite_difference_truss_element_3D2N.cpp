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
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/checks.h"
#include "custom_elements/truss_element_3D2N.hpp"
#include "custom_elements/truss_element_linear_3D2N.hpp"




namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType num_GP = (this->mpPrimalElement->GetGeometry().IntegrationPoints()).size();
    if (rOutput.size() != num_GP) {
        rOutput.resize(num_GP);
    }

    if(rVariable == YOUNG_MODULUS_VARIATIONAL_SENSITIVITY) {
        std::vector<Vector> pseudo_stress;
        this->CalculateOnIntegrationPoints(YOUNG_MODULUS_PSEUDO_PK2_STRESS_VECTOR, pseudo_stress, rCurrentProcessInfo);
        std::vector<Vector> adjoint_strain;
        this->CalculateOnIntegrationPoints(ADJOINT_GL_STRAIN_VECTOR, adjoint_strain, rCurrentProcessInfo);
        const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        rOutput[0] = -1 * pseudo_stress[0][0] * adjoint_strain[0][0] * l_0;
    }
    else if(rVariable == CROSS_AREA_VARIATIONAL_SENSITIVITY) {
        // internal part
        std::vector<Vector> pseudo_stress;
        this->CalculateOnIntegrationPoints(CROSS_AREA_PSEUDO_CAUCHY_STRESS_VECTOR, pseudo_stress, rCurrentProcessInfo);
        std::vector<Vector> adjoint_strain;
        this->CalculateOnIntegrationPoints(ADJOINT_EA_STRAIN_VECTOR, adjoint_strain, rCurrentProcessInfo);
        const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
        rOutput[0] = -1 * pseudo_stress[0][0] * adjoint_strain[0][0] * l;

        // external part
        const Matrix& Ncontainer = this->mpPrimalElement->GetGeometry().ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

        // creating necessary values;
        const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double rho = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);

        double derived_total_mass =  L * rho;
        const SizeType num_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = this->mpPrimalElement->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = num_nodes * dimension;
        Vector body_forces_node = ZeroVector(dimension);
        Vector body_forces_global = ZeroVector(num_dofs);

        // assemble global Vector
        for (IndexType i = 0; i < num_nodes; ++i) {
            body_forces_node =
                derived_total_mass *
                this->mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION) *
                Ncontainer(0, i);

            for (IndexType j = 0; j < dimension; ++j) {
                body_forces_global[(i * dimension) + j] = body_forces_node[j];
            }
        }

        Vector rhs;
        this->mpPrimalElement->CalculateRightHandSide(rhs, rCurrentProcessInfo);

        Vector adjoint_displacements(num_dofs);
        this->GetValuesVector(adjoint_displacements);

        double external_part = 0.00;
        for(IndexType i = 0; i < num_dofs; ++i) {
            external_part += adjoint_displacements[i] * body_forces_global[i];
        }


        rOutput[0] += external_part;
    }
    else if(rVariable == TRUSS_PRESTRESS_PK2_VARIATIONAL_SENSITIVITY) {
        std::vector<Vector> pseudo_stress;
        this->CalculateOnIntegrationPoints(TRUSS_PRESTRESS_PK2_PSEUDO_PK2_STRESS_VECTOR, pseudo_stress, rCurrentProcessInfo);
        std::vector<Vector> adjoint_strain;
        this->CalculateOnIntegrationPoints(ADJOINT_GL_STRAIN_VECTOR, adjoint_strain, rCurrentProcessInfo);
        const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        rOutput[0] = -1 * pseudo_stress[0][0] * adjoint_strain[0][0] * l_0;
    }
    else {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}


template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType num_GP = (this->mpPrimalElement->GetGeometry().IntegrationPoints()).size();
    if (rOutput.size() != num_GP) {
        rOutput.resize(num_GP);
    }
    const SizeType num_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = this->mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs = num_nodes * dimension;

    /*if(rVariable == ADJOINT_GL_STRAIN_VECTOR) {
        const double L = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const array_1d<double, 3> delta_pos =
        this->GetGeometry()[1].GetInitialPosition().Coordinates() -
        this->GetGeometry()[0].GetInitialPosition().Coordinates() +
        this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) -
        this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
        const double l = MathUtils<double>::Norm3(delta_pos);

        KRATOS_ERROR_IF(l <= std::numeric_limits<double>::epsilon())
                << "Element #" << rElement.Id() << " has a current length of zero!" << std::endl;
        const double e = ((l * l - L * L) / (2.00 * L * L));
    }*/
    if(rVariable == ADJOINT_ES_STRAIN_VECTOR || rVariable == ADJOINT_GL_STRAIN_VECTOR || rVariable == ADJOINT_EA_STRAIN_VECTOR) {
        Vector strain = ZeroVector(dimension);

        Vector length_derivative_vector;
        this->CalculateCurrentLengthDisplacementDerivative(length_derivative_vector);

        const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

        double strain_length_derivative = 0.0;
        if( rVariable == ADJOINT_ES_STRAIN_VECTOR) {
            strain_length_derivative = 1.0 / l_0;
        } else if (rVariable == ADJOINT_GL_STRAIN_VECTOR) {
            strain_length_derivative = l / (l_0 * l_0);
        } else {
            strain_length_derivative = 1.0 / l; //(l_0 * l_0) / (l * l * l);
        }

        Vector adjoint_displacements(num_dofs);
        this->GetValuesVector(adjoint_displacements);

        Vector particular_solution = ZeroVector(num_dofs);
        if(this->Has(ADJOINT_PARTICULAR_DISPLACEMENT)) {
            particular_solution = this->GetValue(ADJOINT_PARTICULAR_DISPLACEMENT);
        }

        strain[0] = 0.00;
        for(IndexType i = 0; i < num_dofs; ++i) {
            strain[0] += strain_length_derivative * length_derivative_vector[i] * (adjoint_displacements[i] + particular_solution[i]);
        }

        strain[1] = 0.00;
        strain[2] = 0.00;
        rOutput[0] = strain;
    }
    if(rVariable == YOUNG_MODULUS_PSEUDO_PK2_STRESS_VECTOR) {
        Vector pseudo_stress = ZeroVector(dimension);

        // we assume linear material, i.e., sigma_pk2 = E * epsilon_gl. Hence, the PK2-derivative w.r.t. to E is epsilon_gl.
        std::vector<Vector> gl_strain_vector;
        this->mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, gl_strain_vector, rCurrentProcessInfo);
        pseudo_stress[0] = gl_strain_vector[0][0];

        // we assume pre-integration w.r.t. A in the internal virtual work.
        const double A = this->mpPrimalElement->GetProperties()[CROSS_AREA];
        pseudo_stress[0] *= A;
        pseudo_stress[1] = 0.00;
        pseudo_stress[2] = 0.00;
        rOutput[0] = pseudo_stress;
    }
    if(rVariable == CROSS_AREA_PSEUDO_PK2_STRESS_VECTOR || rVariable == CROSS_AREA_PSEUDO_CAUCHY_STRESS_VECTOR) {
        Vector pseudo_stress = ZeroVector(dimension);
        std::vector<Vector> gl_strain_vector;
        this->mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, gl_strain_vector, rCurrentProcessInfo);
        const double E = this->mpPrimalElement->GetProperties()[YOUNG_MODULUS];
        double prestress = 0.00;
        if (this->mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            prestress = this->mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        // we assume linear material, i.e., sigma_pk2 = E * epsilon_gl + prestress.
        //As we assume pre-integration w.r.t. A in the internal virtual work the pseudo-stress is equal to sigma_pk2.
        pseudo_stress[0] = gl_strain_vector[0][0] * E + prestress;

        pseudo_stress[1] = 0.00;
        pseudo_stress[2] = 0.00;
        rOutput[0] = pseudo_stress;

        if(rVariable == CROSS_AREA_PSEUDO_CAUCHY_STRESS_VECTOR) {
            const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
            const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
            rOutput[0] *= l/L0;
        }
    }
    if(rVariable == TRUSS_PRESTRESS_PK2_PSEUDO_PK2_STRESS_VECTOR) {
        Vector pseudo_stress = ZeroVector(dimension);

        double multiplicator = 0.0;
        if (this->mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            multiplicator = 1.0;
            double prestress = 0.0;
            prestress = this->mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
            if (std::abs(prestress) < 1.0e-12) {
                multiplicator = 0.0;
            }
        }

        const double A = this->mpPrimalElement->GetProperties()[CROSS_AREA];
        pseudo_stress[0] = A * multiplicator;
        pseudo_stress[1] = 0.00;
        pseudo_stress[2] = 0.00;
        rOutput[0] = pseudo_stress;
    }

    KRATOS_CATCH("")
}


template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if(rStressVariable == STRESS_ON_GP)
    {
        const SizeType num_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = this->mpPrimalElement->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs = num_nodes * dimension;
        const SizeType num_GP = (this->mpPrimalElement->GetGeometry().IntegrationPoints()).size();
        if ( (rOutput.size1() != num_dofs) || (rOutput.size2() != num_GP ) )
            rOutput.resize(num_dofs, num_GP);

        double derivative_pre_factor;
        this->GetDerivativePreFactor(derivative_pre_factor, rCurrentProcessInfo);

        Vector length_derivative_vector;
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

template <class TPrimalElement>
int AdjointFiniteDifferenceTrussElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int return_value = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(this->mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!

    KRATOS_ERROR_IF(this->GetGeometry().WorkingSpaceDimension() != 3 || this->GetGeometry().size() != 2)
        << "The truss element works only in 3D and with 2 noded elements" << "" << std::endl;

    this->CheckDofs();
    this->CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this)
         < std::numeric_limits<double>::epsilon())
        << "Element #" << this->Id() << " has a length of zero!" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CheckDofs() const
{
    const SizeType number_of_nodes = this->GetGeometry().size();

    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const auto& r_node = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);
    }
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CheckProperties(const ProcessInfo& rCurrentProcessInfo) const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const PropertiesType & r_properties = this->GetProperties();

    KRATOS_ERROR_IF(r_properties.Has(CROSS_AREA) == false || r_properties[CROSS_AREA] <= numerical_limit)
    << "CROSS_AREA not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF(r_properties.Has(YOUNG_MODULUS) == false || r_properties[YOUNG_MODULUS] <= numerical_limit)
    << "YOUNG_MODULUS not provided for this element" << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT( r_properties.Has(DENSITY) )
    << "DENSITY not provided for this element" << this->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
    << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
    const ConstitutiveLaw::Pointer& cl = r_properties[CONSTITUTIVE_LAW];
    KRATOS_ERROR_IF(cl == nullptr)
    << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
    cl->Check(r_properties ,this->GetGeometry(),rCurrentProcessInfo);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateCurrentLengthDisplacementDerivative(Vector& rDerivativeVector)
{
    KRATOS_TRY;

    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs = number_of_nodes * dimension;

    if (rDerivativeVector.size() != num_dofs)
        rDerivativeVector.resize(num_dofs, false);

    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
    const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
    const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
    const double u1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double u2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    const double v1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double v2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    const double w1 = this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
    const double w2 = this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z);

    rDerivativeVector[0] = (u1 - u2 - dx) / l;
    rDerivativeVector[1] = (v1 - v2 - dy) / l;
    rDerivativeVector[2] = (w1 - w2 - dz) / l;
    rDerivativeVector[3] = -1.0 * rDerivativeVector[0];
    rDerivativeVector[4] = -1.0 * rDerivativeVector[1];
    rDerivativeVector[5] = -1.0 * rDerivativeVector[2];

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::GetDerivativePreFactor(double& rDerivativePreFactor, const ProcessInfo& rCurrentProcessInfo)
{
    TracedStressType traced_stress_type = static_cast<TracedStressType>(this->GetValue(TRACED_STRESS_TYPE));

    switch (traced_stress_type)
    {
        case TracedStressType::FX:
        {
            rDerivativePreFactor = this->CalculateDerivativePreFactorFX(rCurrentProcessInfo);
            break;
        }
        case TracedStressType::PK2:
        {
            rDerivativePreFactor = this->CalculateDerivativePreFactorPK2(rCurrentProcessInfo);
            break;
        }
        default:
            KRATOS_ERROR << "Invalid stress type! Stress type not supported for this element!" << std::endl;
    }
}

template <class TPrimalElement>
double AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateDerivativePreFactorFX(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = this->mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double A = this->mpPrimalElement->GetProperties()[CROSS_AREA];
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    double prestress = 0.0;
    if (this->mpPrimalElement->GetProperties().Has(TRUSS_PRESTRESS_PK2))
        prestress = this->mpPrimalElement->GetProperties()[TRUSS_PRESTRESS_PK2];
    std::vector<Vector> GL_strain;
    this->mpPrimalElement->CalculateOnIntegrationPoints(GREEN_LAGRANGE_STRAIN_VECTOR, GL_strain, rCurrentProcessInfo);
    const double GL_strain_X = GL_strain[0][0]; //one Gauss-Point result is enough due to constant strains.

    double derivative_pre_factor = A / l_0 * (E * GL_strain_X + prestress + E * l * l / (l_0 * l_0));

    KRATOS_DEBUG_ERROR_IF(std::abs(derivative_pre_factor) <= numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

template <class TPrimalElement>
double AdjointFiniteDifferenceTrussElement<TPrimalElement>::CalculateDerivativePreFactorPK2(const ProcessInfo& rCurrentProcessInfo)
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    const double E = this->mpPrimalElement->GetProperties()[YOUNG_MODULUS];
    const double l = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double l_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    double derivative_pre_factor = E * l / (l_0 * l_0);

    KRATOS_DEBUG_ERROR_IF(derivative_pre_factor<=numerical_limit)
        << "Derivative pre-factor of " << this->Id() << "~ 0" << std::endl;

    return derivative_pre_factor;
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceTrussElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

template class AdjointFiniteDifferenceTrussElement<TrussElement3D2N>;
template class AdjointFiniteDifferenceTrussElement<TrussElementLinear3D2N>;

} // namespace Kratos.





