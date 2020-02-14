//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#include "romfemdem3d_element.hpp"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

RomFemDem3DElement::RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : FemDem3DElement(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RomFemDem3DElement::RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : FemDem3DElement(NewId, pGeometry, pProperties)
{
    //BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

RomFemDem3DElement::RomFemDem3DElement(RomFemDem3DElement const &rOther)
    : FemDem3DElement(rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
    //IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

RomFemDem3DElement &RomFemDem3DElement::operator=(RomFemDem3DElement const &rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    FemDem3DElement::operator=(rOther);
    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer RomFemDem3DElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Element::Pointer(new RomFemDem3DElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer RomFemDem3DElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    RomFemDem3DElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

    return Element::Pointer(new RomFemDem3DElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

RomFemDem3DElement::~RomFemDem3DElement()
{
}

void RomFemDem3DElement::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
    //*****************************
    KRATOS_TRY

    //1.-Initialize sizes for the system components:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

    Vector strain_vector(voigt_size);
    noalias(strain_vector) = ZeroVector(voigt_size);
    Vector stress_vector(voigt_size);
    noalias(stress_vector) = ZeroVector(voigt_size);
    Matrix B(voigt_size, dimension * number_of_nodes);
    noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
    Matrix DN_DX(number_of_nodes, dimension);
    noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

    //deffault values for the infinitessimal theory
    double detF = 1;
    Matrix F(dimension, dimension);
    noalias(F) = identity_matrix<double>(dimension);

    //3.-Calculate elemental system:

    //reading integration points
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    //get the shape functions [N] (for the order of the default integration method)
    const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

    //get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
    const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    //calculate delta position (here coincides with the current displacement)
    Matrix delta_position(number_of_nodes, dimension);
    noalias(delta_position) = ZeroMatrix(number_of_nodes, dimension);
    delta_position = this->CalculateDeltaPosition(delta_position);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
    GeometryType::JacobiansType J;
    J.resize(1, false);
    J[0].resize(dimension, dimension, false);
    noalias(J[0]) = ZeroMatrix(dimension, dimension);
    J = GetGeometry().Jacobian(J, mThisIntegrationMethod, delta_position);

    // Loop Over Integration Points
    for (unsigned int integration_point = 0; integration_point < integration_points.size(); integration_point++) {
        Matrix InvJ(dimension, dimension);
        noalias(InvJ) = ZeroMatrix(dimension, dimension);
        double detJ = 0;
        MathUtils<double>::InvertMatrix(J[integration_point], InvJ, detJ);
        KRATOS_ERROR_IF(detJ < 0) << "SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0" << std::endl;

        //compute cartesian derivatives for this integration point  [dN/dx_n]
        noalias(DN_DX) = prod(DN_De[integration_point], InvJ);

        //set shape functions for this integration point
        Vector N = row(Ncontainer, integration_point);

        //b.-compute infinitessimal strainof the composite
        this->CalculateInfinitesimalStrain(strain_vector, DN_DX);
        this->SetValue(STRAIN_VECTOR, strain_vector);

        // Compute predictive stresses for the concrete and steel
        this->CalculatePredictiveStresses(strain_vector);
    }
    KRATOS_CATCH("")
}

void RomFemDem3DElement::CalculatePredictiveStresses(const Vector &rStrainVector)
{
    auto& r_properties = this->GetProperties();
    const double Ec = r_properties[YOUNG_MODULUS];
    const double Es = r_properties[YOUNG_MODULUS_STEEL];
    const double nuc = r_properties[POISSON_RATIO];
    const double nus = r_properties[POISSON_RATIO_STEEL];

    //const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigt_size = dimension * (dimension + 1) / 2;
    Matrix constitutive_matrix_concrete = ZeroMatrix(voigt_size, voigt_size);
    Matrix constitutive_matrix_steel = ZeroMatrix(voigt_size, voigt_size);

    // Elastic rC
    this->CalculateConstitutiveMatrix(constitutive_matrix_concrete, Ec, nuc);
    this->CalculateConstitutiveMatrix(constitutive_matrix_steel, Es, nus);

    Vector r_stress_vector_concrete = prod(constitutive_matrix_concrete, rStrainVector);

    Vector stress_vector_steel;
    if (r_properties[STEEL_VOLUMETRIC_PART] > 0.0) {
        Vector elastic_strain_vector = rStrainVector - mPlasticDeformation; // E-Ep
        stress_vector_steel = prod(constitutive_matrix_steel, elastic_strain_vector);
    } else {
        stress_vector_steel = ZeroVector(voigt_size);
    }

    // Predictive Stresses
    this->SetValue(CONCRETE_STRESS_VECTOR, r_stress_vector_concrete);
    this->SetValue(STEEL_STRESS_VECTOR, stress_vector_steel);
}

void RomFemDem3DElement::CalculateAverageStressOnEdge(Vector &rAverageVector, const std::vector<Element *>& VectorOfElems)
{
    // Only averages the stress over the concrete part!!!!
    Vector current_element_stress = this->GetValue(CONCRETE_STRESS_VECTOR);
    rAverageVector = current_element_stress;
    int counter = 0;

    for (unsigned int elem = 0; elem < VectorOfElems.size(); elem++) {
        rAverageVector += VectorOfElems[elem]->GetValue(CONCRETE_STRESS_VECTOR);
        counter++;
    }
    rAverageVector /= (counter + 1);
}

void RomFemDem3DElement::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    //1.-Initialize sizes for the system components:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size = dimension * (dimension + 1) / 2;
    const unsigned int system_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != system_size)
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);

    if (rRightHandSideVector.size() != system_size)
        rRightHandSideVector.resize(system_size, false);
    noalias(rRightHandSideVector) = ZeroVector(system_size);

    Vector strain_vector(voigt_size);
    noalias(strain_vector) = ZeroVector(voigt_size);
    Vector stress_vector(voigt_size);
    noalias(stress_vector) = ZeroVector(voigt_size);
    Matrix B(voigt_size, dimension * number_of_nodes);
    noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
    Matrix DN_DX(number_of_nodes, dimension);
    noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

    //deffault values for the infinitessimal theory
    double detF = 1;
    Matrix F(dimension, dimension);
    noalias(F) = identity_matrix<double>(dimension);

    //3.-Calculate elemental system:

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    //get the shape functions [N] (for the order of the default integration method)
    const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

    //get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
    const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    //calculate delta position (here coincides with the current displacement)
    Matrix delta_position(number_of_nodes, dimension);
    noalias(delta_position) = ZeroMatrix(number_of_nodes, dimension);
    delta_position = this->CalculateDeltaPosition(delta_position);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
    GeometryType::JacobiansType J;
    J.resize(1, false);
    J[0].resize(dimension, dimension, false);
    noalias(J[0]) = ZeroMatrix(dimension, dimension);
    J = GetGeometry().Jacobian(J, mThisIntegrationMethod, delta_position);

    for (unsigned int integration_point = 0; integration_point < integration_points.size(); integration_point++) {
        const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
        Vector N = row(Ncontainer, integration_point);

        double detJ = 0;
        Matrix InvJ(dimension, dimension);
        noalias(InvJ) = ZeroMatrix(dimension, dimension);
        MathUtils<double>::InvertMatrix(J[integration_point], InvJ, detJ);

        //compute cartesian derivatives for this integration point  [dN/dx_n]
        noalias(DN_DX) = prod(DN_De[integration_point], InvJ);

        const double integration_weight = integration_points[integration_point].Weight() * detJ;
        Vector integrated_stress_concrete = ZeroVector(voigt_size);
        bool is_damaging = false;

        // Loop over edges of the element
        const Vector& characteristic_lengths = this->CalculateCharacteristicLengths();
        for (unsigned int edge = 0; edge < 6; edge++) {
            const std::vector<Element *> edge_neighbours = this->GetEdgeNeighbourElements(edge);

            Vector average_stress_concrete, average_strain_concrete, integrated_stress_vector_on_edge;
            this->CalculateAverageStressOnEdge(average_stress_concrete, edge_neighbours);
            this->CalculateAverageStrainOnEdge(average_strain_concrete, edge_neighbours);

            double damage_edge = mDamages[edge];
            double threshold = mThresholds[edge];

            // Integrate the stress on edge DAMAGE
            this->IntegrateStressDamageMechanics(threshold,
                                                 damage_edge,
                                                 average_strain_concrete, 
                                                 average_stress_concrete, 
                                                 edge, 
                                                 characteristic_lengths[edge],
                                                 is_damaging);
            mNonConvergedDamages[edge] = damage_edge;
            mNonConvergedThresholds[edge] = threshold;
        } // End loop over edges

        // Compute elemental damage
        double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);

        const Vector& r_stress_vector_concrete = this->GetValue(CONCRETE_STRESS_VECTOR);
        noalias(integrated_stress_concrete) = (1.0 - damage_element) * r_stress_vector_concrete;

        // Linear elastic const matrix concrete
        Matrix constitutive_matrix_concrete = ZeroMatrix(voigt_size, voigt_size);
        auto& r_properties = this->GetProperties();
        const double Ec = r_properties[YOUNG_MODULUS];
        const double nuc = r_properties[POISSON_RATIO];
        this->CalculateConstitutiveMatrix(constitutive_matrix_concrete, Ec, nuc);

        // Linear elastic const matrix steel
        Matrix constitutive_matrix_steel = ZeroMatrix(voigt_size, voigt_size);
        const double Es = r_properties[YOUNG_MODULUS_STEEL];
        const double nus = r_properties[POISSON_RATIO_STEEL];
        this->CalculateConstitutiveMatrix(constitutive_matrix_steel, Es, nus);

        Vector volume_force = ZeroVector(dimension);
        volume_force = this->CalculateVolumeForce(volume_force, N);

        // RHS Volumetric load
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            int index = dimension * i;
            for (unsigned int j = 0; j < dimension; j++) {
                rRightHandSideVector[index + j] += integration_weight * N[i] * volume_force[j];
            }
        }

        //compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
        Vector steel_stress_vector = this->GetValue(STEEL_STRESS_VECTOR);
        Vector integrated_steel_stress_vector = ZeroVector(voigt_size);

        // Apply plasticity steel
        this->IntegrateStressPlasticity(integrated_steel_stress_vector, steel_stress_vector, constitutive_matrix_steel);
        this->SetValue(STEEL_STRESS_VECTOR, integrated_steel_stress_vector);

        const double k = r_properties[STEEL_VOLUMETRIC_PART];
        const Vector& composite_stress_vector = k * integrated_steel_stress_vector + (1.0 - k) * integrated_stress_concrete;
        noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), composite_stress_vector);

        this->CalculateDeformationMatrix(B, DN_DX);
        const Matrix& composite_constitutive_matrix = k * constitutive_matrix_steel + (1.0 - k) * (1.0 - damage_element) * constitutive_matrix_concrete;
        noalias(rLeftHandSideMatrix) += prod(trans(B), integration_weight * Matrix(prod(composite_constitutive_matrix, B))); // LHS
    }
    KRATOS_CATCH("")
    //*****************************
}

void RomFemDem3DElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    //1.-Initialize sizes for the system components:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size = dimension * (dimension + 1) / 2;
    const unsigned int system_size = number_of_nodes * dimension;

    if (rLeftHandSideMatrix.size1() != system_size)
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);

    Matrix B(voigt_size, dimension * number_of_nodes);
    noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
    Matrix DN_DX(number_of_nodes, dimension);
    noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

    //deffault values for the infinitessimal theory
    double detF = 1;
    Matrix F(dimension, dimension);
    noalias(F) = identity_matrix<double>(dimension);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    //get the shape functions [N] (for the order of the default integration method)
    const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

    //get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
    const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    //calculate delta position (here coincides with the current displacement)
    Matrix delta_position(number_of_nodes, dimension);
    noalias(delta_position) = ZeroMatrix(number_of_nodes, dimension);
    delta_position = this->CalculateDeltaPosition(delta_position);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
    GeometryType::JacobiansType J;
    J.resize(1, false);
    J[0].resize(dimension, dimension, false);
    noalias(J[0]) = ZeroMatrix(dimension, dimension);
    J = GetGeometry().Jacobian(J, mThisIntegrationMethod, delta_position);

    for (unsigned int integration_point = 0; integration_point < integration_points.size(); integration_point++) {
        const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
        Vector N = row(Ncontainer, integration_point);

        double detJ = 0;
        Matrix InvJ(dimension, dimension);
        noalias(InvJ) = ZeroMatrix(dimension, dimension);
        MathUtils<double>::InvertMatrix(J[integration_point], InvJ, detJ);

        //compute cartesian derivatives for this integration point  [dN/dx_n]
        noalias(DN_DX) = prod(DN_De[integration_point], InvJ);

        const double integration_weight = integration_points[integration_point].Weight() * detJ;
        // Compute elemental damage
        double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);

        // Linear elastic const matrix concrete
        auto& r_properties = this->GetProperties();
        Matrix constitutive_matrix_concrete = ZeroMatrix(voigt_size, voigt_size);
        const double Ec = r_properties[YOUNG_MODULUS];
        const double nuc = r_properties[POISSON_RATIO];
        this->CalculateConstitutiveMatrix(constitutive_matrix_concrete, Ec, nuc);

        // Linear elastic const matrix steel
        Matrix constitutive_matrix_steel = ZeroMatrix(voigt_size, voigt_size);
        const double Es = r_properties[YOUNG_MODULUS_STEEL];
        const double nus = r_properties[POISSON_RATIO_STEEL];
        this->CalculateConstitutiveMatrix(constitutive_matrix_steel, Es, nus);

        const double k = r_properties[STEEL_VOLUMETRIC_PART];
        this->CalculateDeformationMatrix(B, DN_DX);
        const Matrix& composite_constitutive_matrix = k * constitutive_matrix_steel + (1.0 - k) * (1.0 - damage_element) * constitutive_matrix_concrete;
        noalias(rLeftHandSideMatrix) += prod(trans(B), integration_weight * Matrix(prod(composite_constitutive_matrix, B))); // LHS
    }
}

void RomFemDem3DElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //1.-Initialize sizes for the system components:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_size = dimension * (dimension + 1) / 2;
    const unsigned int system_size = number_of_nodes * dimension;

    if (rRightHandSideVector.size() != system_size)
        rRightHandSideVector.resize(system_size, false);
    noalias(rRightHandSideVector) = ZeroVector(system_size);


    Matrix B(voigt_size, dimension * number_of_nodes);
    noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
    Matrix DN_DX(number_of_nodes, dimension);
    noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

    // Linear elastic const matrix steel
    auto& r_properties = this->GetProperties();
    Matrix constitutive_matrix_steel = ZeroMatrix(voigt_size, voigt_size);
    const double Es = r_properties[YOUNG_MODULUS_STEEL];
    const double nus = r_properties[POISSON_RATIO_STEEL];
    this->CalculateConstitutiveMatrix(constitutive_matrix_steel, Es, nus);
    const double k = r_properties[STEEL_VOLUMETRIC_PART];

    //deffault values for the infinitessimal theory
    double detF = 1;
    Matrix F(dimension, dimension);
    noalias(F) = identity_matrix<double>(dimension);

    //3.-Calculate elemental system:

    //reading integration points
    const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    //get the shape functions [N] (for the order of the default integration method)
    const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

    //get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
    const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

    //calculate delta position (here coincides with the current displacement)
    Matrix delta_position(number_of_nodes, dimension);
    noalias(delta_position) = ZeroMatrix(number_of_nodes, dimension);
    delta_position = this->CalculateDeltaPosition(delta_position);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
    GeometryType::JacobiansType J;
    J.resize(1, false);
    J[0].resize(dimension, dimension, false);
    noalias(J[0]) = ZeroMatrix(dimension, dimension);
    J = GetGeometry().Jacobian(J, mThisIntegrationMethod, delta_position);

    for (unsigned int integration_point = 0; integration_point < integration_points.size(); integration_point++) {
        const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
        Vector N = row(Ncontainer, integration_point);

        double detJ = 0;
        Matrix InvJ(dimension, dimension);
        noalias(InvJ) = ZeroMatrix(dimension, dimension);
        MathUtils<double>::InvertMatrix(J[integration_point], InvJ, detJ);

        //compute cartesian derivatives for this integration point  [dN/dx_n]
        noalias(DN_DX) = prod(DN_De[integration_point], InvJ);

        const double integration_weight = integration_points[integration_point].Weight() * detJ;
        Vector integrated_stress_concrete = ZeroVector(voigt_size);

        // Compute elemental damage
        double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);

        const Vector &r_stress_vector_concrete = this->GetValue(CONCRETE_STRESS_VECTOR);
        noalias(integrated_stress_concrete) = (1.0 - damage_element) * r_stress_vector_concrete;

        Vector volume_force = ZeroVector(dimension);
        volume_force = this->CalculateVolumeForce(volume_force, N);

        // RHS Volumetric load
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            int index = dimension * i;
            for (unsigned int j = 0; j < dimension; j++) {
                rRightHandSideVector[index + j] += integration_weight * N[i] * volume_force[j];
            }
        }

        //compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
        Vector steel_stress_vector = this->GetValue(STEEL_STRESS_VECTOR);
        Vector integrated_steel_stress_vector = ZeroVector(voigt_size);

        // Apply plasticity steel
        this->IntegrateStressPlasticity(integrated_steel_stress_vector, steel_stress_vector, constitutive_matrix_steel);
        this->SetValue(STEEL_STRESS_VECTOR, integrated_steel_stress_vector);

        const Vector& composite_stress_vector = k * integrated_steel_stress_vector + (1.0 - k) * integrated_stress_concrete;
        noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), composite_stress_vector);
    }
}

Vector &RomFemDem3DElement::CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rVolumeForce.size() != dimension)
        rVolumeForce.resize(dimension, false);

    noalias(rVolumeForce) = ZeroVector(dimension);

    for (unsigned int j = 0; j < number_of_nodes; j++)
    {
        if (GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION))
        { // it must be checked once at the begining only
            array_1d<double, 3> &VolumeAcceleration = GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
            for (unsigned int i = 0; i < dimension; i++)
                rVolumeForce[i] += rN[j] * VolumeAcceleration[i];
        }
    }
    double k = this->GetProperties()[STEEL_VOLUMETRIC_PART];
    rVolumeForce *= (GetProperties()[DENSITY] * (1.0 - k) + GetProperties()[DENSITY_STEEL] * k);

    return rVolumeForce;

    KRATOS_CATCH("")
}

//     TENSOR VARIABLES
void RomFemDem3DElement::CalculateOnIntegrationPoints(
    const Variable<Matrix> &rVariable,
    std::vector<Matrix> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rOutput[0].size2() != dimension)
        rOutput[0].resize(dimension, dimension, false);

    if (rVariable == CONCRETE_STRESS_TENSOR) {
        rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(CONCRETE_STRESS_VECTOR));
    } else if (rVariable == STEEL_STRESS_TENSOR) {
        if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0) {
            rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(STEEL_STRESS_VECTOR));
        } else {
            Matrix dummy;
            rOutput[0] = dummy;
        }
    } else if (rVariable == STRAIN_TENSOR) {
        rOutput[0] = MathUtils<double>::StrainVectorToTensor(this->GetValue(STRAIN_VECTOR));
    } else if (rVariable == CONCRETE_STRESS_TENSOR_INTEGRATED) {
        rOutput[0] = MathUtils<double>::StressVectorToTensor((1.0 - mDamage) * this->GetValue(CONCRETE_STRESS_VECTOR));
    }
}

// Tensor variables
void RomFemDem3DElement::GetValueOnIntegrationPoints(
    const Variable<Matrix> &rVariable,
    std::vector<Matrix> &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == STRAIN_TENSOR) {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    } else if (rVariable == STEEL_STRESS_TENSOR) {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    } else if (rVariable == CONCRETE_STRESS_TENSOR) {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    } else if (rVariable == CONCRETE_STRESS_TENSOR_INTEGRATED) {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
}

// **** plasticity methods *****
void RomFemDem3DElement::IntegrateStressPlasticity(
    Vector &rIntegratedStress,
    const Vector &rPredictiveStress,
    const Matrix& rC)
{ // ecuua

    double threshold, plastic_dissipation;
    Vector plastic_strain;
    
    // Get Converged values from the prev step
    threshold = this->GetKp();
    plastic_dissipation = this->GetCapap();
    plastic_strain = this->GetPlasticDeformation();

    double uniaxial_stress, plastic_denominator;
    Vector flux_vector = ZeroVector(6), plastic_strain_increment = ZeroVector(6);

    // Compute Plastic variables
    this->CalculatePlasticParameters(rPredictiveStress, uniaxial_stress, threshold,
                                     plastic_denominator, flux_vector, plastic_dissipation,
                                     plastic_strain_increment, rC);
    double F = uniaxial_stress - threshold;

    if (F <= std::abs(1.0e-4 * threshold)) { // Elastic
        rIntegratedStress = rPredictiveStress;
        this->SetNonConvergedKp(threshold);
        this->SetNonConvergedCapap(plastic_dissipation);
        this->SetNonConvergedPlasticDeformation(plastic_strain);

        this->SetValue(EQUIVALENT_STRESS_VM, uniaxial_stress);
    } else { // Plastic case
        noalias(rIntegratedStress) = rPredictiveStress;
        double plastic_consistency_factor;
        Vector DS = ZeroVector(6), DESIG = ZeroVector(6);

        int iteration = 0;
        const int iter_max = 100;
        bool is_converged = false;

        while (is_converged == false && iteration <= iter_max) {
            plastic_consistency_factor = F * plastic_denominator;
            // if (plastic_consistency_factor < 0.0)
            //     plastic_consistency_factor = 0.0;

            noalias(plastic_strain_increment) = plastic_consistency_factor * flux_vector;
            noalias(plastic_strain) += plastic_strain_increment;
            noalias(DS) = prod(rC, plastic_strain_increment);
            noalias(DESIG) -= DS;
            noalias(rIntegratedStress) -= DS;

            this->CalculatePlasticParameters(rIntegratedStress, uniaxial_stress, threshold,
                                             plastic_denominator, flux_vector,
                                             plastic_dissipation, plastic_strain_increment, rC);

            F = uniaxial_stress - threshold;

            if (F < std::abs(1.0e-4 * threshold)) {// Has converged
                is_converged = true;
                // Update Int Vars
                this->SetNonConvergedKp(threshold);
                this->SetNonConvergedCapap(plastic_dissipation);
                this->SetNonConvergedPlasticDeformation(plastic_strain);
                this->SetValue(EQUIVALENT_STRESS_VM, uniaxial_stress);
            } else iteration++;
        }
        KRATOS_WARNING_IF("iteration", iteration == iter_max) << "Reached Max iterations inside Plasticity Loop" << std::endl;
    }
}

void RomFemDem3DElement::VonMisesYieldCriterion(
    const Vector &rStressVector,
    Vector &rDeviator,
    double &rYield,
    double &rJ2)
{
    const double I1 = this->CalculateI1Invariant(rStressVector);

    rDeviator = rStressVector;
    const double Pmean = I1 / 3.0;

    rDeviator[0] -= Pmean;
    rDeviator[1] -= Pmean;
    rDeviator[2] -= Pmean;

    rJ2 = 0.5 * (rDeviator[0] * rDeviator[0] + rDeviator[1] * rDeviator[1] + rDeviator[2] * rDeviator[2]) +
          (rDeviator[3] * rDeviator[3] + rDeviator[4] * rDeviator[4] + rDeviator[5] * rDeviator[5]);

    rYield = std::sqrt(3.0 * rJ2);
}

void RomFemDem3DElement::CalculatePlasticParameters(
    const Vector &rStressVector,
    double &rYield,
    double &rKp,
    double &rPlasticDenominator,
    Vector &rFluxVector,
    double &rPlasticDissipation,
    const Vector &rPlasticStrainIncr,
    const Matrix &rC)
{ // modaux
    Vector Deviator = ZeroVector(6), HCapa = ZeroVector(6);
    double J2 = 0.0, r0 = 0.0, r1 = 0.0, Slope = 0.0, HardeningParam = 0.0;

    this->VonMisesYieldCriterion(rStressVector, Deviator, rYield, J2);
    this->CalculateFluxVector(rStressVector, Deviator, J2, rFluxVector);
    this->CalculateRFactors(rStressVector, r0, r1);
    this->CalculatePlasticDissipation(rStressVector, r0, r1, rPlasticStrainIncr, rPlasticDissipation, HCapa);
    this->CalculateEquivalentStressThreshold(rPlasticDissipation, r0, r1, rKp, Slope);
    this->CalculateHardeningParameter(rFluxVector, Slope, HCapa, HardeningParam);
    this->CalculatePlasticDenominator(rFluxVector, rC, HardeningParam, rPlasticDenominator);
}

void RomFemDem3DElement::CalculateFluxVector(
    const Vector &rStressVector,
    const Vector &rDeviator,
    const double J2,
    Vector &rFluxVector)
{
    // Only valid for Von Mises uniaxial_stress Surf
    Vector auxiliar_vector = ZeroVector(6);
    const double denomJ2 = 1.0 / (2.0 * std::sqrt(J2));

    for (unsigned int i = 0; i < auxiliar_vector.size(); i++) {
        auxiliar_vector[i] = rDeviator[i] * denomJ2;
    }

    auxiliar_vector[3] *= 2.0;
    auxiliar_vector[4] *= 2.0;
    auxiliar_vector[5] *= 2.0;

    rFluxVector = std::sqrt(3.0) * auxiliar_vector;
}

void RomFemDem3DElement::CalculateRFactors(const Vector &rStressVector, double &r0, double &r1)
{
    Vector principal_stresses_vector = ZeroVector(3);
    this->CalculatePrincipalStresses(principal_stresses_vector, rStressVector);

    double suma = 0.0, sumb = 0.0, sumc = 0.0;
    Vector SA = ZeroVector(3), SB = ZeroVector(3), SC = ZeroVector(3);

    for (unsigned int i = 0; i < 3; i++) {
        SA[i] = std::abs(principal_stresses_vector[i]);
        SB[i] = 0.5 * (principal_stresses_vector[i] + SA[i]);
        SC[i] = 0.5 * (-principal_stresses_vector[i] + SA[i]);

        suma += SA[i];
        sumb += SB[i];
        sumc += SC[i];
    }

    if (suma != 0.0) {
        r0 = sumb / suma;
        r1 = sumc / suma;
    } else {
        r0 = sumb;
        r1 = sumc;
    }
}

void RomFemDem3DElement::CalculatePlasticDissipation(
    const Vector &rPredictiveStress,
    const double r0,
    const double r1,
    const Vector &rPlasticStrainIncrement,
    double &rPlasticDissipation,
    Vector &rHCapa)
{
    auto& r_properties = this->GetProperties();
    const double E = r_properties[YOUNG_MODULUS_STEEL];
    const double fc = r_properties[YIELD_STRESS_C_STEEL];
    const double ft = r_properties[YIELD_STRESS_T_STEEL];
    const double n = fc / ft;
    const double Gf = r_properties[FRACTURE_ENERGY_STEEL];
    const double Gfc = Gf * std::pow(n, 2);
    const double Volume = this->GetGeometry().Volume();
    const double l_char = std::pow(Volume, 1.0 / 3.0);

    const double gf = Gf / l_char;
    const double gfc = Gfc / l_char;

    const double hlim = 2.0 * E * gfc / std::pow(fc, 2);
    KRATOS_ERROR_IF(l_char > hlim) << " Characteristic length lower than minimum " << std::endl;

    double Const0 = 0.0, Const1 = 0.0;
    if (gf > 0.000001)
        Const0 = r0 / gf;
    if (gfc > 0.000001)
        Const1 = r1 / gfc;

    const double Const = Const0 + Const1;
    double Dcapa = 0.0;

    for (int i = 0; i < rPredictiveStress.size(); i++) {
        rHCapa[i] = Const * rPredictiveStress[i];
        Dcapa += rHCapa[i] * rPlasticStrainIncrement[i];
    }

    if (Dcapa < 0.0 || Dcapa > 1.0)
        Dcapa = 0.0;
    rPlasticDissipation += Dcapa;

    if (rPlasticDissipation >= 1.0)
        rPlasticDissipation = 0.9999;
}

void RomFemDem3DElement::CalculateEquivalentStressThreshold(
    const double plastic_dissipation,
    const double r0,
    const double r1,
    double &rEquivalentStressThreshold,
    double &rSlope)
{
    auto& r_properties = this->GetProperties();
    const double fc = r_properties[YIELD_STRESS_C_STEEL];
    const double ft = r_properties[YIELD_STRESS_T_STEEL];
    const double n = fc / ft;
    Vector G = ZeroVector(2), EqTrhesholds = ZeroVector(2), Slopes = ZeroVector(2);
    G[0] = r_properties[FRACTURE_ENERGY_STEEL];
    G[1] = n * n * G[0];

    const int HardCurve = r_properties[HARDENING_LAW];

    for (unsigned int i = 0; i < 2; i++) { // tension and compression curves
        switch (HardCurve)
        {
        case 1:
            this->LinearCalculateThreshold(plastic_dissipation, G[i], EqTrhesholds[i], Slopes[i]);
            break;
        case 2:
            this->ExponentialCalculateThreshold(plastic_dissipation, G[i], EqTrhesholds[i], Slopes[i]);
            break;
        case 3:
            this->HardSoftCalculateThreshold(plastic_dissipation, G[i], EqTrhesholds[i], Slopes[i]);
            break;
        default:
            KRATOS_ERROR << "Hardening law not defined..." << std::endl;
        }
    }

    rEquivalentStressThreshold = r0 * EqTrhesholds[0] + r1 * EqTrhesholds[1];
    rSlope = rEquivalentStressThreshold * ((r0 * Slopes[0] / EqTrhesholds[0]) + (r1 * Slopes[1] / EqTrhesholds[1]));
}

void RomFemDem3DElement::LinearCalculateThreshold(
    const double plastic_dissipation,
    const double Gf,
    double &rEqThreshold,
    double &rSlope)
{ // Linear softening case!!
    const double fc = this->GetProperties()[YIELD_STRESS_C_STEEL];

    rEqThreshold = fc * std::sqrt(1.0 - plastic_dissipation);
    rSlope = -0.5 * (std::pow(fc, 2.0) / (rEqThreshold));
}

void RomFemDem3DElement::ExponentialCalculateThreshold(
    const double plastic_dissipation,
    const double Gf,
    double &rEqThreshold,
    double &rSlope)
{ // Exponential softening case!!
    const double fc = this->GetProperties()[YIELD_STRESS_C_STEEL];

    rEqThreshold = fc * (1.0 - plastic_dissipation);
    rSlope = -0.5 * fc;
}

void RomFemDem3DElement::HardSoftCalculateThreshold(
    const double PlasticDissipation,
    const double Gf,
    double &rEqThreshold,
    double &rSlope)
{    // Linear Hardening followed by exp softening
    auto& r_properties = this->GetProperties();
    const double initial_threshold = r_properties[YIELD_STRESS_C_STEEL];        // sikma
    const double peak_stress = r_properties[MAXIMUM_STRESS];                    // sikpi
    const double peak_stress_position = r_properties[MAXIMUM_STRESS_POSITION]; // cappi [0,1]

    if (PlasticDissipation < 1.0) {
        const double Ro = std::sqrt(1.0 - initial_threshold / peak_stress);
        double alpha = std::log((1.0 - (1.0 - Ro) * (1.0 - Ro)) / ((3.0 - Ro) * (1.0 + Ro) * peak_stress_position));
        alpha = std::exp(alpha / (1.0 - peak_stress_position));
        const double Phi = std::pow((1.0 - Ro), 2) + ((3.0 - Ro) * (1.0 + Ro) * PlasticDissipation * (std::pow(alpha, (1.0 - PlasticDissipation))));

        rEqThreshold = peak_stress * (2.0 * std::sqrt(Phi) - Phi);
        
        rSlope = peak_stress * ((1.0 / std::sqrt(Phi)) - 1.0) * (3.0 - Ro) * (1.0 + Ro) * (std::pow(alpha, (1.0 - PlasticDissipation))) *
                 (1.0 - std::log(alpha) * PlasticDissipation);
    } else {
        KRATOS_ERROR << "The Plastic Dissipation is greater that 1.0 ..." << std::endl;
    }
}

void RomFemDem3DElement::CalculateHardeningParameter(
    const Vector &rFluxVector,
    const double SlopeThreshold,
    const Vector &rHCapa,
    double &rHardeningParam)
{
    rHardeningParam = -SlopeThreshold;
    double aux = 0.0;

    for (unsigned int i = 0; i < rFluxVector.size(); i++) {
        aux += rHCapa[i] * rFluxVector[i];
    }
    if (aux != 0.0)
        rHardeningParam *= aux;
}

void RomFemDem3DElement::CalculatePlasticDenominator(
    const Vector &rFluxVector,
    const Matrix &rElasticConstMatrix,
    const double HardeningParam,
    double &rPlasticDenominator)
{ // only for isotropic hardening
    double A1 = 0.0, A2 = 0.0, A3 = 0.0;

    const Vector Dvect = prod(rFluxVector, rElasticConstMatrix);

    for (unsigned int i = 0; i < 6; i++) {
        A1 += Dvect[i] * rFluxVector[i];
    }

    A3 = HardeningParam;
    rPlasticDenominator = 1.0 / (A1 + A2 + A3);
}

void RomFemDem3DElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    //Loop over edges
    for (unsigned int edge = 0; edge < mNumberOfEdges; edge++) {
        mDamages[edge] = mNonConvergedDamages[edge];
        mThresholds[edge] = mNonConvergedThresholds[edge];
    } // End Loop over edges

    mDamage = this->CalculateElementalDamage(mDamages);
    mThreshold = this->CalculateElementalDamage(mThresholds);

    if (mDamage >= 0.98) {
        if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0) {
            if (this->GetCapap() > 0.98) {
                this->Set(ACTIVE, false);
            }
        } else {
            this->Set(ACTIVE, false);
        }
    }

    // plasticity
    this->UpdateAndSaveInternalVariables();
    this->ResetNonConvergedVarsPlast();
}

// Double values
void RomFemDem3DElement::GetValueOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rValues,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DAMAGE_ELEMENT ||
        rVariable == IS_DAMAGED || 
        rVariable == STRESS_THRESHOLD || 
        rVariable == PLASTIC_DISSIPATION_CAPAP) {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }
}

// DOUBLE VARIABLES
void RomFemDem3DElement::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rVariable == DAMAGE_ELEMENT) {
        rOutput.resize(1);
        for (unsigned int integration_point = 0; integration_point < 1; integration_point++) {
            rOutput[integration_point] = double(this->GetValue(DAMAGE_ELEMENT));
        }
    } else if (rVariable == IS_DAMAGED) {
        rOutput.resize(1);
        for (unsigned int integration_point = 0; integration_point < 1; integration_point++) {
            rOutput[integration_point] = double(this->GetValue(IS_DAMAGED));
        }
    } else if (rVariable == STRESS_THRESHOLD) {
        rOutput.resize(1);
        for (unsigned int integration_point = 0; integration_point < 1; integration_point++) {
            rOutput[integration_point] = double(this->GetValue(STRESS_THRESHOLD));
        }
    } else if (rVariable == PLASTIC_DISSIPATION_CAPAP) {
        rOutput.resize(1);
        if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0) {
            for (unsigned int integration_point = 0; integration_point < 1; integration_point++) {
                rOutput[integration_point] = double(this->GetCapap());
            }
        } else {
            double dummy = 0.0;
            rOutput[0] = dummy;
        }
    }
}
} // namespace Kratos