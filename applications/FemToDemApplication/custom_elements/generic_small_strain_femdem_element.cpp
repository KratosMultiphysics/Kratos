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

// System includes

// External includes

// Project includes
#include "custom_elements/generic_small_strain_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry)
{
    BaseType::(NewId, pGeometry);
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry, pProperties)
{
    BaseType::(NewId, pGeometry, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericSmallStrainFemDemElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericSmallStrainFemDemElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim,TyieldSurf>::~GenericSmallStrainFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    BaseType::Clone(NewId, rThisNodes);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType variables;
    this->InitializeElementData(variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions = values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    for (SizeType point_number = 0; point_number < r_integration_points.size(); point_number++) {
        //compute element kinematic variables B, F, DN_DX ...
        this->CalculateKinematics(variables,point_number);

        //calculate material response
        this->CalculateMaterialResponse(variables, values, point_number);
    }
    this->SetValue(STRESS_VECTOR, values.GetStressVector());
    this->SetValue(STRAIN_VECTOR, values.GetStrainVector());

    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateTangentTensor(
    Matrix& rTangentTensor,
    const Vector& rStrainVectorGP,
    const Vector& rStressVectorGP,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double number_components = rStrainVectorGP.size();
    rTangentTensor.resize(number_components, number_components);
    Vector perturbed_stress, perturbed_strain;
    perturbed_strain.resize(number_components);
    perturbed_stress.resize(number_components);
    
    for (unsigned int component = 0; component < number_components; component++) {
        double perturbation;
        this->CalculatePerturbation(rStrainVectorGP, perturbation, component);
        this->PerturbateStrainVector(perturbed_strain, rStrainVectorGP, perturbation, component);
        this->IntegratePerturbedStrain(perturbed_stress, perturbed_strain, rElasticMatrix, rValues);
        const Vector& r_delta_stress = perturbed_stress - rStressVectorGP;
        this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress, perturbation, component);
    }
}


/***********************************************************************************/
/***********************************************************************************/

// template class GenericSmallStrainFemDemElement<2,0>;
template class GenericSmallStrainFemDemElement<2,1>;
// template class GenericSmallStrainFemDemElement<2,2>;
// template class GenericSmallStrainFemDemElement<2,3>;
// template class GenericSmallStrainFemDemElement<2,4>;
// template class GenericSmallStrainFemDemElement<2,5>;
// template class GenericSmallStrainFemDemElement<2,6>;
// template class GenericSmallStrainFemDemElement<3,0>;
// template class GenericSmallStrainFemDemElement<3,1>;
// template class GenericSmallStrainFemDemElement<3,2>;
// template class GenericSmallStrainFemDemElement<3,3>;
// template class GenericSmallStrainFemDemElement<3,4>;
// template class GenericSmallStrainFemDemElement<3,5>;
// template class GenericSmallStrainFemDemElement<3,6>;
} // namespace Kratos