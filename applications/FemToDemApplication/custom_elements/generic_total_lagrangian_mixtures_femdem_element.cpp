//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   
// System includes

// External includes

// Project includes
#include "custom_elements/generic_total_lagrangian_mixtures_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
    if (this->mThresholds.size() != NumberOfEdges)
        this->mThresholds.resize(NumberOfEdges);
    noalias(this->mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (this->mDamages.size() != NumberOfEdges)
        this->mDamages.resize(NumberOfEdges);
    noalias(this->mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry, pProperties)
{
    // DO NOT ADD DOFS HERE!!!
    if (this->mThresholds.size() != NumberOfEdges)
        this->mThresholds.resize(NumberOfEdges);
    noalias(this->mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (this->mDamages.size() != NumberOfEdges)
        this->mDamages.resize(NumberOfEdges);
    noalias(this->mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>( NewId, this->GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::~GenericTotalLagrangianMixturesFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Pointer p_new_elem = Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/


template<unsigned int TDim, unsigned int TyieldSurf>
Vector GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>::IntegrateSmoothedConstitutiveLaw(
    const std::string& rYieldSurface,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveVariables& rThisConstVars,
    const KinematicVariables& rKinVariables,
    Vector& rStrainVector,
    double& rDamageElement,
    bool& rIsDamaging,
    const double CharacteristicLength,
    const bool SaveIntVars
    )
{
    Vector damages_edges = ZeroVector(NumberOfEdges);
    if (rYieldSurface != "Elastic") {
        // Loop over edges of the element...
        Vector average_stress_edge(VoigtSize);
        Vector average_strain_edge(VoigtSize);

        for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
            noalias(average_stress_edge) = rThisConstVars.StressVector;
            noalias(average_strain_edge) = rThisConstVars.StrainVector;
            this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
            this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

            if (!SaveIntVars) {
                damages_edges[edge] = mDamages[edge];
                double threshold = mThresholds[edge];
            
                this->IntegrateStressDamageMechanics(threshold, damages_edges[edge], average_strain_edge, 
                                                        average_stress_edge, edge, CharacteristicLength, rValues, 
                                                        rIsDamaging);
                rDamageElement = this->CalculateElementalDamage(damages_edges);         
            } else {
                this->IntegrateStressDamageMechanics(mThresholds[edge], mDamages[edge], average_strain_edge, 
                                        average_stress_edge, edge, CharacteristicLength, rValues, 
                                        rIsDamaging);
                mDamage = this->CalculateElementalDamage(mDamages);
                rDamageElement = mDamage;
            }
        } // Loop over edges
    }

    this->CalculateGreenLagrangeStrainVector(rStrainVector, rKinVariables.F);
    const Vector& r_stress_vector = rThisConstVars.StressVector;
    return (1.0 - rDamageElement)*r_stress_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericTotalLagrangianMixturesFemDemElement<2,0>;
template class GenericTotalLagrangianMixturesFemDemElement<2,1>;
template class GenericTotalLagrangianMixturesFemDemElement<2,2>;
template class GenericTotalLagrangianMixturesFemDemElement<2,3>;
template class GenericTotalLagrangianMixturesFemDemElement<2,4>;
template class GenericTotalLagrangianMixturesFemDemElement<2,5>;
template class GenericTotalLagrangianMixturesFemDemElement<2,6>;
template class GenericTotalLagrangianMixturesFemDemElement<3,0>;
template class GenericTotalLagrangianMixturesFemDemElement<3,1>;
template class GenericTotalLagrangianMixturesFemDemElement<3,2>;
template class GenericTotalLagrangianMixturesFemDemElement<3,3>;
template class GenericTotalLagrangianMixturesFemDemElement<3,4>;
template class GenericTotalLagrangianMixturesFemDemElement<3,5>;
template class GenericTotalLagrangianMixturesFemDemElement<3,6>;
} // namespace Kratos