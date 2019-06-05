//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "generic_large_displacement_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{


//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: GenericSmallStrainFemDemElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: GenericSmallStrainFemDemElement(NewId, pGeometry, pProperties)
{
	// BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

	// Each component == Each edge
    mNonConvergedThresholds.resize(NumberOfEdges);
	noalias(mNonConvergedThresholds) = ZeroVector(NumberOfEdges);   // Equivalent stress
    mThresholds.resize(NumberOfEdges);
	noalias(mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge
    mDamages.resize(NumberOfEdges);
	noalias(mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
    mNonConvergedDamages.resize(NumberOfEdges);
	noalias(mNonConvergedDamages) = ZeroVector(NumberOfEdges); // mDamages on edges of "i" iteration
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::GenericLargeDisplacementFemDemElement(GenericLargeDisplacementFemDemElement const& rOther)
	: GenericSmallStrainFemDemElement(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf> &GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::operator=(GenericLargeDisplacementFemDemElement const& rOther)
{
	GenericSmallStrainFemDemElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
{
	return Element::Pointer(new GenericLargeDisplacementFemDemElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
	GenericLargeDisplacementFemDemElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
	return Element::Pointer(new GenericLargeDisplacementFemDemElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::~GenericLargeDisplacementFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericLargeDisplacementFemDemElement<2,0>;
template class GenericLargeDisplacementFemDemElement<2,1>;
template class GenericLargeDisplacementFemDemElement<2,2>;
template class GenericLargeDisplacementFemDemElement<2,3>;
template class GenericLargeDisplacementFemDemElement<2,4>;
template class GenericLargeDisplacementFemDemElement<2,5>;
template class GenericLargeDisplacementFemDemElement<3,0>;
template class GenericLargeDisplacementFemDemElement<3,1>;
template class GenericLargeDisplacementFemDemElement<3,2>;
template class GenericLargeDisplacementFemDemElement<3,3>;
template class GenericLargeDisplacementFemDemElement<3,4>;
template class GenericLargeDisplacementFemDemElement<3,5>;
} // namespace Kratos