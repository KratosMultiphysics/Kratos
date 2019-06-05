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
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    )
	: GenericSmallStrainFemDemElement(NewId, pGeometry, pProperties)
{
	// BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

	// Each component == Each edge
	if (mNonConvergedThresholds.size() != NumberOfEdges)
    	mNonConvergedThresholds.resize(NumberOfEdges);
	noalias(mNonConvergedThresholds) = ZeroVector(NumberOfEdges);   // Equivalent stress

	if (mThresholds.size() != NumberOfEdges)
    	mThresholds.resize(NumberOfEdges);
	noalias(mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

	if (mNonConvergedDamages.size() != NumberOfEdges)
    	mNonConvergedDamages.resize(NumberOfEdges);
	noalias(mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge

	if (mNonConvergedDamages.size() != NumberOfEdges)
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
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Create(
    IndexType NewId, 
    NodesArrayType const& rThisNodes, 
    PropertiesType::Pointer pProperties
    ) const
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

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{

}


/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateLocalSystem(
	MatrixType& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, 
    ProcessInfo& rCurrentProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericLargeDisplacementFemDemElement<2,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}


template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateB2D(
    Matrix& rB, 
    const Matrix& rF, 
    const Matrix& rDN_DX
    )
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for (SizeType i = 0; i < number_of_nodes; i++) {
        unsigned int index = 2 * i;
        rB( 0, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 0 );
        rB( 0, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 0 );
        rB( 1, index + 0 ) = rF( 0, 1 ) * rDN_DX( i, 1 );
        rB( 1, index + 1 ) = rF( 1, 1 ) * rDN_DX( i, 1 );
        rB( 2, index + 0 ) = rF( 0, 0 ) * rDN_DX( i, 1 ) + rF( 0, 1 ) * rDN_DX( i, 0 );
        rB( 2, index + 1 ) = rF( 1, 0 ) * rDN_DX( i, 1 ) + rF( 1, 1 ) * rDN_DX( i, 0 );

    }
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateB3D(
    Matrix& rB, 
    const Matrix& rF, 
    const Matrix& rDN_DX
    )
{
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
    }
}
/***********************************************************************************/
/***********************************************************************************/

































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