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

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_FEMDEM_ELEMENT_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_FEMDEM_ELEMENT_H_INCLUDED


// System includes


// External include

// Project includes

#include "custom_elements/solid_elements/small_displacement_element.hpp"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class GenericSmallStrainFemDemElement
 * @ingroup FemToDemApplication
 * @brief Small Displacement element for the 2D and 3D cases
 * @author Alejandro Cornejo
 */
template<unsigned int TDim, unsigned int TyieldSurf>
class GenericSmallStrainFemDemElement 
    : public SmallDisplacementElement // Derived Element from SolidMechanics
{
public:
    ///@name Type Definitions
    ///@{

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef SmallDisplacementElement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of GenericSmallStrainFemDemElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GenericSmallStrainFemDemElement);

    ///@}
    ///@name Life Cycle
    ///@{

	/// Default constructors
	GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry);
	GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

	///Copy constructor
	GenericSmallStrainFemDemElement(GenericSmallStrainFemDemElement const &rOther);

	/// Destructor.
	virtual ~GenericSmallStrainFemDemElement();

	/// Assignment operator.
	GenericSmallStrainFemDemElement &operator=(GenericSmallStrainFemDemElement const &rOther);
	Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;
	Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

	GenericSmallStrainFemDemElement()
	{
	}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

	int mNumberOfEdges;
	Vector mNonConvergedThresholds;     // Equivalent stress
	Vector mThresholds;                 // Stress mThreshold on edge
	Vector mDamages;                    // Converged Damage on each edge
	Vector mNonConvergedDamages;        // Damages at edges of "i" iteration
	double mThreshold = 0.0;            // Converged Threshold
	double mDamage = 0.0;               // Converged Damage

    // Vector to storage the neigh elements sharing a certain edge
    std::vector<std::vector<Element*>> mEdgeNeighboursContainer;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
}; // Class TotalLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_H_INCLUDED  defined