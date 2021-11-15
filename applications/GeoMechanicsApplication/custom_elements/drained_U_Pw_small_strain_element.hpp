// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_DRAINED_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_DRAINED_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/comparison_utilities.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) DrainedUPwSmallStrainElement
 : public UPwSmallStrainElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( DrainedUPwSmallStrainElement );

    typedef std::size_t IndexType;
    typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    /// The definition of the sizetype
    typedef std::size_t SizeType;
    using UPwSmallStrainElement<TDim,TNumNodes>::mConstitutiveLawVector;
    using UPwSmallStrainElement<TDim,TNumNodes>::mStressVector;
    using UPwSmallStrainElement<TDim,TNumNodes>::mStateVariablesFinalized;
    typedef typename UPwSmallStrainElement<TDim,TNumNodes>::ElementVariables ElementVariables;
    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    DrainedUPwSmallStrainElement(IndexType NewId = 0) : UPwSmallStrainElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    DrainedUPwSmallStrainElement(IndexType NewId,
                                 const NodesArrayType& ThisNodes) :
                                 UPwSmallStrainElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    DrainedUPwSmallStrainElement(IndexType NewId,
                                 GeometryType::Pointer pGeometry) :
                                 UPwSmallStrainElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    DrainedUPwSmallStrainElement(IndexType NewId,
                                 GeometryType::Pointer pGeometry,
                                 PropertiesType::Pointer pProperties) :
                                 UPwSmallStrainElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~DrainedUPwSmallStrainElement() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Element::Pointer Create( IndexType NewId,
                             NodesArrayType const& ThisNodes,
                             PropertiesType::Pointer pProperties ) const override;

    Element::Pointer Create( IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties ) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) override;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint) override;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    /// Assignment operator.
    DrainedUPwSmallStrainElement & operator=(DrainedUPwSmallStrainElement const& rOther);

    /// Copy constructor.
    DrainedUPwSmallStrainElement(DrainedUPwSmallStrainElement const& rOther);

}; // Class DrainedUPwSmallStrainElement

} // namespace Kratos

#endif // KRATOS_GEO_DRAINED_U_PW_SMALL_STRAIN_ELEMENT_H_INCLUDED  defined
