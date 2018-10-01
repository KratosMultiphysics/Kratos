//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_U_PW_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PW_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_element.hpp"
#include "custom_elements/U_Pw_small_strain_interface_element.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPwSmallStrainLinkInterfaceElement : public UPwSmallStrainInterfaceElement<TDim,TNumNodes>
{

public:
    
    KRATOS_CLASS_POINTER_DEFINITION( UPwSmallStrainLinkInterfaceElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPwElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPwElement<TDim,TNumNodes>::mConstitutiveLawVector;
    typedef typename UPwSmallStrainInterfaceElement<TDim,TNumNodes>::SFGradAuxVariables SFGradAuxVariables;
    typedef typename UPwSmallStrainInterfaceElement<TDim,TNumNodes>::InterfaceElementVariables InterfaceElementVariables;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPwSmallStrainLinkInterfaceElement() : UPwSmallStrainInterfaceElement<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPwSmallStrainLinkInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPwSmallStrainInterfaceElement<TDim,TNumNodes>( NewId, pGeometry ) {}
    
    // Constructor 2
    UPwSmallStrainLinkInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPwSmallStrainInterfaceElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    // Destructor
    ~UPwSmallStrainLinkInterfaceElement() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo ) override;
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }


}; // Class UPwSmallStrainLinkInterfaceElement

} // namespace Kratos

#endif // KRATOS_U_PW_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED  defined 
