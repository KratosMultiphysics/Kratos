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


#if !defined(KRATOS_U_PL_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED )
#define  KRATOS_U_PL_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/one-phase_flow/U_Pl_element.hpp"
#include "custom_elements/one-phase_flow/U_Pl_small_strain_interface_element.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(POROMECHANICS_APPLICATION) UPlSmallStrainLinkInterfaceElement : public UPlSmallStrainInterfaceElement<TDim,TNumNodes>
{

public:
    
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UPlSmallStrainLinkInterfaceElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    using UPlElement<TDim,TNumNodes>::mThisIntegrationMethod;
    using UPlElement<TDim,TNumNodes>::mConstitutiveLawVector;
    typedef typename UPlSmallStrainInterfaceElement<TDim,TNumNodes>::SFGradAuxVariables SFGradAuxVariables;
    typedef typename UPlSmallStrainInterfaceElement<TDim,TNumNodes>::InterfaceElementVariables InterfaceElementVariables;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    UPlSmallStrainLinkInterfaceElement() : UPlSmallStrainInterfaceElement<TDim,TNumNodes>() {}
    
    // Constructor 1
    UPlSmallStrainLinkInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry) : UPlSmallStrainInterfaceElement<TDim,TNumNodes>( NewId, pGeometry ) {}
    
    // Constructor 2
    UPlSmallStrainLinkInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : UPlSmallStrainInterfaceElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    // Destructor
    ~UPlSmallStrainLinkInterfaceElement() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3>>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    
    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;

    void CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo ) override;
        
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Serialization
    
    friend class Serializer;
    
    void save(Serializer& rSerializer) const override
    {
        typedef UPlSmallStrainInterfaceElement<TDim,TNumNodes> BaseElement;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseElement )
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPlSmallStrainInterfaceElement<TDim,TNumNodes> BaseElement;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseElement )
    }


}; // Class UPlSmallStrainLinkInterfaceElement

} // namespace Kratos

#endif // KRATOS_U_PL_SMALL_STRAIN_LINK_INTERFACE_ELEMENT_H_INCLUDED  defined 
