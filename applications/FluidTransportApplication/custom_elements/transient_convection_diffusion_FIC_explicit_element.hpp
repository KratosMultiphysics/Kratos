//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_TRANSIENT_CONVECTION_DIFFUSION_FIC_EXPLICIT_ELEMENT_H_INCLUDED )
#define  KRATOS_TRANSIENT_CONVECTION_DIFFUSION_FIC_EXPLICIT_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/transient_convection_diffusion_FIC_element.hpp"
#include "includes/convection_diffusion_settings.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/element_size_calculator.h"

#include "fluid_transport_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(FLUID_TRANSPORT_APPLICATION) TransientConvectionDiffusionFICExplicitElement : public TransientConvectionDiffusionFICElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( TransientConvectionDiffusionFICExplicitElement );

    typedef std::size_t IndexType;
	typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef typename Element::DofsVectorType DofsVectorType;
    typedef typename Element::EquationIdVectorType EquationIdVectorType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;
    typedef typename SteadyConvectionDiffusionFICElement<TDim,TNumNodes>::ElementVariables ElementVariables;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    TransientConvectionDiffusionFICExplicitElement(IndexType NewId = 0) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    TransientConvectionDiffusionFICExplicitElement(IndexType NewId, const NodesArrayType& ThisNodes) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    TransientConvectionDiffusionFICExplicitElement(IndexType NewId, GeometryType::Pointer pGeometry) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    TransientConvectionDiffusionFICExplicitElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    virtual ~TransientConvectionDiffusionFICExplicitElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeElementVariables(ElementVariables& rVariables, const GeometryType& Geom, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    /// Assignment operator.
    TransientConvectionDiffusionFICExplicitElement & operator=(TransientConvectionDiffusionFICExplicitElement const& rOther);

    /// Copy constructor.
    TransientConvectionDiffusionFICExplicitElement(TransientConvectionDiffusionFICExplicitElement const& rOther);

}; // Class TransientConvectionDiffusionFICExplicitElement

} // namespace Kratos

#endif // KRATOS_TRANSIENT_CONVECTION_DIFFUSION_FIC_EXPLICIT_ELEMENT_H_INCLUDED  defined
