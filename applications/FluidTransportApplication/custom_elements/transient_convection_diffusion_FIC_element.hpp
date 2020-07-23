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

#if !defined(KRATOS_TRANSIENT_CONVECTION_DIFFUSION_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_TRANSIENT_CONVECTION_DIFFUSION_FIC_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/steady_convection_diffusion_FIC_element.hpp"
#include "includes/convection_diffusion_settings.h"
#include "custom_utilities/element_utilities.hpp"
#include "utilities/element_size_calculator.h"

#include "fluid_transport_application_variables.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(FLUID_TRANSPORT_APPLICATION) TransientConvectionDiffusionFICElement : public SteadyConvectionDiffusionFICElement<TDim,TNumNodes>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( TransientConvectionDiffusionFICElement );

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
    TransientConvectionDiffusionFICElement(IndexType NewId = 0) : SteadyConvectionDiffusionFICElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    TransientConvectionDiffusionFICElement(IndexType NewId, const NodesArrayType& ThisNodes) : SteadyConvectionDiffusionFICElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    TransientConvectionDiffusionFICElement(IndexType NewId, GeometryType::Pointer pGeometry) : SteadyConvectionDiffusionFICElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    TransientConvectionDiffusionFICElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : SteadyConvectionDiffusionFICElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    virtual ~TransientConvectionDiffusionFICElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateDiffusivityVariables(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateHVector(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

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
    TransientConvectionDiffusionFICElement & operator=(TransientConvectionDiffusionFICElement const& rOther);

    /// Copy constructor.
    TransientConvectionDiffusionFICElement(TransientConvectionDiffusionFICElement const& rOther);

}; // Class TransientConvectionDiffusionFICElement

} // namespace Kratos

#endif // KRATOS_TRANSIENT_CONVECTION_DIFFUSION_FIC_ELEMENT_H_INCLUDED  defined
