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
//                   Miguel Maso Sotomayor
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_TRANSIENT_CONVECTION_DIFFUSION_PFEM2_FIC_ELEMENT_H_INCLUDED )
#define  KRATOS_TRANSIENT_CONVECTION_DIFFUSION_PFEM2_FIC_ELEMENT_H_INCLUDED

// Project includes

// Application includes
#include "custom_elements/transient_convection_diffusion_FIC_element.hpp"


namespace Kratos
{
///@addtogroup FluidTransportApplication
///@{

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

/// Short class definition.
/** Detail class definition.
*/
template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(FLUID_TRANSPORT_APPLICATION) TransientConvectionDiffusionPFEM2FICElement : public TransientConvectionDiffusionFICElement<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TransientConvectionDiffusionPFEM2FICElement
    KRATOS_CLASS_POINTER_DEFINITION(TransientConvectionDiffusionPFEM2FICElement);

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor
    TransientConvectionDiffusionPFEM2FICElement(IndexType NewId = 0) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    TransientConvectionDiffusionPFEM2FICElement(IndexType NewId, const NodesArrayType& ThisNodes) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    TransientConvectionDiffusionPFEM2FICElement(IndexType NewId, GeometryType::Pointer pGeometry) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    TransientConvectionDiffusionPFEM2FICElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : TransientConvectionDiffusionFICElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor.
    virtual ~TransientConvectionDiffusionPFEM2FICElement() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateDiffusivityVariables(ElementVariables& rVariables, const PropertiesType& Prop, const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAndAddAdvectionMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) override {}

    void CalculateAndAddRHSAdvection(VectorType& rRightHandSideVector, ElementVariables& rVariables) override {}

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

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
    TransientConvectionDiffusionPFEM2FICElement & operator=(TransientConvectionDiffusionPFEM2FICElement const& rOther);

    /// Copy constructor.
    TransientConvectionDiffusionPFEM2FICElement(TransientConvectionDiffusionPFEM2FICElement const& rOther);

}; // Class TransientConvectionDiffusionPFEM2FICElement

} // namespace Kratos

#endif // KRATOS_TRANSIENT_CONVECTION_DIFFUSION_PFEM2_FIC_ELEMENT_H_INCLUDED  defined
