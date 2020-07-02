//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Riccardo Tosi
//

#ifndef KRATOS_SYMBOLIC_DYNAMIC_EULERIAN_CONVECTION_DIFFUSION_EXPLICIT_H
#define KRATOS_SYMBOLIC_DYNAMIC_EULERIAN_CONVECTION_DIFFUSION_EXPLICIT_H

// System includes


// External includes


// Project includes
#include "symbolic_quasi_static_eulerian_convection_diffusion_explicit.h"

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

template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class SymbolicDynamicEulerianConvectionDiffusionExplicit : public SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

        typedef SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes> BaseType;
        typedef typename BaseType::ElementVariables ElementVariables;
        typedef Node < 3 > NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
        typedef Vector VectorType;
        typedef Matrix MatrixType;
        typedef std::size_t IndexType;
        typedef std::vector<std::size_t> EquationIdVectorType;
        typedef std::vector< Dof<double>::Pointer > DofsVectorType;
        typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Pointer definition of SymbolicDynamicEulerianConvectionDiffusionExplicit
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicDynamicEulerianConvectionDiffusionExplicit);

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    SymbolicDynamicEulerianConvectionDiffusionExplicit(IndexType NewId, GeometryType::Pointer pGeometry);
    SymbolicDynamicEulerianConvectionDiffusionExplicit(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~SymbolicDynamicEulerianConvectionDiffusionExplicit();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize(
        const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SymbolicDynamicEulerianConvectionDiffusionExplicitElement #";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << this->Id();
    }

    ///@}

protected:

    ///@name Protected member Variables
    ///@{

    VectorType mUnknownSubScale;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void ComputeGaussPointContribution(
        ElementVariables& rVariables,
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector);

    void ComputeOSSGaussPointContribution(
        ElementVariables& rVariables,
        VectorType& rRightHandSideVector);

    void UpdateUnknownSubgridScaleGaussPoint(
        ElementVariables& rVariables,
        unsigned int g);

    void CalculateTau(
        ElementVariables& rVariables);

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    // Protected default constructor necessary for serialization
    SymbolicDynamicEulerianConvectionDiffusionExplicit() : SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>()
    {
    }

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SymbolicDynamicEulerianConvectionDiffusionExplicit& operator=(SymbolicDynamicEulerianConvectionDiffusionExplicit const& rOther);

    /// Copy constructor.
    SymbolicDynamicEulerianConvectionDiffusionExplicit(SymbolicDynamicEulerianConvectionDiffusionExplicit const& rOther);

    ///@}


}; // Class SymbolicDynamicEulerianConvectionDiffusionExplicit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::istream& operator >>(std::istream& rIStream,
                                 SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SymbolicDynamicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_DYNAMIC_EULERIAN_CONVECTION_DIFFUSION_EXPLICIT_H
