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

#ifndef KRATOS_SYMBOLIC_D_CONVECTION_DIFFUSION_EXPLICIT_H
#define KRATOS_SYMBOLIC_D_CONVECTION_DIFFUSION_EXPLICIT_H

// System includes


// External includes


// Project includes
#include "symbolic_qs_convection_diffusion_explicit.h"

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
 * @class SymbolicDConvectionDiffusionExplicit
 * @ingroup ConvectionDiffusionApplication
 * @brief This element solves the convection-diffusion equation, stabilized with
 * algebraic subgrid scale or orthogonal subgrid scale.
 * @details This element solves the convection-diffusion equation:
 * $ \frac{\partial \phi}{\partial t} + v \cdot  \nabla \phi + \phi \nabla \cdot v - \nabla \cdot k \nabla \phi = f $
 * where $ \phi $ is the scalar unknown, $ v $ the convective velocity, $ k > 0 $ the diffusivity coefficient,
 * $ f $ the forcing term.
 * Dynamic algebraic subgrid scale and dynamic orthogonal subgrid scale methods are exploited for stabilization.
 * The element is designed to use an explicit integration method.
 * @author Riccardo Tosi
 */
template< unsigned int TDim, unsigned int TNumNodes>
class SymbolicDConvectionDiffusionExplicit : public SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

        typedef SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes> BaseType;
        typedef typename BaseType::ElementData ElementData;
        typedef Node < 3 > NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
        typedef Vector VectorType;
        typedef Matrix MatrixType;
        typedef std::size_t IndexType;
        typedef std::vector<std::size_t> EquationIdVectorType;
        typedef std::vector< Dof<double>::Pointer > DofsVectorType;
        typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Pointer definition of SymbolicDConvectionDiffusionExplicit
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicDConvectionDiffusionExplicit);

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.
    SymbolicDConvectionDiffusionExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry);
    SymbolicDConvectionDiffusionExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Default constuctor.

    /// Destructor.
    virtual ~SymbolicDConvectionDiffusionExplicit();

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
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit Convection-Diffusion element.";
        KRATOS_CATCH("");
    }

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit Convection-Diffusion element. Call the CalculateRightHandSideInternal() instead.";
        KRATOS_CATCH("");
    }

    void AddExplicitContribution(
        const ProcessInfo &rCurrentProcessInfo) override;

    void Initialize(
        const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<double>& rVariable,
        double& Output,
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
        return "SymbolicDConvectionDiffusionExplicitElement #";
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

    BoundedVector<double, TNumNodes> mUnknownSubScale;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

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
    SymbolicDConvectionDiffusionExplicit() : SymbolicQSConvectionDiffusionExplicit<TDim,TNumNodes>()
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

    /// Assignment operator
    SymbolicDConvectionDiffusionExplicit& operator=(SymbolicDConvectionDiffusionExplicit const& rOther) = delete;
    /// Copy constructor
    SymbolicDConvectionDiffusionExplicit(SymbolicDConvectionDiffusionExplicit const& rOther) = delete;

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateRightHandSideInternal(
        BoundedVector<double, TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateOrthogonalSubgridScaleRHSInternal(
        BoundedVector<double, TNumNodes>& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo);

    void UpdateUnknownSubgridScaleGaussPoint(
        ElementData& rData,
        unsigned int g);

    void CalculateTau(
        ElementData& rData);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SymbolicDConvectionDiffusionExplicit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::istream& operator >>(
    std::istream& rIStream,
    SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const SymbolicDConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_D_CONVECTION_DIFFUSION_EXPLICIT_H
