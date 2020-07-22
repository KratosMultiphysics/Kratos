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

#ifndef KRATOS_SYMBOLIC_QUASI_STATIC_EULERIAN_CONVECTION_DIFFUSION_EXPLICIT_H
#define KRATOS_SYMBOLIC_QUASI_STATIC_EULERIAN_CONVECTION_DIFFUSION_EXPLICIT_H

// System includes


// External includes


// Project includes
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

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
 * @class SymbolicQuasiStaticEulerianConvectionDiffusionExplicit
 * @ingroup ConvectionDiffusionApplication
 * @brief This element solves the convection-diffusion equation, stabilized with
 * algebraic subgrid scale or orthogonal subgrid scale.
 * @details This element solves the convection-diffusion equation:
 * $ \frac{\partial \phi}{\partial t} + v \cdot  \nabla \phi + \phi \nabla \cdot v - \nabla \cdot k \nabla \phi = f $
 * where $ \phi $ is the scalar unknown, $ v $ the convective velocity, $ k > 0 $ the diffusivity coefficient,
 * $ f $ the forcing term.
 * Quasi-static algebraic subgrid scale and quasi-static orthogonal subgrid scale methods are exploited for stabilization.
 * The element is designed to use an explicit integration method.
 * @author Riccardo Tosi
 */
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class SymbolicQuasiStaticEulerianConvectionDiffusionExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

        typedef Element BaseType;
        typedef Node < 3 > NodeType;
        typedef Geometry<NodeType> GeometryType;

    /// Pointer definition of SymbolicQuasiStaticEulerianConvectionDiffusionExplicit
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicQuasiStaticEulerianConvectionDiffusionExplicit);

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    SymbolicQuasiStaticEulerianConvectionDiffusionExplicit(IndexType NewId, GeometryType::Pointer pGeometry);
    SymbolicQuasiStaticEulerianConvectionDiffusionExplicit(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~SymbolicQuasiStaticEulerianConvectionDiffusionExplicit();

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

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
        DofsVectorType& ElementalDofList,
        ProcessInfo& CurrentProcessInfo) override;

    void AddExplicitContribution(
        ProcessInfo &rCurrentProcessInfo) override;

    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    void Calculate(
        const Variable<double>& rVariable,
        double& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOrthogonalSubgridScaleSystem(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo);

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
        return "SymbolicQuasiStaticEulerianConvectionDiffusionExplicitElement #";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    ///@}

protected:

    ///@name Protected member Variables
    ///@{

    struct ElementVariables
    {
        // scalars
        double diffusivity;
        double lumping_factor;
        double weight;
        double delta_time;
        double RK_time_coefficient;
        double dynamic_tau;
        double unknown_subscale;
	    double volume;
        // arrays
	    array_1d<double,TNumNodes> tau;
        array_1d<double,TNumNodes> forcing;
        array_1d<double,TNumNodes> unknown;
        array_1d<double,TNumNodes> unknown_old;
        array_1d<double,TNumNodes> oss_projection;
        // matrices
        BoundedMatrix<double,TNumNodes,3> convective_velocity;
        // auxiliary containers for the symbolically-generated matrices
        BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs;
        array_1d<double,TNumNodes*(TDim+1)> rhs;
        // auxiliary containers for the symbolically-generated variables for Gauss integration
        array_1d<double,TNumNodes> N;
        BoundedMatrix<double,TNumNodes,TNumNodes> N_gausspoint;
	    BoundedMatrix<double,TNumNodes,TDim> DN_DX;
    };

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void InitializeEulerianElement(
        ElementVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystemInternal(
        ElementVariables& rVariables,
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector);

    void CalculateOrthogonalSubgridScaleSystemInternal(
        ElementVariables& rVariables,
        VectorType& rRightHandSideVector);

    double ComputeH(
        BoundedMatrix<double,TNumNodes,TDim>& rDN_DX);

    void CalculateTau(
        ElementVariables& rVariables);

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{

    IntegrationMethod GetIntegrationMethod() const override;

    ///@}
    ///@name Protected LifeCycle
    ///@{

    // Protected default constructor necessary for serialization
    SymbolicQuasiStaticEulerianConvectionDiffusionExplicit() : Element()
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
    SymbolicQuasiStaticEulerianConvectionDiffusionExplicit& operator=(SymbolicQuasiStaticEulerianConvectionDiffusionExplicit const& rOther);

    /// Copy constructor.
    SymbolicQuasiStaticEulerianConvectionDiffusionExplicit(SymbolicQuasiStaticEulerianConvectionDiffusionExplicit const& rOther);

    ///@}


}; // Class SymbolicQuasiStaticEulerianConvectionDiffusionExplicit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::istream& operator >>(std::istream& rIStream,
                                 SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SymbolicQuasiStaticEulerianConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_QUASI_STATIC_EULERIAN_CONVECTION_DIFFUSION_EXPLICIT_H
