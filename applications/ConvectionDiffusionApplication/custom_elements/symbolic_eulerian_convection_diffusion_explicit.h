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

#ifndef KRATOS_SYMBOLIC_CONVECTION_DIFFUSION_EXPLICIT_H
#define KRATOS_SYMBOLIC_CONVECTION_DIFFUSION_EXPLICIT_H

// System includes


// External includes


// Project includes
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "includes/convection_diffusion_settings.h"
#include "geometries/geometry.h"

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
class SymbolicEulerianConvectionDiffusionExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SymbolicEulerianConvectionDiffusionExplicit
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicEulerianConvectionDiffusionExplicit);

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    SymbolicEulerianConvectionDiffusionExplicit(IndexType NewId, GeometryType::Pointer pGeometry);
    SymbolicEulerianConvectionDiffusionExplicit(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~SymbolicEulerianConvectionDiffusionExplicit();

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

    virtual void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

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
        return "SymbolicEulerianConvectionDiffusionExplicitElement #";
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
        double density;
        double specific_heat;
        double lumping_factor;
        double tau;
        double weight;
        // arrays
        array_1d<double,TNumNodes> forcing;
        array_1d<double,TNumNodes> unknown;
        // matrices
        BoundedMatrix<double,TNumNodes,3> convective_velocity;
        // auxiliary containers for the symbolically-generated matrices
        BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)> lhs;
        array_1d<double,TNumNodes*(TDim+1)> rhs;
        // auxiliary containers for the symbolically-generated variables for Gauss integration
        array_1d<double,TNumNodes> N;
        BoundedMatrix<double,TNumNodes,TDim> DN;
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

    void ComputeGaussPointContribution(
        ElementVariables& rVariables,
        MatrixType& rLeftHandSideMatrix,
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
    SymbolicEulerianConvectionDiffusionExplicit() : Element()
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
    SymbolicEulerianConvectionDiffusionExplicit& operator=(SymbolicEulerianConvectionDiffusionExplicit const& rOther);

    /// Copy constructor.
    SymbolicEulerianConvectionDiffusionExplicit(SymbolicEulerianConvectionDiffusionExplicit const& rOther);

    ///@}


}; // Class SymbolicEulerianConvectionDiffusionExplicit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::istream& operator >>(std::istream& rIStream,
                                 SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SymbolicEulerianConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_CONVECTION_DIFFUSION_EXPLICIT_H
