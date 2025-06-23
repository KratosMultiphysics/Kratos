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

#ifndef KRATOS_QS_CONVECTION_DIFFUSION_EXPLICIT_H
#define KRATOS_QS_CONVECTION_DIFFUSION_EXPLICIT_H

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
#include "convection_diffusion_application_variables.h"

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
 * @class QSConvectionDiffusionExplicit
 * @ingroup ConvectionDiffusionApplication
 * @brief This element solves the convection-diffusion equation, stabilized with
 * algebraic subgrid scale or orthogonal subgrid scale.
 * @details This element solves the convection-diffusion equation:
 * $ \frac{\partial \phi}{\partial t} + v \cdot  \nabla \phi + \phi \nabla \cdot v - \nabla \cdot alpha \nabla \phi = f $
 * where $ \phi $ is the scalar unknown, $ v $ the convective velocity, $ alpha > 0 $ the diffusivity coefficient,
 * $ f $ the forcing term.
 * Quasi-static algebraic subgrid scale and quasi-static orthogonal subgrid scale methods are exploited for stabilization.
 * The element is designed to use an explicit integration method.
 * The formulation is described in https://github.com/KratosMultiphysics/Documentation/blob/master/Resources_files/convection_diffusion_explicit_elements/Eulerian_convection_diffusion_explicit_element.pdf
 * @author Riccardo Tosi
 */
template< unsigned int TDim, unsigned int TNumNodes>
class QSConvectionDiffusionExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

        typedef Element BaseType;
        typedef Node NodeType;
        typedef Geometry<NodeType> GeometryType;

    /// Pointer definition of QSConvectionDiffusionExplicit
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(QSConvectionDiffusionExplicit);

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.
    QSConvectionDiffusionExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    QSConvectionDiffusionExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        Properties::Pointer pProperties);

    /// Default constructor.

    /// Destructor.
    virtual ~QSConvectionDiffusionExplicit();

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
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateLumpedMassVector(
        VectorType& rLumpedMassVector,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    void Calculate(
        const Variable<double>& rVariable,
        double& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "QSConvectionDiffusionExplicitElement #";
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

    struct ElementData
    {
        // scalars
        double diffusivity;
        double lumping_factor;
        double weight;
        double delta_time;
        double explicit_step_coefficient;
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
        BoundedMatrix<double,TNumNodes,TNumNodes> lhs;
        array_1d<double,TNumNodes> rhs;
        // auxiliary containers for the symbolically-generated data for Gauss integration
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
        ElementData& rData,
        const ProcessInfo& rCurrentProcessInfo);

    double ComputeH(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX);

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
    QSConvectionDiffusionExplicit() : Element()
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
    QSConvectionDiffusionExplicit& operator=(QSConvectionDiffusionExplicit const& rOther) = delete;

    /// Copy constructor
    QSConvectionDiffusionExplicit(QSConvectionDiffusionExplicit const& rOther) = delete;

    ///@}
    ///@name Private Operations
    ///@{

    void QSCalculateRightHandSideInternal(
        BoundedVector<double, TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

    void QSCalculateOrthogonalSubgridScaleRHSInternal(
        BoundedVector<double, TNumNodes>& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo);

    void QSCalculateTau(ElementData& rData);

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

}; // Class QSConvectionDiffusionExplicit

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
    QSConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1>
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const QSConvectionDiffusionExplicit<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_QS_CONVECTION_DIFFUSION_EXPLICIT_H
