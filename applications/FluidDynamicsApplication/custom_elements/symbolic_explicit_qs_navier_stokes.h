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

#ifndef KRATOS_SYMBOLIC_EXPLICIT_QS_NAVIER_STOKES_H
#define KRATOS_SYMBOLIC_EXPLICIT_QS_NAVIER_STOKES_H

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

///@addtogroup FluidDynamicsApplication
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

template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
class SymbolicExplicitQSNavierStokes : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SymbolicExplicitQSNavierStokes
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SymbolicExplicitQSNavierStokes);

    struct ElementDataStruct
    {
        array_1d<double, TDim> forcing;
        array_1d<double, TNumNodes > N;
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        double h;           // Element size
        double volume;      // In 2D: element area. In 3D: element volume
        double mu;          // Dynamic viscosity
        double nu;          // Kinematic viscosity
        bool UseOSS;        // Use orthogonal subscales
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    SymbolicExplicitQSNavierStokes(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    SymbolicExplicitQSNavierStokes(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~SymbolicExplicitQSNavierStokes() override = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive< SymbolicExplicitQSNavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive< SymbolicExplicitQSNavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * This determines the elemental equation ID vector for all elemental DOFs.
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * Determines the elemental list of DOFs.
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * This is called during the assembling process in order to
     * calculate all elemental contributions to the global system
     * matrix and the right hand side
     * Note that this is explicitly forbidden as this element is
     * conceived to only work with explicit time integration schemes
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit Navier-Stokes element.";

        KRATOS_CATCH("");
    }

    /**
     * This is called during the assembling process in order
     * to calculate the elemental right hand side vector only.
     * Note that this is explicitly forbidden as this element is
     * conceived to work with bounded arrays for the sake of efficiency.
     * A CalculateRightHandSideInternal() method is implemented instead.
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit Navier-Stokes element. Call the CalculateRightHandSideInternal() instead.";

        KRATOS_CATCH("");
    }

    /**
     * This is called during the assembling process in order
     * to calculate the elemental contribution in explicit calculation.
     * NodalData is modified inside the function, so the
     * "AddEXplicit" FUNCTIONS ARE THE ONLY FUNCTIONS IN WHICH AN ELEMENT
     * IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
      * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * This is called during the assembling process in order
     * to calculate the elemental mass matrix.
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    // void Calculate(
    //     const Variable<double>& rVariable,
    //     double& Output,
    //     const ProcessInfo& rCurrentProcessInfo) override;

    // void CalculateOrthogonalSubgridScaleSystem(
    //     VectorType& rRightHandSideVector,
    //     const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string
    std::string Info() const override
    {
        return "SymbolicExplicitQSNavierStokes #";
    }

    /// Print information about this object
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}

protected:

    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Internal CalculateRightHandSide() method
     * This auxiliary RHS calculated method is created to bypass the element API
     * In this way bounded vectors can be used in the explicit residual calculation
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateRightHandSideInternal(
        BoundedVector<double, (TDim+1) * TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

    // void CalculateOrthogonalSubgridScaleSystemInternal(
    //     ElementDataStruct& rData,
    //     VectorType& rRightHandSideVector);

    /**
     * @brief Fill element data
     * Auxiliary function to fill the element data structure
     * @param rData Reference to the element data structure to be filled
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void FillElementData(
        ElementDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate the element size
     * This function calculates and returns the element size from the shape function gradients
     * @param rDN_DX Reference to the shape functions container
     * @return double The computed element size
     */
    double CalculateElementSize(
        const BoundedMatrix<double,TNumNodes, TDim>& rDN_DX);


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
    SymbolicExplicitQSNavierStokes() : Element()
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

    /// Assignment operator
    SymbolicExplicitQSNavierStokes& operator=(SymbolicExplicitQSNavierStokes const& rOther);

    /// Copy constructor
    SymbolicExplicitQSNavierStokes(SymbolicExplicitQSNavierStokes const& rOther);

    ///@}


}; // Class SymbolicExplicitQSNavierStokes

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
inline std::istream& operator >>(std::istream& rIStream,
                                 SymbolicExplicitQSNavierStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const SymbolicExplicitQSNavierStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_SYMBOLIC_EXPLICIT_QS_NAVIER_STOKES_H
