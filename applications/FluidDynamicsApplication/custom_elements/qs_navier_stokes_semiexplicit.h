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

#ifndef QS_NAVIER_STOKES_SEMIEXPLICIT_H
#define QS_NAVIER_STOKES_SEMIEXPLICIT_H

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

template< unsigned int TDim, unsigned int TNumNodes >
class QSNavierStokesSemiExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of QSNavierStokesSemiExplicit
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(QSNavierStokesSemiExplicit);

    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        BoundedMatrix<double, TNumNodes, TDim > forcing;
        BoundedMatrix<double, TNumNodes, TDim > fractional_velocity;
        BoundedMatrix<double, TNumNodes, TDim > fractional_convective_velocity;
        BoundedMatrix<double, TNumNodes, TNumNodes> N_gausspoint;
        BoundedMatrix<double, TNumNodes, TDim > velocity;
        BoundedMatrix<double, TNumNodes, TDim > velocity_convective;
        BoundedMatrix<double, TNumNodes, TDim > velocity_old;

        array_1d<double, TNumNodes > N;
        array_1d<double,TNumNodes> pressure;
        array_1d<double,TNumNodes> pressure_old;

        double dt;             // Time step
        double gamma = 1;      // Fractional step splitting parameter
        double h;              // Element size
        double lumping_factor; // Factor accounting for dimension
        double mu;             // Dynamic viscosity
        double nu;             // Kinematic viscosity
        double rho;            // Density
        bool UseOSS;           // Use orthogonal subscales
        double volume;         // In 2D: element area. In 3D: element volume
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    QSNavierStokesSemiExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry);

    QSNavierStokesSemiExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties);

    /// Destructor.
    ~QSNavierStokesSemiExplicit() override = default;

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
        return Kratos::make_intrusive< QSNavierStokesSemiExplicit < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        return Kratos::make_intrusive< QSNavierStokesSemiExplicit < TDim, TNumNodes > >(NewId, pGeom, pProperties);
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
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This is called during the assembling process in order
     * to calculate the elemental right hand side vector only.
     * Note that this is explicitly forbidden as this element is
     * conceived to work with bounded arrays for the sake of efficiency.
     * Internal methods are implemented instead.
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the semi-explicit Navier-Stokes element. Call the internal method instead.";

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
        return "QSNavierStokesSemiExplicit #";
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
     * @brief Internal CalculateRightHandSide() method for momentum equation
     * This auxiliary RHS calculated method is created to bypass the element API
     * In this way bounded vectors can be used in the explicit residual calculation
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateLocalFractionalVelocitySystem(
        BoundedVector<double, (TDim) * TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Internal CalculateRightHandSide() method for mass equation
     * This auxiliary RHS calculated method is created to bypass the element API
     * @param rLeftHandSideBoundedMatrix Reference to the auxiliary LHS matrix
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateLocalPressureSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Internal CalculateRightHandSide() method for end-of-step equation
     * This auxiliary RHS calculated method is created to bypass the element API
     * @param rLeftHandSideBoundedMatrix Reference to the auxiliary LHS matrix
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateLocalEndOfStepSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
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

    void FractionalVelocityEquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const ;

    void VelocityEquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const ;

    void PressureEquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const ;

    void GetFractionalVelocityDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const ;

    void GetVelocityDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const ;

    void GetPressureDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const ;

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
    QSNavierStokesSemiExplicit() : Element()
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
    QSNavierStokesSemiExplicit& operator=(QSNavierStokesSemiExplicit const& rOther);

    /// Copy constructor
    QSNavierStokesSemiExplicit(QSNavierStokesSemiExplicit const& rOther);

    ///@}


}; // Class QSNavierStokesSemiExplicit

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
inline std::istream& operator >>(std::istream& rIStream,
                                 QSNavierStokesSemiExplicit<TDim,TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const QSNavierStokesSemiExplicit<TDim,TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // QS_NAVIER_STOKES_SEMIEXPLICIT_H
