
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla (based on Elisa Magliozzi previous work)
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED)
#define  KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED

// System includes

// External includes

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
 * @brief Compressible Navier-Stokes explicit element
 * This element implements a compressible Navier-Stokes explicit formulation.
 * The formulation is written in conservative form so the element unknowns are
 * the DENSITY, MOMENTUM and TOTAL_ENERGY variables.
 * This element is intended to work with the Kratos explicit DOF based strategy.
 * Hence, the explicit residual is written in the corresponding REACTION variables.
 * @tparam TDim The space dimension (2 or 3)
 * @tparam TNumNodes The number of nodes
 */
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1, unsigned int TBlockSize = TDim + 2 >
class CompressibleNavierStokesExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesExplicit);

    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, TBlockSize> U;
        BoundedMatrix<double, TNumNodes, TDim> f_ext;
        array_1d<double, TNumNodes> r; // At the moment considering all parameters as constant in the domain (mu, nu, etc...)
        array_1d<double, TDim> f_gauss;
        double r_gauss;

        array_1d<double, TNumNodes > N;
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;

        double h;           // Element size
        double volume;      // In 2D: element area. In 3D: element volume
        double mu;          // Dynamic viscosity
        double nu;          // Kinematic viscosity
        double nu_sc;       // Kinematic viscosity (shock capturing)
        double lambda;      // Heat conductivity
        double lambda_sc;   // Heat conductivity (shock capturing)
        double c_v;         // Heat capacity at constant volume
        double gamma;       // Heat capacity ratio
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CompressibleNavierStokesExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    CompressibleNavierStokesExplicit(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~CompressibleNavierStokesExplicit() override = default;

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
        KRATOS_TRY
        return Kratos::make_intrusive< CompressibleNavierStokesExplicit < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< CompressibleNavierStokesExplicit < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

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
        KRATOS_TRY

        KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit compressible Navier-Stokes element.";

        KRATOS_CATCH("")
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
        KRATOS_TRY

        KRATOS_ERROR << "Calling the CalculateRightHandSide() method for the explicit compressible Navier-Stokes element. Call the CalculateRightHandSideInternal() instead.";

        KRATOS_CATCH("")
    }

    /**
     * This is called during the assembling process in order
     * to calculate the elemental contribution in explicit calculation.
     * NodalData is modified Inside the function, so the
     * The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT
     * IS ALLOWED TO WRITE ON ITS NODES.
     * the caller is expected to ensure thread safety hence
     * SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
      * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * This is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateMassMatrix(
        MatrixType &rMassMatrix,
        const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CompressibleNavierStokesExplicit #";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType &rResult,
        ProcessInfo &rCurrentProcessInfo) override;

    /**
     * Determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType &ElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const override;

    /**
     * @brief Calculates the shock capturing values
     * This function is intended to calculate the shock capturing values
     * These are the shock capturing viscosity and thermal diffusivity
     * @param rData Reference to the element data container
     */
    void CalculateShockCapturingValues(ElementDataStruct &rData) const;

    ///@}
    ///@name Protected Operators
    ///@{

    CompressibleNavierStokesExplicit() = default;

    ///@}
    ///@name Protected Operations
    ///@{

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
    double CalculateElementSize(const BoundedMatrix<double,TNumNodes, TDim>& rDN_DX);

    /**
     * @brief Internal CalculateRightHandSide() method
     * This auxiliary RHS calculated method is created to bypass the element API
     * In this way bounded vectors can be used in the explicit residual calculation
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Refeecen to the current process inf
     */
    void CalculateRightHandSideInternal(
        BoundedVector<double, TBlockSize * TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

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

    ///@}
};
///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED  defined
