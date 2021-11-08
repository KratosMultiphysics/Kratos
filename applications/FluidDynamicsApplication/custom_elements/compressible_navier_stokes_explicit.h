
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
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
template< unsigned int TDim, unsigned int TNumNodes>
class CompressibleNavierStokesExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Block size
    constexpr static unsigned int BlockSize = TDim + 2;

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesExplicit);

    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, BlockSize> U;
        BoundedMatrix<double, TNumNodes, BlockSize> dUdt;
        BoundedMatrix<double, TNumNodes, BlockSize> ResProj;
        BoundedMatrix<double, TNumNodes, TDim> f_ext;
        array_1d<double, TNumNodes> m_ext;
        array_1d<double, TNumNodes> r_ext;
        array_1d<double, TNumNodes> nu_sc_node;
        array_1d<double, TNumNodes> alpha_sc_node;
        array_1d<double, TNumNodes> mu_sc_nodes;
        array_1d<double, TNumNodes> beta_sc_nodes;
        array_1d<double, TNumNodes> lamb_sc_nodes;

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

        bool UseOSS;         // Use orthogonal subscales
        bool ShockCapturing; // Activate shock capturing
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
     * @brief Calculate the lumped mass vector
     * This is called during the assembling process in order
     * to calculate the elemental lumped mass vector
     * @param rLumpedMassVector the elemental lumped mass vector
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateLumpedMassVector(
        VectorType& rLumpedMassVector,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Calculate(
        const Variable<double>& rVariable,
        double& Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<Matrix>& rVariable,
        Matrix & Output,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double,3>>& rVariable,
        std::vector<array_1d<double,3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

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
        const ProcessInfo &rCurrentProcessInfo) const override;

    /**
     * Determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType &ElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const override;

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
     * @brief Internal CalculateRightHandSide() method
     * This auxiliary RHS calculated method is created to bypass the element API
     * In this way bounded vectors can be used in the explicit residual calculation
     * @param rRightHandSideBoundedVector Reference to the auxiliary RHS vector
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateRightHandSideInternal(
        BoundedVector<double, BlockSize * TNumNodes>& rRightHandSideBoundedVector,
        const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate the momentum projection
     * Auxiliary method to calculate the momentum projections for the OSS.
     * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
     * The division by the lumped mass matrix values requires to be done at the strategy level.
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate the density projection
     * Auxiliary method to calculate the denstiy projections for the OSS.
     * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
     * The division by the lumped mass matrix values requires to be done at the strategy level.
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo);

    /**
     * @brief Calculate the total energy projection
     * Auxiliary method to calculate the total energy projections for the OSS.
     * Note that this method threadsafe adds the elemental RHS values of the L2 projection to the nodes.
     * The division by the lumped mass matrix values requires to be done at the strategy level.
     * @param rCurrentProcessInfo Reference to the current process info
     */
    void CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo);

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

    /**
     * @brief Calculate the midpoint velocity divergence
     * This method calculates the velocity divergence in the midpoint of the element
     * @return double Velocity divergence in the midpoint
     */
    double CalculateMidPointVelocityDivergence() const;

    /**
     * @brief Calculate the midpoint sound velocity
     * This method calculates the speed of sound velocity in the midpoint of the element
     * @return double Speed of sound velocity in the midpoint
     */
    double CalculateMidPointSoundVelocity() const;

    /**
     * @brief Calculate the midpoint density gradient
     * This method calculates the gradient of the density in the midpoint of the element
     * @return array_1d<double,3> Density gradient in the midpoint
     */
    array_1d<double,3> CalculateMidPointDensityGradient() const;

    /**
     * @brief Calculate the midpoint temperature gradient
     * This method calculates the gradient of the temperature in the midpoint of the element
     * @return array_1d<double,3> Temperature gradient in the midpoint
     */
    array_1d<double,3> CalculateMidPointTemperatureGradient() const;

    /**
     * @brief Calculate the midpoint velocity rotational
     * This method calculates the rotational of the velocity in the midpoint of the element
     * @return array_1d<double,3> Velocity rotational in the midpoint
     */
    array_1d<double,3> CalculateMidPointVelocityRotational() const;

    /**
     * @brief Calculate the midpoint velocity gradient
     * This method calculates the gradient of the velocity in the midpoint of the element
     * @return BoundedMatrix<double, 3, 3> Velocity gradient in the midpoint
     */
    BoundedMatrix<double, 3, 3> CalculateMidPointVelocityGradient() const;

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
