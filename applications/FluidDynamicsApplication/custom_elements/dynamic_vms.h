//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_DYNAMIC_VMS_H_INCLUDED )
#define  KRATOS_DYNAMIC_VMS_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

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

/// A stabilized element for the incompressible Navier-Stokes equations.
/**
 * @todo Rewrite this documentation
 * This class implements a stabilized formulation based on the
 * Variational Multiscale framework. The the subscales can be modeled
 * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
 * In the case of OSS, the projection terms are treated explicitly (computed
 * using the results of the previous iteration) and the subscales are not
 * tracked in time. The choice of subscale model is made based on the ProcessInfo
 * variable OSS_SWITCH (OSS if 1, ASGS otherwise).
 * This class implements both the 2D and 3D versions of the element.
 *
 * The ASGS implementation follows Ramon Codina, A stabilized finite element
 * method for generalized stationary incompressible flows, Computer Methods in
 * Applied Mechanics and Engineering. Vol. 190 (2001), 2681-2706.
 *
 * The OSS implementation corresponds to the case identified as explicit, quasi-
 * static orthogonal subscales in Ramon Codina, Stabilized finite element approximation
 * of transient incompressible flows using orthogonal subscales, Computer Methods
 * in Applied Mechanics and Engineering. Vol. 191 (2002), 4295-4321.
 *
 * In addition to the stabilization, this element implements the Smagorinsky
 * model of turbulence. This turbulent term is only activated if the elemental
 * value C_SMAGORINSKY is set to something other than zero.
 *
 * This class requires at least the following variables:\n
 * On each Node, as solution step variables VELOCITY, PRESSURE, ACCELERATION, MESH_VELOCITY.\n
 * On ProcessInfo OSS_SWITCH, DELTA_TIME.\n
 * If OSS is used, the nodes also require NODAL_AREA, ADVPROJ and DIVPROJ as solution step variables.\n
 * If Smagorinsky is used, C_SMAGORINSKY has to be defined on the elements.\n
 * Error estimation stores ERROR_RATIO on the elements.\n
 * Some additional variables can be used to print results on the element: TAUONE, TAUTWO, MU, VORTICITY.
 *
 * @note Unlike VMS, this class does not use the DYNAMIC_TAU ProcessInfo value, which is always
 * assumed 1.0
 *
 * @see ResidualBasedEliminationBuilderAndSolver compatible monolithic solution strategy.
 * @see PressureSplittingBuilderAndSolver compatible segregated solution strategy.
 * @see TrilinosPressureSplittingBuilderAndSolver compatible mpi strategy.
 * @see DynamicSmagorinskyUtils to set the Smagorinsky parameter dynamically.
 * @see ResidualBasedPredictorCorrectorVelocityBossakScheme time scheme that can use
 * OSS stabilization.
 */
template< unsigned int TDim >
class DynamicVMS : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DynamicVMS
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(DynamicVMS);

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for shape function derivatives container
    typedef Kratos::Matrix ShapeDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    DynamicVMS(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constructor using a geometry object and a custom integration rule.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param ThisIntegrationMethod Integration method to use in all element integrals.
     * @note Only GI_GAUSS_1 (exact integrals for linear terms) or GI_GAUSS_2 (exact integrals up to quadratic terms) make sense as integration method.
     */
    DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, const GeometryData::IntegrationMethod ThisIntegrationMethod);

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * @brief Constructor using geometry, properties and custom integration rule.
     * @param NewId Index of the new element.
     * @param pGeometry Pointer to a geometry object.
     * @param pProperties Pointer to the element's properties.
     * @param ThisIntegrationMethod Integration method to use in all element integrals.
     * @note Only GI_GAUSS_1 (exact integrals for linear terms) or GI_GAUSS_2 (exact integrals up to quadratic terms) make sense as integration method.
     */
    DynamicVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties, const GeometryData::IntegrationMethod ThisIntegrationMethod);

    /// Destructor.
    ~DynamicVMS() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type.
    /**
     * Returns a pointer to a new DynamicVMS element, created using given input.
     * The new element will use the same number of integration points as the original.
     * @param NewId the ID of the new element.
     * @param ThisNodes the nodes of the new element.
     * @param pProperties the properties assigned to the new element.
     * @return a Pointer to the new element.
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                    PropertiesType::Pointer pProperties) const override;

    /// Create a new element of this type.
	/**
	 @param NewId Index of the new element
     @param pGeom A pointer to the geometry of the new element
	 @param pProperties Pointer to the element's properties
	 */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    /// Initialize containters for subscales on integration points.
    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Prepare the element for a new solution step.
     * Update the values on the subscales and evaluate elemental shape functions.
     * @param rCurrentProcessInfo. ProcessInfo instance (unused).
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /// Calculate a new value for the velocity subscale.
    /**
     * @param rCurrentProcessInfo ProcessInfo instance containig the time step as DELTA_TIME
     */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;


    void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                      VectorType &rRightHandSideVector,
                                     const ProcessInfo &rCurrentProcessInfo) override;


    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;


    /// Computes the local contribution associated to 'new' velocity and pressure values
    /**
     * Provides local contributions to the system associated to the velocity and
     * pressure terms (convection, diffusion, pressure gradient/velocity divergence
     * and stabilization).
     * @param rDampingMatrix Will be filled with the velocity-proportional "damping" matrix
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType &rMassMatrix, const ProcessInfo &rCurrentProcessInfo) override;


    /// Provides the global indices for each one of this element's local rows
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    /// Returns ACCELERATION_X, ACCELERATION_Y, (ACCELERATION_Z,) 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;


//    virtual void Calculate(const Variable<double>& rVariable,
//                           double& rOutput,
//                           const ProcessInfo& rCurrentProcessInfo);

    void Calculate(const Variable< array_1d<double,3> >& rVariable,
                           array_1d<double,3>& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;


    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
            std::vector<array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
            std::vector<array_1d<double, 6 > >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
            std::vector<Vector>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {}

    /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {}

    void SetValuesOnIntegrationPoints(const Variable<double> &rVariable, const std::vector<double> &rValues, const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    /// Accessor to the integration method.
    /** GiDIO uses it to determine the number of Gauss Points on each element.
     */
    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream &rOStream) const override;

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

    /// Initialize shape functions derivatives and calculate the determinant of the element's Jacobian.
    virtual void CalculateGeometryData();

    /// Calculate the value of the subscale velocity for next iteration
    /**
     * @param rCurrentProcessInfo ProcessInfo instance containing OSS_SWITCH and the time step as DELTA_TIME.
     */
    virtual void UpdateSubscale(const ProcessInfo& rCurrentProcessInfo);

    virtual void LinearUpdateSubscale(const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateASGSVelocityContribution(MatrixType& rDampingMatrix,
                                                   VectorType& rRightHandSideVector,
                                                   const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateOSSVelocityContribution(MatrixType& rDampingMatrix,
                                                  VectorType& rRightHandSideVector,
                                                  const ProcessInfo& rCurrentProcessInfo);

    virtual void LumpedMassMatrix(MatrixType& rMassMatrix);

    virtual void ConsistentMassMatrix(MatrixType& rMassMatrix);


    template< typename TValueType >
    void EvaluateInPoint(TValueType& rValue,
                         const Kratos::Variable< TValueType >& rVariable,
                         const ShapeFunctionsType& rN)
    {
        GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();

        rValue = rN[0] * rGeom[0].FastGetSolutionStepValue(rVariable);
        for (unsigned int i = 1; i < NumNodes; i++)
            rValue += rN[i] * rGeom[i].FastGetSolutionStepValue(rVariable);
    }

    template< typename TValueType >
    void EvaluateInPoint(TValueType& rValue,
                         const Kratos::Variable< TValueType >& rVariable,
                         const ShapeFunctionsType& rN,
                         const unsigned int Step)
    {
        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumNodes = rGeom.PointsNumber();

        rValue = rN[0] * rGeom[0].FastGetSolutionStepValue(rVariable,Step);
        for (unsigned int i = 1; i < NumNodes; i++)
            rValue += rN[i] * rGeom[i].FastGetSolutionStepValue(rVariable,Step);
    }

    virtual void EvaluateConvVelocity(array_1d<double,3>& rConvection,
                                      const ShapeFunctionsType& rN);

    virtual void EvaluateConvVelocity(array_1d<double,3>& rConvection,
                                      const array_1d<double,3>& rSubscaleVel,
                                      const ShapeFunctionsType& rN);

    /**
     * @brief Evaluate kinematic viscosity at the given area coordinates.
     * This function is intended to be used for the implementation of Smagorinsky in derived classes.
     * @param rViscosity Container for result
     * @param rN Array of area coordinates.
     */
    virtual void EvaluateViscosity(double& rViscosity,
                                   const ShapeFunctionsType& rN);

    virtual void ConvectionOperator(Vector& rResult,
                                    const array_1d<double,3>& rConvVel);

    /// Calculate velocity subscale stabilization parameter.
    /**
     * @param Density Density at integration point.
     * @param Viscosity Kinematic viscosity at integration point.
     * @param ConvVel Convective velocity norm (including subscale contribution).
     * @return TauOne
     */
    virtual double TauOne(const double Density,
                          const double Viscosity,
                          const double ConvVel);

    /// Calculate pressure subscale stabilization parameter
    /**
     * @param Density Density at integration point.
     * @param Viscosity Kinematic viscosity at integration point.
     * @param ConvVel Convective velocity norm (including subscale contribution).
     * @return TauTwo
     */
    virtual double TauTwo(const double Density,
                          const double Viscosity,
                          const double ConvVel);

    /// Calculate 1 / ( rho/dt + 1/TauOne )
    /**
     * @param Density Density at integration point.
     * @param Viscosity Kinematic viscosity at integration point.
     * @param ConvVel Convective velocity norm (including subscale contribution).
     * @param Dt Time Step.
     * @return TauTime
     */
    virtual double TauTime(const double Density,
                           const double Viscosity,
                           const double ConvVel,
                           const double Dt);

    /// Calculate 1 / TauTime
    virtual double InvTauTime(const double Density,
                              const double Viscosity,
                              const double ConvVel,
                              const double Dt);


    virtual void ASGSMomentumResidual(array_1d<double,3>& rResult,
                                      const double Density,
                                      const array_1d<double,3>& rConvVel,
                                      const ShapeFunctionsType& rN);

    virtual void OSSMomentumResidual(array_1d<double,3>& rResult,
                                     const double Density,
                                     const array_1d<double,3>& rConvVel,
                                     const ShapeFunctionsType& rN);

    virtual void MassResidual(double& rResult);

    virtual void AddViscousTerm(MatrixType& rDampingMatrix,
                                const double Weight,
                                const ShapeDerivativesType& rDN_DX);

    virtual void DenseSystemSolve(const Matrix& rA,
                                  const Vector& rB,
                                  Vector& rX);

    void EvaluateVorticity(array_1d<double,3>& rVorticity,
                           const ShapeDerivativesType& rDN_DX);


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

    /// Tolerance for subscale iterations (squared).
    static const double mSubscaleTol;

    /// Maximum RHS norm for subscale iterations (squared).
    static const double mSubscaleRHSTol;

    ///@}
    ///@name Member Variables
    ///@{

    GeometryData::IntegrationMethod mIntegrationMethod;

    ShapeDerivativesType mDN_DX;

    double mDetJ;

    double mElemSize;

    /// Subscale velocity evaluated on the integration point
    std::vector< array_1d<double,3> > mSubscaleVel;

    /// Subscale velocity on integration point obtained on last time step
    std::vector< array_1d<double,3> > mOldSubscaleVel;

    /// Iteration count for the non-linear velocity subscale loop
    std::vector< unsigned int > mIterCount;

//    std::vector< array_1d<double,3> > mExtraTerm;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    /// Default constructor
    DynamicVMS();

    void save(Serializer& rSerializer) const override;


    void load(Serializer& rSerializer) override;


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
    DynamicVMS & operator=(DynamicVMS const& rOther);

    /// Copy constructor.
    DynamicVMS(DynamicVMS const& rOther);

    ///@}

}; // Class DynamicVMS

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >>(std::istream& rIStream,
                                 DynamicVMS<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator <<(std::ostream& rOStream,
                                 const DynamicVMS<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} // Fluid Dynamics Application group


} // namespace Kratos.

#endif // KRATOS_DYNAMIC_VMS__H_INCLUDED  defined
