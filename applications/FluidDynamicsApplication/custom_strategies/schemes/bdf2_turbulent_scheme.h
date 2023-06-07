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

#if !defined(KRATOS_BDF2_TURBULENT_SCHEME_H_INCLUDED )
#define  KRATOS_BDF2_TURBULENT_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "includes/define.h"
// #include "includes/serializer.h"
#include "includes/dof.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"
#include "utilities/coordinate_transformation_utilities.h"

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

/// A scheme for BDF2 time integration.
/**
 */
template<class TSparseSpace,class TDenseSpace>
class BDF2TurbulentScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BDF2TurbulentScheme
    KRATOS_CLASS_POINTER_DEFINITION(BDF2TurbulentScheme);
    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;
    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Dof<TDataType> TDofType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef CoordinateTransformationUtils<LocalSystemMatrixType, LocalSystemVectorType, double> RotationToolType;
    typedef typename RotationToolType::UniquePointer RotationToolPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BDF2TurbulentScheme()
    : Scheme<TSparseSpace, TDenseSpace>()
    , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}

    /// Constructor to use the formulation combined with a turbulence model.
    /**
     * The turbulence model is assumed to be implemented as a Kratos::Process.
     * The model's Execute() method wil be called at the start of each
     * non-linear iteration.
     * @param pTurbulenceModel pointer to the turbulence model
     */
    BDF2TurbulentScheme(Process::Pointer pTurbulenceModel)
        : Scheme<TSparseSpace, TDenseSpace>()
        , mpTurbulenceModel(pTurbulenceModel)
        , mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}

    /// Constructor for periodic boundary conditions.
    /**
     * @param rPeriodicVar the variable used to store periodic pair indices.
     */
    BDF2TurbulentScheme(const Kratos::Variable<int>& rPeriodicVar)
        : Scheme<TSparseSpace, TDenseSpace>()
        , mrPeriodicIdVar(rPeriodicVar)
    {}


    /// Destructor.
    ~BDF2TurbulentScheme() override
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Check input data for errors.
    /**
     * @param rModelPart The fluid's ModelPart
     * @return 0 if no errors were found
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        // Base scheme check
        int error_code = BaseType::Check(rModelPart);
        if (error_code != 0) {
            return error_code;
        }

        // Check buffer size
        KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 3)
            << "Insufficient buffer size for BDF2, should be at least 3, got " << rModelPart.GetBufferSize() << std::endl;

        return 0;

        KRATOS_CATCH("");
    }

    void Initialize(ModelPart& rModelPart) override
    {
        // Set up the rotation tool pointer
        const auto& r_proces_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size = r_proces_info[DOMAIN_SIZE];
        auto p_aux = Kratos::make_unique<RotationToolType>(domain_size, domain_size + 1, SLIP);
        mpRotationTool.swap(p_aux);

        // Base initialize call
        BaseType::Initialize(rModelPart);
    }

    /// Set the time iteration coefficients
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        this->SetTimeCoefficients(rModelPart.GetProcessInfo());

        // Base function initializes elements and conditions
        BaseType::InitializeSolutionStep(rModelPart,A,Dx,b);

        // Recalculate mesh velocity (to account for variable time step)
        const double tol = 1.0e-12;
        const double Dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        const double OldDt = rModelPart.GetProcessInfo().GetPreviousSolutionStepInfo(1)[DELTA_TIME];
        if(std::abs(Dt - OldDt) > tol) {
            const int n_nodes = rModelPart.NumberOfNodes();
            const Vector& BDFcoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];

#pragma omp parallel for
            for(int i_node = 0; i_node < n_nodes; ++i_node) {
                auto it_node = rModelPart.NodesBegin() + i_node;
                auto& rMeshVel = it_node->FastGetSolutionStepValue(MESH_VELOCITY);
                const auto& rDisp0 = it_node->FastGetSolutionStepValue(DISPLACEMENT);
                const auto& rDisp1 = it_node->FastGetSolutionStepValue(DISPLACEMENT,1);
                const auto& rDisp2 = it_node->FastGetSolutionStepValue(DISPLACEMENT,2);
                rMeshVel = BDFcoefs[0] * rDisp0 + BDFcoefs[1] * rDisp1 + BDFcoefs[2] * rDisp2;
            }
        }
    }

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (mpTurbulenceModel != 0) mpTurbulenceModel->Execute();

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //if orthogonal subscales are computed
        if (CurrentProcessInfo[OSS_SWITCH] == 1.0)
        {
            this->LumpedProjection(rModelPart);
            //this->FullProjection(rModelPart);
        }

    }

    /// Start the iteration by providing a first approximation to the solution.
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        const int n_nodes = rModelPart.NumberOfNodes();
        const Vector& BDFcoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];

#pragma omp parallel for
        for(int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;
            auto& rVel0 = it_node->FastGetSolutionStepValue(VELOCITY);
            const auto& rVel1 = it_node->FastGetSolutionStepValue(VELOCITY,1);
            const auto& rVel2 = it_node->FastGetSolutionStepValue(VELOCITY,2);
            auto& rAcceleration = it_node->FastGetSolutionStepValue(ACCELERATION);

            // Predict velocities
            if(!it_node->IsFixed(VELOCITY_X))
                rVel0[0] = 2.00 * rVel1[0] - rVel2[0];
            if(!it_node->IsFixed(VELOCITY_Y))
                rVel0[1] = 2.00 * rVel1[1] - rVel2[1];
            if(!it_node->IsFixed(VELOCITY_Z))
                rVel0[2] = 2.00 * rVel1[2] - rVel2[2];

            // Predict acceleration
            rAcceleration = BDFcoefs[0] * rVel0 + BDFcoefs[1] * rVel1 + BDFcoefs[2] * rVel2;
        }

        KRATOS_CATCH("")
    }

    /// Store the iteration results as solution step variables and update acceleration after a Newton-Raphson iteration.
    /**
     * @param rModelPart fluid ModelPart
     * @param rDofSet DofSet containing the Newton-Raphson system degrees of freedom.
     * @param A Newton-Raphson system matrix (unused)
     * @param Dx Newton-Raphson iteration solution
     * @param b Newton-Raphson right hand side (unused)
     */
    virtual void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        mpRotationTool->RotateVelocities(rModelPart);

        mpDofUpdater->UpdateDofs(rDofSet,Dx);

        mpRotationTool->RecoverVelocities(rModelPart);

        const Vector& BDFCoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];
        this->UpdateAcceleration(rModelPart,BDFCoefs);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Get Equation Id
        rCurrentElement.EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rCurrentElement.CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentElement.CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->CombineLHSContributions(LHS_Contribution,Mass,Damp,rCurrentProcessInfo);
        this->AddDynamicRHSContribution<Kratos::Element>(rCurrentElement,RHS_Contribution,Mass,rCurrentProcessInfo);

        // Apply slip condition
        mpRotationTool->Rotate(LHS_Contribution, RHS_Contribution, rCurrentElement.GetGeometry());
        mpRotationTool->ApplySlipCondition(LHS_Contribution, RHS_Contribution, rCurrentElement.GetGeometry());

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType &RHS_Contribution,
        Element::EquationIdVectorType &rEquationId,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Get Equation Id
        rCurrentElement.EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentElement.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        rCurrentElement.CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentElement.CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->AddDynamicRHSContribution<Kratos::Element>(rCurrentElement,RHS_Contribution,Mass,rCurrentProcessInfo);

        // Apply slip condition
        mpRotationTool->Rotate(RHS_Contribution, rCurrentElement.GetGeometry());
        mpRotationTool->ApplySlipCondition(RHS_Contribution, rCurrentElement.GetGeometry());

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Get Equation Id
        rCurrentCondition.EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rCurrentCondition.CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentCondition.CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->CombineLHSContributions(LHS_Contribution,Mass,Damp,rCurrentProcessInfo);
        this->AddDynamicRHSContribution<Kratos::Condition>(rCurrentCondition,RHS_Contribution,Mass,rCurrentProcessInfo);

        // Apply slip condition
        mpRotationTool->Rotate(LHS_Contribution, RHS_Contribution, rCurrentCondition.GetGeometry());
        mpRotationTool->ApplySlipCondition(LHS_Contribution, RHS_Contribution, rCurrentCondition.GetGeometry());

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Condition &rCurrentCondition,
        LocalSystemVectorType &RHS_Contribution,
        Element::EquationIdVectorType &rEquationId,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Get Equation Id
        rCurrentCondition.EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        rCurrentCondition.CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentCondition.CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->AddDynamicRHSContribution<Kratos::Condition>(rCurrentCondition,RHS_Contribution,Mass,rCurrentProcessInfo);

        // Apply slip condition
        mpRotationTool->Rotate(RHS_Contribution, rCurrentCondition.GetGeometry());
        mpRotationTool->ApplySlipCondition(RHS_Contribution, rCurrentCondition.GetGeometry());

        KRATOS_CATCH("")
    }

    /// Free memory allocated by this object.
    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

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
        std::stringstream buffer;
        buffer << "BDF2TurbulentScheme";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

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

    /// Calculate the coefficients for time iteration.
    /**
     * @param rCurrentProcessInfo ProcessInfo instance from the fluid ModelPart. Must contain DELTA_TIME and BDF_COEFFICIENTS variables.
     */
    void SetTimeCoefficients(ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        //calculate the BDF coefficients
        double Dt = rCurrentProcessInfo[DELTA_TIME];
        double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

        double Rho = OldDt / Dt;
        double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

        Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
        BDFcoeffs.resize(3, false);

        BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
        BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
        BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)

        KRATOS_CATCH("");
    }

    /// Update Dof values after a Newton-Raphson iteration.
    /**
     * @param rDofSet Container for the Degrees of freedom in the system
     * @param Dx Solution of the linear system
     */
    virtual void UpdateDofs(
        DofsArrayType& rDofSet,
        TSystemVectorType& Dx)
    {
        KRATOS_TRY

        const int n_dof = rDofSet.size();

#pragma omp parallel for
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) {
            auto it_dof = rDofSet.begin() + i_dof;
            if (it_dof->IsFree()) {
                it_dof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx, it_dof->EquationId());
            }
        }

        KRATOS_CATCH("")
    }

    /// Update Dof values after a Newton-Raphson iteration
    /**
     * @param rModelPart fluid ModelPart
     * @param rBDFcoefs Time stepping coefficients for this iteration.
     */
    void UpdateAcceleration(
        ModelPart& rModelPart,
        const Vector& rBDFcoefs)
    {
        KRATOS_TRY

        const double Coef0 = rBDFcoefs[0];
        const double Coef1 = rBDFcoefs[1];
        const double Coef2 = rBDFcoefs[2];
        const int n_nodes = rModelPart.NumberOfNodes();

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = rModelPart.NodesBegin() + i_node;
            const auto& rVel0 = it_node->FastGetSolutionStepValue(VELOCITY);
            const auto& rVel1 = it_node->FastGetSolutionStepValue(VELOCITY,1);
            const auto& rVel2 = it_node->FastGetSolutionStepValue(VELOCITY,2);
            auto& rAcceleration = it_node->FastGetSolutionStepValue(ACCELERATION);

            rAcceleration = Coef0 * rVel0 + Coef1 * rVel1 + Coef2 * rVel2;
        }

        KRATOS_CATCH("")
    }

    void CombineLHSContributions(
        LocalSystemMatrixType& rLHS,
        LocalSystemMatrixType& rMass,
        LocalSystemMatrixType& rDamp,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const double Coef0 = rCurrentProcessInfo.GetValue(BDF_COEFFICIENTS)[0];
        if (rMass.size1() != 0) noalias(rLHS) += Coef0 * rMass;
        if (rDamp.size1() != 0) noalias(rLHS) += rDamp;
    }

    template<class TObject>
    void AddDynamicRHSContribution(
        TObject& rObject,
        LocalSystemVectorType& rRHS,
        LocalSystemMatrixType& rMass,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rMass.size1() != 0)
        {
            const Vector& rCoefs = rCurrentProcessInfo.GetValue(BDF_COEFFICIENTS);
            const auto& r_const_obj_ref = rObject;
            LocalSystemVectorType Acc;
            r_const_obj_ref.GetFirstDerivativesVector(Acc);
            Acc *= rCoefs[0];

            for(unsigned int n = 1; n < 3; ++n)
            {
                LocalSystemVectorType rVel;
                r_const_obj_ref.GetFirstDerivativesVector(rVel,n);
                noalias(Acc) += rCoefs[n] * rVel;
            }

            noalias(rRHS) -= prod(rMass,Acc);
        }
    }

    void FullProjection(ModelPart& rModelPart)
    {
        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // Initialize containers
        const int n_nodes = rModelPart.NumberOfNodes();
        const int n_elems = rModelPart.NumberOfElements();
        const array_1d<double,3> zero_vect = ZeroVector(3);
#pragma omp parallel for firstprivate(zero_vect)
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto ind = rModelPart.NodesBegin() + i_node;
            noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = zero_vect; // "x"
            ind->FastGetSolutionStepValue(DIVPROJ) = 0.0; // "x"
            ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0; // "Ml"
        }

        // Newton-Raphson parameters
        const double RelTol = 1e-4 * rModelPart.NumberOfNodes();
        const double AbsTol = 1e-6 * rModelPart.NumberOfNodes();
        const unsigned int MaxIter = 100;

        // iteration variables
        unsigned int iter = 0;
        array_1d<double,3> dMomProj = zero_vect;
        double dMassProj = 0.0;

        double RelMomErr = 1000.0 * RelTol;
        double RelMassErr = 1000.0 * RelTol;
        double AbsMomErr = 1000.0 * AbsTol;
        double AbsMassErr = 1000.0 * AbsTol;

        while( ( (AbsMomErr > AbsTol && RelMomErr > RelTol) || (AbsMassErr > AbsTol && RelMassErr > RelTol) ) && iter < MaxIter)
        {
            // Reinitialize RHS
#pragma omp parallel for firstprivate(zero_vect)
            for (int i_node = 0; i_node < n_nodes; ++i_node)
            {
                auto ind = rModelPart.NodesBegin() + i_node;
                noalias(ind->GetValue(ADVPROJ)) = zero_vect; // "b"
                ind->GetValue(DIVPROJ) = 0.0; // "b"
                ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0; // Reset because Calculate will overwrite it
            }

            // Reinitialize errors
            RelMomErr = 0.0;
            RelMassErr = 0.0;
            AbsMomErr = 0.0;
            AbsMassErr = 0.0;

            // Compute new values
            array_1d<double, 3 > output;
#pragma omp parallel for private(output)
            for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
                auto it_elem = rModelPart.ElementsBegin() + i_elem;
                it_elem->Calculate(SUBSCALE_VELOCITY, output, rCurrentProcessInfo);
            }

            rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
            rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);

            // Update iteration variables
#pragma omp parallel for
            for (int i_node = 0; i_node < n_nodes; ++i_node) {
                auto ind = rModelPart.NodesBegin() + i_node;
                const double Area = ind->FastGetSolutionStepValue(NODAL_AREA); // Ml dx = b - Mc x
                dMomProj = ind->GetValue(ADVPROJ) / Area;
                dMassProj = ind->GetValue(DIVPROJ) / Area;

                RelMomErr += sqrt( dMomProj[0]*dMomProj[0] + dMomProj[1]*dMomProj[1] + dMomProj[2]*dMomProj[2]);
                RelMassErr += fabs(dMassProj);

                auto& rMomRHS = ind->FastGetSolutionStepValue(ADVPROJ);
                double& rMassRHS = ind->FastGetSolutionStepValue(DIVPROJ);
                rMomRHS += dMomProj;
                rMassRHS += dMassProj;

                AbsMomErr += sqrt( rMomRHS[0]*rMomRHS[0] + rMomRHS[1]*rMomRHS[1] + rMomRHS[2]*rMomRHS[2]);
                AbsMassErr += fabs(rMassRHS);
            }

            if(AbsMomErr > 1e-10)
                RelMomErr /= AbsMomErr;
            else // If residual is close to zero, force absolute convergence to avoid division by zero errors
                RelMomErr = 1000.0;

            if(AbsMassErr > 1e-10)
                RelMassErr /= AbsMassErr;
            else
                RelMassErr = 1000.0;

            iter++;
        }

        KRATOS_INFO("BDF2TurbulentScheme") << "Performed OSS Projection in " << iter << " iterations" << std::endl;
    }

    void LumpedProjection(ModelPart& rModelPart)
    {
        const int n_nodes = rModelPart.NumberOfNodes();
        const int n_elems = rModelPart.NumberOfElements();
        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        const array_1d<double,3> zero_vect = ZeroVector(3);
#pragma omp parallel for firstprivate(zero_vect)
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto itNode = rModelPart.NodesBegin() + i_node;
            noalias(itNode->FastGetSolutionStepValue(ADVPROJ)) = zero_vect;
            itNode->FastGetSolutionStepValue(DIVPROJ) = 0.0;
            itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        }

        array_1d<double, 3 > Out;
#pragma omp parallel for private(Out)
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            auto itElem = rModelPart.ElementsBegin() + i_elem;
            itElem->Calculate(ADVPROJ, Out, rCurrentProcessInfo);
        }

        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
        rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

        // Correction for periodic conditions
        if (mrPeriodicIdVar.Key() != 0) {
            this->PeriodicConditionProjectionCorrection(rModelPart);
        }

        const double zero_tol = 1.0e-12;
#pragma omp parallel for firstprivate(zero_tol)
        for (int i_node = 0; i_node < n_nodes; ++i_node){
            auto iNode = rModelPart.NodesBegin() + i_node;
            if (iNode->FastGetSolutionStepValue(NODAL_AREA) < zero_tol) {
                iNode->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
            }
            const double Area = iNode->FastGetSolutionStepValue(NODAL_AREA);
            iNode->FastGetSolutionStepValue(ADVPROJ) /= Area;
            iNode->FastGetSolutionStepValue(DIVPROJ) /= Area;
        }

        KRATOS_INFO("BDF2TurbulentScheme") << "Computing OSS projections" << std::endl;
    }

    /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
     * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
     * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
     * 2- The non-historical containers are added across processes, communicating the right value from the condition owner to all partitions.\n
     * 3- The value on all periodic nodes is replaced by the one received in step 2.
     */
    void PeriodicConditionProjectionCorrection(ModelPart& rModelPart)
    {
        const int num_nodes = rModelPart.NumberOfNodes();
        const int num_conditions = rModelPart.NumberOfConditions();

        #pragma omp parallel for
        for (int i = 0; i < num_nodes; i++) {
            auto it_node = rModelPart.NodesBegin() + i;

            it_node->SetValue(NODAL_AREA,0.0);
            it_node->SetValue(ADVPROJ,ZeroVector(3));
            it_node->SetValue(DIVPROJ,0.0);
        }

        #pragma omp parallel for
        for (int i = 0; i < num_conditions; i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;

            if(it_cond->Is(PERIODIC)) {
                this->AssemblePeriodicContributionToProjections(it_cond->GetGeometry());
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

        #pragma omp parallel for
        for (int i = 0; i < num_nodes; i++) {
            auto it_node = rModelPart.NodesBegin() + i;
            this->CorrectContributionsOnPeriodicNode(*it_node);
        }
    }

    void AssemblePeriodicContributionToProjections(Geometry< Node >& rGeometry)
    {
        unsigned int nodes_in_cond = rGeometry.PointsNumber();

        double nodal_area = 0.0;
        array_1d<double,3> momentum_projection = ZeroVector(3);
        double mass_projection = 0.0;
        for ( unsigned int i = 0; i < nodes_in_cond; i++ )
        {
            auto& r_node = rGeometry[i];
            nodal_area += r_node.FastGetSolutionStepValue(NODAL_AREA);
            noalias(momentum_projection) += r_node.FastGetSolutionStepValue(ADVPROJ);
            mass_projection += r_node.FastGetSolutionStepValue(DIVPROJ);
        }

        for ( unsigned int i = 0; i < nodes_in_cond; i++ )
        {
            auto& r_node = rGeometry[i];
            /* Note that this loop is expected to be threadsafe in normal conditions,
             * since each node should belong to a single periodic link. However, I am
             * setting the locks for openmp in case that we try more complicated things
             * in the future (like having different periodic conditions for different
             * coordinate directions).
             */
            r_node.SetLock();
            r_node.GetValue(NODAL_AREA) = nodal_area;
            noalias(r_node.GetValue(ADVPROJ)) = momentum_projection;
            r_node.GetValue(DIVPROJ) = mass_projection;
            r_node.UnSetLock();
        }
    }

    void CorrectContributionsOnPeriodicNode(Node& rNode)
    {
        //TODO: This needs to be done in another manner as soon as we start using non-historical NODAL_AREA
        if (rNode.GetValue(NODAL_AREA) != 0.0) // Only periodic nodes will have a non-historical NODAL_AREA set.
        {
            rNode.FastGetSolutionStepValue(NODAL_AREA) = rNode.GetValue(NODAL_AREA);
            noalias(rNode.FastGetSolutionStepValue(ADVPROJ)) = rNode.GetValue(ADVPROJ);
            rNode.FastGetSolutionStepValue(DIVPROJ) = rNode.GetValue(DIVPROJ);
        }
    }

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

    /// Pointer to a turbulence model
    Process::Pointer mpTurbulenceModel = nullptr;

    RotationToolPointerType mpRotationTool = nullptr;

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

    const Kratos::Variable<int>& mrPeriodicIdVar;

//        ///@}
//        ///@name Serialization
//        ///@{
//
//        friend class Serializer;
//
//        virtual void save(Serializer& rSerializer) const
//        {
//            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
//            rSerializer.save("mpTurbulenceModel",mpTurbulenceModel);
//        }
//
//        virtual void load(Serializer& rSerializer)
//        {
//            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType );
//            rSerializer.load("mpTurbulenceModel",mpTurbulenceModel);
//        }

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
    BDF2TurbulentScheme & operator=(BDF2TurbulentScheme const& rOther)
    {}

    /// Copy constructor.
    BDF2TurbulentScheme(BDF2TurbulentScheme const& rOther)
    {}

    ///@}

}; // Class BDF2TurbulentScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<class TSparseSpace,class TDenseSpace>
inline std::istream& operator >>(std::istream& rIStream,BDF2TurbulentScheme<TSparseSpace,TDenseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpace,class TDenseSpace>
inline std::ostream& operator <<(std::ostream& rOStream,const BDF2TurbulentScheme<TSparseSpace,TDenseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_BDF2_TURBULENT_SCHEME_H_INCLUDED  defined
