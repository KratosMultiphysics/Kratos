//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_GEAR_SCHEME_H_INCLUDED )
#define  KRATOS_GEAR_SCHEME_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "solving_strategies/schemes/scheme.h"
#include "includes/define.h"
//#include "includes/serializer.h"
#include "includes/dof.h"
//#include "includes/variables.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "containers/pointer_vector_set.h"
#include "utilities/openmp_utils.h"


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
class GearScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GearScheme
    KRATOS_CLASS_POINTER_DEFINITION(GearScheme);
    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;
    typedef typename TSparseSpace::DataType TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Dof<TDataType> TDofType;
    typedef typename BaseType::DofsArrayType DofsArrayType;



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GearScheme()
        :
        Scheme<TSparseSpace, TDenseSpace>(),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}

    /// Constructor to use the formulation combined with a turbulence model.
    /**
     * The turbulence model is assumed to be implemented as a Kratos::Process.
     * The model's Execute() method wil be called at the start of each
     * non-linear iteration.
     * @param pTurbulenceModel pointer to the turbulence model
     */
    GearScheme(Process::Pointer pTurbulenceModel)
        :
        Scheme<TSparseSpace, TDenseSpace>(),
        mpTurbulenceModel(pTurbulenceModel),
        mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
    {}

    /// Constructor for periodic boundary conditions.
    /**
     * @param rPeriodicVar the variable used to store periodic pair indices.
     */
    GearScheme(const Kratos::Variable<int>& rPeriodicVar)
        :
        Scheme<TSparseSpace, TDenseSpace>(),
        mrPeriodicIdVar(rPeriodicVar)
    {}


    /// Destructor.
    ~GearScheme() override
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

        int ErrorCode = BaseType::Check(rModelPart);
        if (ErrorCode != 0) return ErrorCode;

//            const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // Check buffer size
        if (rModelPart.GetBufferSize() < 3)
            KRATOS_THROW_ERROR(std::logic_error, "GearScheme error: Insufficient buffer size for BDF2, should be at least 3, got ",rModelPart.GetBufferSize());

        // Check that all required variables were registered
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"TIME_STEP Key is 0. Check if all applications were correctly registered.","");
        if(BDF_COEFFICIENTS.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"BDF_COEFFICIENTS Key is 0. Check if all applications were correctly registered.","");
        if(OSS_SWITCH.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"OSS_SWITCH Key is 0. Check if all applications were correctly registered.","");

        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISPLACEMENT Key is 0. Check if all applications were correctly registered.","");
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if all applications were correctly registered.","");
        if(MESH_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if all applications were correctly registered.","");
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check if all applications were correctly registered.","");

//            // Check that the ModelPart's ProcessInfo contains the required variables
//            if(rCurrentProcessInfo.Has(DELTA_TIME) != true)
//                KRATOS_THROW_ERROR(std::invalid_argument,"No value of DELTA_TIME defined in ProcessInfo for a model part passed to GearScheme","");
//            if(rCurrentProcessInfo.Has(BDF_COEFFICIENTS) != true)
//                KRATOS_THROW_ERROR(std::invalid_argument,"No value of BDF_COEFFICIENTS defined in ProcessInfo for a model part passed to GearScheme","");
//            if(rCurrentProcessInfo.Has(OSS_SWITCH) != true)
//                KRATOS_THROW_ERROR(std::invalid_argument,"No value of OSS_SWITCH defined in ProcessInfo for a model part passed to GearScheme","");

        return 0;
        KRATOS_CATCH("");
    }

    /// Set the time iteration coefficients
    void InitializeSolutionStep(ModelPart& rModelPart,
                                        TSystemMatrixType& A,
                                        TSystemVectorType& Dx,
                                        TSystemVectorType& b) override
    {
        this->SetTimeCoefficients(rModelPart.GetProcessInfo());

        // Base function initializes elements and conditions
        BaseType::InitializeSolutionStep(rModelPart,A,Dx,b);

        // Recalculate mesh velocity (to account for variable time step)
        const double Dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        const double OldDt = rModelPart.GetProcessInfo().GetPreviousSolutionStepInfo(1)[DELTA_TIME];
        if(Dt != OldDt)
        {
            const Vector& BDFcoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];

            // OpenMP partition
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.NumberOfNodes(),NumThreads,NodePartition);

            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();
                ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k+1];

                for(ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                {
                    array_1d<double,3>& rMeshVel = itNode->FastGetSolutionStepValue(MESH_VELOCITY);
                    const array_1d<double,3>& rDisp0 = itNode->FastGetSolutionStepValue(DISPLACEMENT);
                    const array_1d<double,3>& rDisp1 = itNode->FastGetSolutionStepValue(DISPLACEMENT,1);
                    const array_1d<double,3>& rDisp2 = itNode->FastGetSolutionStepValue(DISPLACEMENT,2);

                    rMeshVel = BDFcoefs[0] * rDisp0 + BDFcoefs[1] * rDisp1 + BDFcoefs[2] * rDisp2;
                }
            }
        }
    }

    void InitializeNonLinIteration(ModelPart& rModelPart,
                                           TSystemMatrixType& A,
                                           TSystemVectorType& Dx,
                                           TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (mpTurbulenceModel != 0) mpTurbulenceModel->Execute();

        KRATOS_CATCH("")
    }

    void FinalizeNonLinIteration(ModelPart &rModelPart,
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
    void Predict(ModelPart& rModelPart,
                         DofsArrayType& rDofSet,
                         TSystemMatrixType& A,
                         TSystemVectorType& Dx,
                         TSystemVectorType& b) override
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfNodes(), NumThreads, NodePartition);

        const Vector& BDFcoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];
        return;

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
            ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k+1];

            for(ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                array_1d<double,3>& rVel0 = itNode->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double,3>& rVel1 = itNode->FastGetSolutionStepValue(VELOCITY,1);
                const array_1d<double,3>& rVel2 = itNode->FastGetSolutionStepValue(VELOCITY,2);
                array_1d<double,3>& rAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);

                // Predict velocities
                if(!itNode->IsFixed(VELOCITY_X))
                    rVel0[0] = 2.00 * rVel1[0] - rVel2[0];
                if(!itNode->IsFixed(VELOCITY_Y))
                    rVel0[1] = 2.00 * rVel1[1] - rVel2[1];
                if(!itNode->IsFixed(VELOCITY_Z))
                    rVel0[2] = 2.00 * rVel1[2] - rVel2[2];

                // Predict acceleration
                rAcceleration = BDFcoefs[0] * rVel0 + BDFcoefs[1] * rVel1 + BDFcoefs[2] * rVel2;
            }
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
    void Update(ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        TSystemMatrixType& A,
                        TSystemVectorType& Dx,
                        TSystemVectorType& b) override
    {
        KRATOS_TRY

        mpDofUpdater->UpdateDofs(rDofSet,Dx);

        const Vector& BDFCoefs = rModelPart.GetProcessInfo()[BDF_COEFFICIENTS];

        this->UpdateAcceleration(rModelPart,BDFCoefs);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
            LocalSystemMatrixType& LHS_Contribution,
            LocalSystemVectorType& RHS_Contribution,
            Element::EquationIdVectorType& rEquationId,
            ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Initialize element
        rCurrentElement->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Get Equation Id
        rCurrentElement->EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentElement->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rCurrentElement->CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentElement->CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->CombineLHSContributions(LHS_Contribution,Mass,Damp,rCurrentProcessInfo);
        this->AddDynamicRHSContribution<Kratos::Element>(rCurrentElement,RHS_Contribution,Mass,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                    LocalSystemVectorType& RHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Initialize element
        rCurrentElement->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Get Equation Id
        rCurrentElement->EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentElement->CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        rCurrentElement->CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentElement->CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->AddDynamicRHSContribution<Kratos::Element>(rCurrentElement,RHS_Contribution,Mass,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
            LocalSystemMatrixType& LHS_Contribution,
            LocalSystemVectorType& RHS_Contribution,
            Element::EquationIdVectorType& rEquationId,
            ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Initialize element
        rCurrentCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Get Equation Id
        rCurrentCondition->EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentCondition->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);
        rCurrentCondition->CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentCondition->CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->CombineLHSContributions(LHS_Contribution,Mass,Damp,rCurrentProcessInfo);
        this->AddDynamicRHSContribution<Kratos::Condition>(rCurrentCondition,RHS_Contribution,Mass,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }


    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
            LocalSystemVectorType& RHS_Contribution,
            Element::EquationIdVectorType& rEquationId,
            ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemMatrixType Mass;
        LocalSystemMatrixType Damp;

        // Initialize element
        rCurrentCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Get Equation Id
        rCurrentCondition->EquationIdVector(rEquationId,rCurrentProcessInfo);

        // Get matrix contributions
        rCurrentCondition->CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);
        rCurrentCondition->CalculateMassMatrix(Mass,rCurrentProcessInfo);
        rCurrentCondition->CalculateLocalVelocityContribution(Damp,RHS_Contribution,rCurrentProcessInfo);

        // Add the dynamic contributions to the local system using BDF2 coefficients
        this->AddDynamicRHSContribution<Kratos::Condition>(rCurrentCondition,RHS_Contribution,Mass,rCurrentProcessInfo);

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GearScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GearScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
    virtual void UpdateDofs(DofsArrayType& rDofSet,
                            TSystemVectorType& Dx)
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];

            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {
                if (itDof->IsFree())
                    itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx, itDof->EquationId());
            }
        }

        KRATOS_CATCH("")
    }

    /// Update Dof values after a Newton-Raphson iteration
    /**
     * @param rModelPart fluid ModelPart
     * @param rBDFcoefs Time stepping coefficients for this iteration.
     */
    void UpdateAcceleration(ModelPart& rModelPart,
                            const Vector& rBDFcoefs)
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector NodePartition;
        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfNodes(), NumThreads, NodePartition);

        const double Coef0 = rBDFcoefs[0];
        const double Coef1 = rBDFcoefs[1];
        const double Coef2 = rBDFcoefs[2];

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
            ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k+1];

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                const array_1d<double,3>& rVel0 = itNode->FastGetSolutionStepValue(VELOCITY);
                const array_1d<double,3>& rVel1 = itNode->FastGetSolutionStepValue(VELOCITY,1);
                const array_1d<double,3>& rVel2 = itNode->FastGetSolutionStepValue(VELOCITY,2);
                array_1d<double,3>& rAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);

                rAcceleration = Coef0 * rVel0 + Coef1 * rVel1 + Coef2 * rVel2;
            }
        }

        KRATOS_CATCH("")
    }

    void CombineLHSContributions(LocalSystemMatrixType& rLHS,
                                 LocalSystemMatrixType& rMass,
                                 LocalSystemMatrixType& rDamp,
                                 const ProcessInfo& rCurrentProcessInfo)
    {
        const double Coef0 = rCurrentProcessInfo.GetValue(BDF_COEFFICIENTS)[0];
        if (rMass.size1() != 0) noalias(rLHS) += Coef0 * rMass;
        if (rDamp.size1() != 0) noalias(rLHS) += rDamp;
    }

    template<class TObject>
    void AddDynamicRHSContribution(typename TObject::Pointer pObject,
                                   LocalSystemVectorType& rRHS,
                                   LocalSystemMatrixType& rMass,
                                   const ProcessInfo& rCurrentProcessInfo)
    {
        if (rMass.size1() != 0)
        {
            const Vector& rCoefs = rCurrentProcessInfo.GetValue(BDF_COEFFICIENTS);

            LocalSystemVectorType Acc;
            pObject->GetFirstDerivativesVector(Acc);
            Acc *= rCoefs[0];

            for(unsigned int n = 1; n < 3; ++n)
            {
                LocalSystemVectorType rVel;
                pObject->GetFirstDerivativesVector(rVel,n);
                noalias(Acc) += rCoefs[n] * rVel;
            }

            noalias(rRHS) -= prod(rMass,Acc);
        }
    }


    void FullProjection(ModelPart& rModelPart)
    {
        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        // Initialize containers
        for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++)
        {
            noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3); // "x"
            ind->FastGetSolutionStepValue(DIVPROJ) = 0.0; // "x"
            ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0; // "Ml"
        }

        // Newton-Raphson parameters
        const double RelTol = 1e-4 * rModelPart.NumberOfNodes();
        const double AbsTol = 1e-6 * rModelPart.NumberOfNodes();
        const unsigned int MaxIter = 100;

        // iteration variables
        unsigned int iter = 0;
        array_1d<double,3> dMomProj = ZeroVector(3);
        double dMassProj = 0.0;

        double RelMomErr = 1000.0 * RelTol;
        double RelMassErr = 1000.0 * RelTol;
        double AbsMomErr = 1000.0 * AbsTol;
        double AbsMassErr = 1000.0 * AbsTol;

        while( ( (AbsMomErr > AbsTol && RelMomErr > RelTol) || (AbsMassErr > AbsTol && RelMassErr > RelTol) ) && iter < MaxIter)
        {
            // Reinitialize RHS
            for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++)
            {
                noalias(ind->GetValue(ADVPROJ)) = ZeroVector(3); // "b"
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

            for (typename ModelPart::ElementsContainerType::iterator elem = rModelPart.ElementsBegin(); elem != rModelPart.ElementsEnd(); elem++)
            {
                elem->Calculate(SUBSCALE_VELOCITY, output, rCurrentProcessInfo);
            }

            rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
            rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);
            rModelPart.GetCommunicator().AssembleNonHistoricalData(ADVPROJ);

            // Update iteration variables
            for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++)
            {
                const double Area = ind->FastGetSolutionStepValue(NODAL_AREA); // Ml dx = b - Mc x
                dMomProj = ind->GetValue(ADVPROJ) / Area;
                dMassProj = ind->GetValue(DIVPROJ) / Area;

                RelMomErr += sqrt( dMomProj[0]*dMomProj[0] + dMomProj[1]*dMomProj[1] + dMomProj[2]*dMomProj[2]);
                RelMassErr += fabs(dMassProj);

                array_1d<double,3>& rMomRHS = ind->FastGetSolutionStepValue(ADVPROJ);
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

        if (rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "Performed OSS Projection in " << iter << " iterations" << std::endl;
    }

    void LumpedProjection(ModelPart& rModelPart)
    {
        const ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        for (typename ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
        {
            noalias(itNode->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);
            itNode->FastGetSolutionStepValue(DIVPROJ) = 0.0;
            itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
        }

        array_1d<double, 3 > Out;

        for (typename ModelPart::ElementsContainerType::iterator itElem = rModelPart.ElementsBegin(); itElem != rModelPart.ElementsEnd(); itElem++)
        {
            itElem->Calculate(ADVPROJ, Out, rCurrentProcessInfo);
        }

        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
        rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);

        // Correction for periodic conditions
        if (mrPeriodicIdVar.Key() != 0) {
            this->PeriodicConditionProjectionCorrection(rModelPart);
        }

        for (typename ModelPart::NodesContainerType::iterator iNode = rModelPart.NodesBegin(); iNode != rModelPart.NodesEnd(); iNode++)
        {
            if (iNode->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
            {
                iNode->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
            }
            const double Area = iNode->FastGetSolutionStepValue(NODAL_AREA);
            iNode->FastGetSolutionStepValue(ADVPROJ) /= Area;
            iNode->FastGetSolutionStepValue(DIVPROJ) /= Area;
        }

        if (rModelPart.GetCommunicator().MyPID() == 0)
            std::cout << "Performed OSS Projection" << std::endl;
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

    void AssemblePeriodicContributionToProjections(Geometry< Node<3> >& rGeometry)
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

    void CorrectContributionsOnPeriodicNode(Node<3>& rNode)
    {
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

    /// Poiner to a turbulence model
    Process::Pointer mpTurbulenceModel = nullptr;

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
    GearScheme & operator=(GearScheme const& rOther)
    {}

    /// Copy constructor.
    GearScheme(GearScheme const& rOther)
    {}

    ///@}

}; // Class GearScheme

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpace,class TDenseSpace>
inline std::istream& operator >>(std::istream& rIStream,GearScheme<TSparseSpace,TDenseSpace>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpace,class TDenseSpace>
inline std::ostream& operator <<(std::ostream& rOStream,const GearScheme<TSparseSpace,TDenseSpace>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_GEAR_SCHEME_H_INCLUDED  defined
