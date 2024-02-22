//
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_PORO_NEWMARK_DYNAMIC_U_PL_SCHEME )
#define  KRATOS_PORO_NEWMARK_DYNAMIC_U_PL_SCHEME

// Application includes
#include "custom_strategies/schemes/one-phase_flow/poro_newmark_quasistatic_U_Pl_scheme.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class PoroNewmarkDynamicUPlScheme : public PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( PoroNewmarkDynamicUPlScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    using PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::mDeltaTime;
    using PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::mBeta;
    using PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::mGamma;
    using PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::mTheta_u;
    using PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::mTheta_p;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoroNewmarkDynamicUPlScheme(double theta_u,
                                double theta_p,
                                double beta,
                                double gamma) :
                                PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>(
                                    theta_u,
                                    theta_p,
                                    beta,
                                    gamma)
    {
        // Here mBeta and mGamma are used for the GN22 scheme
        mBeta = beta;
        mGamma = gamma;
        // Here mTheta_u is not used
        mTheta_u = 1.0;

        //Allocate auxiliary memory
        int NumThreads = ParallelUtilities::GetNumThreads();
        mMassMatrix.resize(NumThreads);
        mAccelerationVector.resize(NumThreads);
        mDampingMatrix.resize(NumThreads);
        mVelocityVector.resize(NumThreads);
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoroNewmarkDynamicUPlScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        // Base class checks
        int ierr = PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::Check(r_model_part);
        if(ierr != 0) return ierr;

        //check for variables keys (verify that the variables are correctly initialized)
        if(ACCELERATION.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION Key is 0. Check if all applications were correctly registered.", "" )

        //check that variables are correctly allocated
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin(); it!=r_model_part.NodesEnd(); it++)
        {
            if(it->SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_THROW_ERROR( std::logic_error, "ACCELERATION variable is not allocated for node ", it->Id() )
        }

        // Check beta, gamma
        if(mBeta <= 0.0 || mGamma< 0.0)
            KRATOS_THROW_ERROR( std::invalid_argument,"Some of the scheme variables: beta or gamma has an invalid value ", "" )

        return ierr;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        PoroNewmarkQuasistaticUPlScheme<TSparseSpace,TDenseSpace>::Initialize(r_model_part);

        // Note: VELOCITY_COEFFICIENT is updated according to the GN22 scheme used here
        r_model_part.GetProcessInfo().SetValue(VELOCITY_COEFFICIENT,mGamma/(mBeta*mDeltaTime));

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        // Predict Displacements on free nodes and update Acceleration, Velocity and DtPressure

        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            array_1d<double,3>& CurrentDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);

            const array_1d<double,3>& PreviousDisplacement = itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            if (itNode -> IsFixed(ACCELERATION_X))
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime * PreviousVelocity[0] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[0] + mBeta * CurrentAcceleration[0] );
            }
            else if (itNode -> IsFixed(VELOCITY_X))
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[0]-PreviousVelocity[0])+PreviousVelocity[0]);
            }
            else if (itNode -> IsFixed(DISPLACEMENT_X) == false)
            {
                CurrentDisplacement[0] = PreviousDisplacement[0] + mDeltaTime * PreviousVelocity[0] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[0];
            }

            if (itNode -> IsFixed(ACCELERATION_Y))
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime * PreviousVelocity[1] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[1] + mBeta * CurrentAcceleration[1] );
            }
            else if (itNode -> IsFixed(VELOCITY_Y))
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[1]-PreviousVelocity[1])+PreviousVelocity[1]);
            }
            else if (itNode -> IsFixed(DISPLACEMENT_Y) == false)
            {
                CurrentDisplacement[1] = PreviousDisplacement[1] + mDeltaTime * PreviousVelocity[1] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[1];
            }

            // For 3D cases
            if (itNode -> HasDofFor(DISPLACEMENT_Z))
            {
                if (itNode -> IsFixed(ACCELERATION_Z))
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime * PreviousVelocity[2] + std::pow(mDeltaTime, 2) * ( ( 0.5 - mBeta) * PreviousAcceleration[2] + mBeta * CurrentAcceleration[2] );
                }
                else if (itNode -> IsFixed(VELOCITY_Z))
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime*(mBeta/mGamma*(CurrentVelocity[2]-PreviousVelocity[2])+PreviousVelocity[2]);
                }
                else if (itNode -> IsFixed(DISPLACEMENT_Z) == false)
                {
                    CurrentDisplacement[2] = PreviousDisplacement[2] + mDeltaTime * PreviousVelocity[2] + 0.5 * std::pow(mDeltaTime, 2) * PreviousAcceleration[2];
                }
            }

            noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(LIQUID_PRESSURE) - itNode->FastGetSolutionStepValue(LIQUID_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta_p*mDeltaTime)*(DeltaPressure - (1.0-mTheta_p)*mDeltaTime*PreviousDtPressure);
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread],rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],rCurrentProcessInfo);

        this->AddDynamicsToLHS (LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        this->AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],rCurrentProcessInfo);

        this->AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Note: this is in a parallel loop

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLeftHandSide(LHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(mMassMatrix[thread],rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(mDampingMatrix[thread],rCurrentProcessInfo);

        this->AddDynamicsToLHS (LHS_Contribution, mMassMatrix[thread], mDampingMatrix[thread], rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    std::vector< Matrix > mMassMatrix;
    std::vector< Vector > mAccelerationVector;
    std::vector< Matrix > mDampingMatrix;
    std::vector< Vector > mVelocityVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    inline void UpdateVariablesDerivatives(ModelPart& r_model_part) override
    {
        KRATOS_TRY

        //Update Acceleration, Velocity and DtPressure

        array_1d<double,3> DeltaDisplacement;
        double DeltaPressure;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for private(DeltaDisplacement,DeltaPressure)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            array_1d<double,3>& CurrentAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION);
            array_1d<double,3>& CurrentVelocity = itNode->FastGetSolutionStepValue(VELOCITY);
            noalias(DeltaDisplacement) = itNode->FastGetSolutionStepValue(DISPLACEMENT) - itNode->FastGetSolutionStepValue(DISPLACEMENT, 1);
            const array_1d<double,3>& PreviousAcceleration = itNode->FastGetSolutionStepValue(ACCELERATION, 1);
            const array_1d<double,3>& PreviousVelocity = itNode->FastGetSolutionStepValue(VELOCITY, 1);

            noalias(CurrentAcceleration) = 1.0/(mBeta*mDeltaTime*mDeltaTime)*(DeltaDisplacement - mDeltaTime*PreviousVelocity - (0.5-mBeta)*mDeltaTime*mDeltaTime*PreviousAcceleration);
            noalias(CurrentVelocity) = PreviousVelocity + (1.0-mGamma)*mDeltaTime*PreviousAcceleration + mGamma*mDeltaTime*CurrentAcceleration;

            double& CurrentDtPressure = itNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE);
            DeltaPressure = itNode->FastGetSolutionStepValue(LIQUID_PRESSURE) - itNode->FastGetSolutionStepValue(LIQUID_PRESSURE, 1);
            const double& PreviousDtPressure = itNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE, 1);

            CurrentDtPressure = 1.0/(mTheta_p*mDeltaTime)*(DeltaPressure - (1.0-mTheta_p)*mDeltaTime*PreviousDtPressure);
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,LocalSystemMatrixType& M,LocalSystemMatrixType& C,const ProcessInfo& rCurrentProcessInfo)
    {
        // adding mass contribution
        if (M.size1() != 0)
            noalias(LHS_Contribution) += 1.0/(mBeta*mDeltaTime*mDeltaTime)*M;

        // adding damping contribution
        if (C.size1() != 0)
            noalias(LHS_Contribution) += mGamma/(mBeta*mDeltaTime)*C;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void AddDynamicsToRHS(Element& rCurrentElement,LocalSystemVectorType& RHS_Contribution,
                            LocalSystemMatrixType& M,LocalSystemMatrixType& C,const ProcessInfo& rCurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement.GetSecondDerivativesVector(mAccelerationVector[thread], 0);

            noalias(RHS_Contribution) -= prod(M, mAccelerationVector[thread]);
        }

        //adding damping contribution
        if (C.size1() != 0)
        {
            rCurrentElement.GetFirstDerivativesVector(mVelocityVector[thread], 0);

            noalias(RHS_Contribution) -= prod(C, mVelocityVector[thread]);
        }
    }

}; // Class PoroNewmarkDynamicUPlScheme
}  // namespace Kratos

#endif // KRATOS_PORO_NEWMARK_DYNAMIC_U_PL_SCHEME defined
