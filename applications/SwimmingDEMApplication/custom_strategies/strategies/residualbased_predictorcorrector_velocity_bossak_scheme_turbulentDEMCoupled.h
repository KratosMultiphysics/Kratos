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


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_DEM_COUPLED_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_DEM_COUPLED_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/dof_updater.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "processes/process.h"
#include "../FluidDynamicsApplication/custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent.h"

namespace Kratos {

    /**@name Kratos Globals */
    /*@{ */


    /*@} */
    /**@name Type Definitions */
    /*@{ */

    /*@} */


    /**@name  Enum's */
    /*@{ */


    /*@} */
    /**@name  Functions */
    /*@{ */



    /*@} */
    /**@name Kratos Classes */
    /*@{ */

    /// Bossak time scheme for the incompressible flow problem.
    /** This class provides a second order time scheme of the generalized-alpha Newmark
        family of methods. It also includes code required to implement slip conditions
        on the incompressible flow problem and provides the possibility of using a RANS
        model by passing a turbulence model as an argument to the constructor.
        This time scheme is intended to be used in combination with elements of type
        ASGS2D, ASGS3D, VMS or derived classes.
        To use the slip condition, set the SLIP flag on slip wall nodes. To use
        a wall law in combination with the slip condition, use MonolithicWallCondition to
        mesh the boundary
        @see ASGS2D, ASGS3D, VMS, MonolithicWallConditon
     */
    template<class TSparseSpace,
    class TDenseSpace //= DenseSpace<double>
    >
    class ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled : public ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace> {
    public:
        /**@name Type Definitions */
        /*@{ */

        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled);

        typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename Element::DofsVectorType DofsVectorType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

        typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

        typedef Element::GeometryType  GeometryType;


        /*@} */
        /**@name Life Cycle
         */
        /*@{ */

        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize)
        :
          ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize),
          mRotationTool(DomainSize,DomainSize+1,SLIP), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
            mGammaNewmark = 0.5 - mAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;


            //Allocate auxiliary memory
            int NumThreads = OpenMPUtils::GetNumThreads();
            mMass.resize(NumThreads);
            mDamp.resize(NumThreads);
            mvel.resize(NumThreads);
            macc.resize(NumThreads);
            maccold.resize(NumThreads);}


        /** Constructor without a turbulence model with periodic conditions
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled(
            double NewAlphaBossak,
            unsigned int DomainSize,
            const Variable<int>& rPeriodicIdVar)
        :
          ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, DomainSize, rPeriodicIdVar),
          mRotationTool(DomainSize,DomainSize+1,SLIP), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(rPeriodicIdVar)
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
            mGammaNewmark = 0.5 - mAlphaBossak;
            mMeshVelocity = 0.0;


            //Allocate auxiliary memory
            int NumThreads = OpenMPUtils::GetNumThreads();
            mMass.resize(NumThreads);
            mDamp.resize(NumThreads);
            mvel.resize(NumThreads);
            macc.resize(NumThreads);
            maccold.resize(NumThreads);
          }


        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Kratos::Flags& rSlipFlag)
        :
          ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, rSlipFlag),
          mRotationTool(DomainSize,DomainSize+1,rSlipFlag), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject())
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
            mGammaNewmark = 0.5 - mAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;


            //Allocate auxiliary memory
            int NumThreads = OpenMPUtils::GetNumThreads();
            mMass.resize(NumThreads);
            mDamp.resize(NumThreads);
            mvel.resize(NumThreads);
            macc.resize(NumThreads);
            maccold.resize(NumThreads);
          }

        /** Constructor with a turbulence model
         */
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled(
            double NewAlphaBossak,
            double MoveMeshStrategy,
            unsigned int DomainSize,
            Process::Pointer pTurbulenceModel)
        :
          ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, pTurbulenceModel),
          mRotationTool(DomainSize,DomainSize+1,SLIP), // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs
          mrPeriodicIdVar(Kratos::Variable<int>::StaticObject()),
          mpTurbulenceModel(pTurbulenceModel)
          {
            //default values for the Newmark Scheme
            mAlphaBossak = NewAlphaBossak;
            mBetaNewmark = 0.25 * pow((1.00 - mAlphaBossak), 2);
            mGammaNewmark = 0.5 - mAlphaBossak;
            mMeshVelocity = MoveMeshStrategy;


            //Allocate auxiliary memory
            int NumThreads = OpenMPUtils::GetNumThreads();
            mMass.resize(NumThreads);
            mDamp.resize(NumThreads);
            mvel.resize(NumThreads);
            macc.resize(NumThreads);
            maccold.resize(NumThreads);
          }

        /** Destructor.
         */
        ~ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled() override {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */
        /** this function is designed to be called in the builder and solver
        to introduce
        the selected time integration scheme. It "asks" the matrix needed to
        the element and
        performs the operations needed to introduce the seected time
        integration scheme.

          this function calculates at the same time the contribution to the
        LHS and to the RHS
          of the system
         */

        //*************************************************************************************
        //*************************************************************************************

        void CalculateFluidFraction(ModelPart& r_model_part, ProcessInfo& CurrentProcessInfo)
        {

            double DeltaTime = CurrentProcessInfo[DELTA_TIME];
            double delta_time_inv = 1.0 / DeltaTime;

            //#pragma omp parallel for
            for(int k = 0; k<static_cast<int>(r_model_part.Nodes().size()); k++)
            {
                auto itNode = r_model_part.NodesBegin() + k;
                const double & FluidFraction = (itNode)->FastGetSolutionStepValue(FLUID_FRACTION);
                const double & FluidFractionRate = delta_time_inv * ((itNode)->FastGetSolutionStepValue(FLUID_FRACTION) - (itNode)->FastGetSolutionStepValue(FLUID_FRACTION_OLD));
                (itNode)->FastGetSolutionStepValue(FLUID_FRACTION_RATE) = FluidFractionRate;
                (itNode)->SetValue(FLUID_FRACTION, FluidFraction);
                (itNode)->SetValue(FLUID_FRACTION_RATE, FluidFractionRate);

            }
        }

        void InitializeSolutionStep(ModelPart& r_model_part,
                                    TSystemMatrixType& A,
                                    TSystemVectorType& Dx,
                                    TSystemVectorType& b) override
        {
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);
            this->CalculateFluidFraction(r_model_part, CurrentProcessInfo);

        }


        //*************************************************************************************
        //*************************************************************************************

        void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
            Element::EquationIdVectorType EquationId;
            LocalSystemVectorType RHS_Contribution;
            LocalSystemMatrixType LHS_Contribution;
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            #pragma omp parallel for
            for(int k = 0; k<static_cast<int>(rModelPart.Nodes().size()); k++)
            {
                auto itNode = rModelPart.NodesBegin() + k;
                (itNode->FastGetSolutionStepValue(REACTION)).clear();

                // calculating relaxed acceleration
                const array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 0);
                const array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                const array_1d<double, 3> relaxed_acceleration = (1 - mAlphaBossak) * CurrentAcceleration
                                                                    + mAlphaBossak * OldAcceleration;

                (itNode)->SetValue(RELAXED_ACCELERATION, relaxed_acceleration);
            }

            //for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin(); itElem != rModelPart.Elements().ptr_end(); ++itElem)

            #pragma omp parallel for firstprivate(EquationId,RHS_Contribution,LHS_Contribution)
            for(int k = 0; k<static_cast<int>(rModelPart.Elements().size()); k++)
            {
                auto itElem = rModelPart.Elements().ptr_begin()+k;
                int thread_id = OpenMPUtils::ThisThread();

                (*itElem)->InitializeNonLinearIteration(CurrentProcessInfo);

                //basic operations for the element considered
                (*itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

                (*itElem)->CalculateMassMatrix(mMass[thread_id], CurrentProcessInfo);

                (*itElem)->CalculateLocalVelocityContribution(mDamp[thread_id], RHS_Contribution, CurrentProcessInfo);

                (*itElem)->EquationIdVector(EquationId, CurrentProcessInfo);

                //adding the dynamic contributions (statics is already included)
                this->AddDynamicsToLHS(LHS_Contribution, mDamp[thread_id], mMass[thread_id], CurrentProcessInfo);
                this->AddDynamicsToRHS(*(*itElem), RHS_Contribution, mDamp[thread_id], mMass[thread_id], CurrentProcessInfo);

                GeometryType& rGeom = (*itElem)->GetGeometry();
                unsigned int NumNodes = rGeom.PointsNumber();
                unsigned int Dimension = rGeom.WorkingSpaceDimension();

                unsigned int index = 0;
                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    auto& reaction = rGeom[i].FastGetSolutionStepValue(REACTION);

                    double& target_value0 = reaction[0];
                    const double& origin_value0 = RHS_Contribution[index++];
                    #pragma omp atomic
                    target_value0 -= origin_value0;

                    double& target_value1 = reaction[1];
                    const double& origin_value1 = RHS_Contribution[index++];
                    #pragma omp atomic
                    target_value1 -= origin_value1;

                    if (Dimension == 3)
                    {
                      double& target_value2 = reaction[2];
                      const double& origin_value2 = RHS_Contribution[index++];
                      #pragma omp atomic
                      target_value2 -= origin_value2;
                    }
            //        rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
             //          rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
            //        if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
                    index++; // skip pressure dof
                }
            }

            rModelPart.GetCommunicator().AssembleCurrentData(REACTION);

            #pragma omp parallel for
            for(int k = 0; k<static_cast<int>(rModelPart.Nodes().size()); k++)
            {
                auto itNode = rModelPart.NodesBegin() + k;
                (itNode)->FastGetSolutionStepValue(FLUID_FRACTION_OLD) = (itNode)->FastGetSolutionStepValue(FLUID_FRACTION);
            }

            // Base scheme calls FinalizeSolutionStep method of elements and conditions
            Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, A, Dx, b);
        }

        //************************************************************************************************
        //************************************************************************************************

        /// Free memory allocated by this object.
        void Clear() override
        {
            this->mpDofUpdater->Clear();
        }

        /*@} */
        /**@name Operations */
        /*@{ */


        /*@} */
        /**@name Access */
        /*@{ */


        /*@} */
        /**@name Inquiry */
        /*@{ */


        /*@} */
        /**@name Friends */
        /*@{ */


        /*@} */

    protected:
        /**@name Protected static Member Variables */
        /*@{ */


        /*@} */
        /**@name Protected member Variables */
        /*@{ */
        double mAlphaBossak;
        double mBetaNewmark;
        double mGammaNewmark;
        double mMeshVelocity;

        double ma0;
        double ma1;
        double ma2;
        double ma3;
        double ma4;
        double ma5;
        double mam;

        std::vector< Matrix > mMass;
        std::vector< Matrix > mDamp;
        std::vector< Vector > mvel;
        std::vector< Vector > macc;
        std::vector< Vector > maccold;

        /*@} */
        /**@name Protected Operators*/
        /*@{ */

        /** On periodic boundaries, the nodal area and the values to project need to take into account contributions from elements on
         * both sides of the boundary. This is done using the conditions and the non-historical nodal data containers as follows:\n
         * 1- The partition that owns the PeriodicCondition adds the values on both nodes to their non-historical containers.\n
         * 2- The non-historical containers are added across processes, communicating the right value from the condition owner to all partitions.\n
         * 3- The value on all periodic nodes is replaced by the one received in step 2.
         */

        //****************************************************************************

        /**
        Kdyn = am*M + D + a1*K
         */

        //****************************************************************************

        /// Add Bossak contributions from the inertial term to the RHS vector.
        /** This essentially performs bdyn = b - M*acc for the current element.
         *  Note that viscous/pressure contributions to the RHS are expected to be added by the element itself.
         *  @param[in] rCurrentElement The fluid element we are assembling.
         *  @param[in/out] rRHS_Contribution The right hand side term where the contribution will be added.
         *  @param[in] rD The elemental velocity/pressure LHS matrix.
         *  @param[in] rM The elemental acceleration LHS matrix.
         *  @param[in] rCurrentProcessInfo ProcessInfo instance for the containing ModelPart.
         */

        /*@} */
        /**@name Protected Operations*/
        /*@{ */


        /*@} */
        /**@name Protected  Access */
        /*@{ */


        /*@} */
        /**@name Protected Inquiry */
        /*@{ */


        /*@} */
        /**@name Protected LifeCycle */
        /*@{ */



        /*@} */

    private:
        /**@name Static Member Variables */
        /*@{ */


        /*@} */
        /**@name Member Variables */
        /*@{ */

        CoordinateTransformationUtils<LocalSystemMatrixType,LocalSystemVectorType,double> mRotationTool;

        const Variable<int>& mrPeriodicIdVar;

        Process::Pointer mpTurbulenceModel;

        typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

        /*@} */
        /**@name Private Operators*/
        /*@{ */

        /*@} */
        /**@name Private Operations*/
        /*@{ */


        /*@} */
        /**@name Private  Access */
        /*@{ */


        /*@} */
        /**@name Private Inquiry */
        /*@{ */


        /*@} */
        /**@name Un accessible methods */
        /*@{ */


        /*@} */

    }; /* Class Scheme */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */