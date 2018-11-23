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

#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BDF_TURBULENT_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BDF_TURBULENT_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/openmp_utils.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "processes/process.h"

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

    /// BDF2 time scheme for the incompressible flow problem.
    /** This scheme implements update operations and the calculation of the BDF coefficients for variable time step sizes.
     *
     * WARNING: this scheme assumes that the element internally implements the BDF2 scheme and is hence NOT compatible with the
     * elements ASGS2D, ASGS3D, VMS, MonolithicWallConditon
     *
     * the compatible element so far is
     *   @see TwoFluidVMS
     *
     * note also that in the prediction step only the velocity, and NOT the pressure is extrapolated in time.
     */
    template<class TSparseSpace,
    class TDenseSpace //= DenseSpace<double>
    >
    class ResidualBasedPredictorCorrectorBDFSchemeTurbulent : public Scheme<TSparseSpace, TDenseSpace> {
    public:
        /**@name Type Definitions */
        /*@{ */

        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedPredictorCorrectorBDFSchemeTurbulent);

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
        ResidualBasedPredictorCorrectorBDFSchemeTurbulent(
            unsigned int DomainSize)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0) // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
        {}

        /** Constructor without a turbulence model
         */
        ResidualBasedPredictorCorrectorBDFSchemeTurbulent(
            unsigned int DomainSize,
            Variable<double>& rSlipVar)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mRotationTool(DomainSize,DomainSize+1,rSlipVar,0.0) // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
        {}

        /** Constructor with a turbulence model
         */
        ResidualBasedPredictorCorrectorBDFSchemeTurbulent(
            unsigned int DomainSize,
            Process::Pointer pTurbulenceModel)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mpTurbulenceModel(pTurbulenceModel),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0) // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs
        {}

        /** Destructor.
         */
        ~ResidualBasedPredictorCorrectorBDFSchemeTurbulent() override {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

        /**
                Performing the update of the solution.
         */
        //***************************************************************************

        void Update(ModelPart& r_model_part,
                            DofsArrayType& rDofSet,
                            TSystemMatrixType& A,
                            TSystemVectorType& Dv,
                            TSystemVectorType& b) override
        {
            KRATOS_TRY;

            mRotationTool.RotateVelocities(r_model_part);

            mpDofUpdater->UpdateDofs(rDofSet,Dv);

            mRotationTool.RecoverVelocities(r_model_part);

            KRATOS_CATCH("")
        }

        //***************************************************************************
        //predicts the solution at the current step as
        // v = vold
        void Predict(ModelPart& rModelPart,
                             DofsArrayType& rDofSet,
                             TSystemMatrixType& A,
                             TSystemVectorType& Dv,
                             TSystemVectorType& b) override
        {
            ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            if(Dt != 0.0 && OldDt != 0)
            {
                //estimate acceleration from velocity in the past and predict the future. Note that pressure is NOT predicted
                const ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();

                array_1d<double,3> dv;

                //in the next loop we do for each node
                //vn+1 = vn + dt*(vn - vn-1)/oldDt
                #pragma omp parallel for private(dv)
                for(int i=0; i< static_cast<int>(rModelPart.Nodes().size()); i++)
                {

                    ModelPart::NodesContainerType::iterator it = it_begin + i;




                    const array_1d<double,3>& aux = it->FastGetSolutionStepValue(VELOCITY,1);

                    noalias(dv) = aux;
                    noalias(dv) -= it->FastGetSolutionStepValue(VELOCITY,2);

                    array_1d<double,3>& v = it->FastGetSolutionStepValue(VELOCITY);

                    const double dt_ratio = Dt/OldDt;
                    if(it->IsFixed(VELOCITY_X) == false) v[0] = aux[0] + dt_ratio*dv[0];
                    if(it->IsFixed(VELOCITY_Y) == false) v[1] = aux[1] + dt_ratio*dv[1];
                    if(it->IsFixed(VELOCITY_Z) == false) v[2] = aux[2] + dt_ratio*dv[2];

                    //noalias(v) = aux;
                    //noalias( v ) += () * dv;
                }
            }
            else
            {
                if (rModelPart.GetCommunicator().MyPID() == 0)
                    std::cout << "predict is doing nothing since OldDt = " << OldDt << "and Dt = " << Dt << std::endl;
            }


            // if (rModelPart.GetCommunicator().MyPID() == 0)
            //     std::cout << "end of prediction" << std::endl;
        }


        //***************************************************************************

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
        void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                          LocalSystemMatrixType& LHS_Contribution,
                                          LocalSystemVectorType& RHS_Contribution,
                                          Element::EquationIdVectorType& EquationId,
                                          ProcessInfo& CurrentProcessInfo) override
        {
            KRATOS_TRY

            //Initializing the non linear iteration for the current element
            (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentElement)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentElement->GetGeometry());

            KRATOS_CATCH("")
        }

        void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                        LocalSystemVectorType& RHS_Contribution,
                                        Element::EquationIdVectorType& EquationId,
                                        ProcessInfo& CurrentProcessInfo) override
        {
            //Initializing the non linear iteration for the current element
            (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentElement)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            // If there is a slip condition, apply it on a rotated system of coordinates
            mRotationTool.Rotate(RHS_Contribution,rCurrentElement->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentElement->GetGeometry());


        }

        /** functions totally analogous to the precedent but applied to
        the "condition" objects
         */
        void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                            LocalSystemMatrixType& LHS_Contribution,
                                                            LocalSystemVectorType& RHS_Contribution,
                                                            Element::EquationIdVectorType& EquationId,
                                                            ProcessInfo& CurrentProcessInfo) override
        {
            KRATOS_TRY

            //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
            (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

            (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }

        void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                                          LocalSystemVectorType& RHS_Contribution,
                                                          Element::EquationIdVectorType& EquationId,
                                                          ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;

            //KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
            //Initializing the non linear iteration for the current condition
            (rCurrentCondition) -> InitializeNonLinearIteration(rCurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentCondition)->CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

            (rCurrentCondition)->EquationIdVector(EquationId,rCurrentProcessInfo);

            // Rotate contributions (to match coordinates for slip conditions)
            mRotationTool.Rotate(RHS_Contribution,rCurrentCondition->GetGeometry());
            mRotationTool.ApplySlipCondition(RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("");
        }
        //*************************************************************************************
        //*************************************************************************************

        void InitializeSolutionStep(ModelPart& r_model_part,
                                    TSystemMatrixType& A,
                                    TSystemVectorType& Dx,
                                    TSystemVectorType& b) override
        {
            ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();

            if (r_model_part.GetBufferSize() != 3)
                KRATOS_THROW_ERROR(std::logic_error, "wrong buffer size. Expects 3, currently: ", r_model_part.GetBufferSize());

            //calculate the BDF coefficients
            double Dt = rCurrentProcessInfo[DELTA_TIME];
            double OldDt = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            if(OldDt == 0.0)
                KRATOS_THROW_ERROR(std::logic_error,"found an OldDt = 0.0 in InitializeSolutionStep","");

            double Rho = OldDt / Dt;
            double TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho);

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            if(BDFcoeffs.size() != 3)
                BDFcoeffs.resize(3, false);

            BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho); //coefficient for step n+1 (3/2Dt if Dt is constant)
            BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0); //coefficient for step n (-4/2Dt if Dt is constant)
            BDFcoeffs[2] = TimeCoeff; //coefficient for step n-1 (1/2Dt if Dt is constant)



            Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);


        }
        //*************************************************************************************
        //*************************************************************************************

        void InitializeNonLinIteration(ModelPart& r_model_part,
                                               TSystemMatrixType& A,
                                               TSystemVectorType& Dx,
                                               TSystemVectorType& b) override
        {
            KRATOS_TRY

            if (mpTurbulenceModel != 0) // If not null
                mpTurbulenceModel->Execute();

            KRATOS_CATCH("")
        }

        void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //if orthogonal subscales are computed
            if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {
                KRATOS_INFO_IF("ResidualBasedPredictorCorrectorBDFSchemeTurbulent", rModelPart.GetCommunicator().MyPID() == 0)
                    << "Computing OSS projections" << std::endl;
                for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++) {

                    noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

                    ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

                    ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


                }//end of loop over nodes

                //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
                array_1d<double, 3 > output;


                for (typename ModelPart::ElementsContainerType::iterator elem = rModelPart.ElementsBegin(); elem != rModelPart.ElementsEnd(); elem++)
                {
                    elem->Calculate(ADVPROJ, output, CurrentProcessInfo);
                }

                rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
                rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
                rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);


                for (typename ModelPart::NodesContainerType::iterator ind = rModelPart.NodesBegin(); ind != rModelPart.NodesEnd(); ind++)
                {
                    if (ind->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
                    {
                        ind->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
                        //KRATOS_WATCH("*********ATTENTION: NODAL AREA IS ZERRROOOO************");
                    }
                    const double Area = ind->FastGetSolutionStepValue(NODAL_AREA);
                    ind->FastGetSolutionStepValue(ADVPROJ) /= Area;
                    ind->FastGetSolutionStepValue(DIVPROJ) /= Area;
                }
            }
        }

        void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
			ComputeReactions(rModelPart, A, Dx, b);
            //Element::EquationIdVectorType EquationId;
            //LocalSystemVectorType RHS_Contribution;
            //LocalSystemMatrixType LHS_Contribution;
            //ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //
            //ModelPart::NodeIterator itnodes_begin = rModelPart.NodesBegin();
            //const int nnodes = static_cast<int>(rModelPart.Nodes().size());
            //#pragma omp parallel for firstprivate(nnodes, itnodes_begin)
            //for(int i=0; i<nnodes; i++)
            //{
            //    ModelPart::NodeIterator itNode = itnodes_begin + i;
            //    (itNode->FastGetSolutionStepValue(REACTION)).clear();
            //}
            //
            //
            //ModelPart::ElementsContainerType::iterator itelem_begin = rModelPart.ElementsBegin();
            //const int nelems = static_cast<int>(rModelPart.Elements().size());
            // #pragma omp parallel for firstprivate(nelems, itelem_begin)
            //for(int i=0; i<nelems; i++)
            //{
            //    ModelPart::ElementsContainerType::iterator itElem = itelem_begin + i;
            //
            //    (itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
            //    (itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo); //TODO: call CalculateRHS instead

            //    GeometryType& rGeom = (itElem)->GetGeometry();
            //    const int NumNodes = static_cast<int>(rGeom.PointsNumber());
            //    unsigned int Dimension = rGeom.WorkingSpaceDimension();
            //    unsigned int index = 0;

            //
            //    for (int i = 0; i < NumNodes; i++)
            //    {
            //
            //        array_1d<double,3>& rReaction = rGeom[i].FastGetSolutionStepValue(REACTION);
            //        rGeom[i].SetLock();
            //        rReaction[0] -= RHS_Contribution[index++];
            //        rReaction[1] -= RHS_Contribution[index++];
            //        if (Dimension == 3) rReaction[2] -= RHS_Contribution[index++];
            //        rGeom[i].UnSetLock();
            //        index++; // skip pressure dof
            //
            //    }
            //}
            //
            //rModelPart.GetCommunicator().AssembleCurrentData(REACTION);
            //
            //#pragma omp parallel for firstprivate(nelems, itelem_begin)
            //for(int i=0; i<nelems; i++)
            //{
            //    ModelPart::ElementsContainerType::iterator itElem = itelem_begin + i;
            //    (itElem)->FinalizeSolutionStep(CurrentProcessInfo);
            //}

        }
		virtual void ComputeReactions(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b)
		{
			Element::EquationIdVectorType EquationId;
			LocalSystemVectorType RHS_Contribution;
			LocalSystemMatrixType LHS_Contribution;
			ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();


			ModelPart::NodeIterator itnodes_begin = rModelPart.NodesBegin();
			const int nnodes = static_cast<int>(rModelPart.Nodes().size());
#pragma omp parallel for firstprivate(nnodes, itnodes_begin)
			for (int i = 0; i<nnodes; i++)
			{
				ModelPart::NodeIterator itNode = itnodes_begin + i;
				(itNode->FastGetSolutionStepValue(REACTION)).clear();
			}


			ModelPart::ElementsContainerType::iterator itelem_begin = rModelPart.ElementsBegin();
			const int nelems = static_cast<int>(rModelPart.Elements().size());
#pragma omp parallel for firstprivate(nelems, itelem_begin)
			for (int i = 0; i<nelems; i++)
			{
				ModelPart::ElementsContainerType::iterator itElem = itelem_begin + i;

				(itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
				(itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo); //TODO: call CalculateRHS instead

				GeometryType& rGeom = (itElem)->GetGeometry();
				const int NumNodes = static_cast<int>(rGeom.PointsNumber());
				unsigned int Dimension = rGeom.WorkingSpaceDimension();
				unsigned int index = 0;


				for (int i = 0; i < NumNodes; i++)
				{

					array_1d<double, 3>& rReaction = rGeom[i].FastGetSolutionStepValue(REACTION);
					rGeom[i].SetLock();
					rReaction[0] -= RHS_Contribution[index++];
					rReaction[1] -= RHS_Contribution[index++];
					if (Dimension == 3) rReaction[2] -= RHS_Contribution[index++];
					rGeom[i].UnSetLock();
					index++; // skip pressure dof

				}
			}

			rModelPart.GetCommunicator().AssembleCurrentData(REACTION);

#pragma omp parallel for firstprivate(nelems, itelem_begin)
			for (int i = 0; i<nelems; i++)
			{
				ModelPart::ElementsContainerType::iterator itElem = itelem_begin + i;
				(itElem)->FinalizeSolutionStep(CurrentProcessInfo);
			}
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
		Process::Pointer mpTurbulenceModel;

		CoordinateTransformationUtils<LocalSystemMatrixType, LocalSystemVectorType, double> mRotationTool;

        /*@} */
        /**@name Protected Operators*/
        /*@{ */

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

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BDF_TURBULENT_SCHEME  defined */
