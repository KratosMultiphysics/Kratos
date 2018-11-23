/*
==============================================================================
KratosFluidDynamicsApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#if !defined(KRATOS_PFEM2_MONOLITHIC_SLIP_SCHEME )
#define  KRATOS_PFEM2_MONOLITHIC_SLIP_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"

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

    template<class TSparseSpace,
    class TDenseSpace //= DenseSpace<double>
    >
    class PFEM2MonolithicSlipScheme : public Scheme<TSparseSpace, TDenseSpace> {
    public:
        /**@name Type Definitions */
        /*@{ */

        //typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;

        KRATOS_CLASS_POINTER_DEFINITION(PFEM2MonolithicSlipScheme);

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
        PFEM2MonolithicSlipScheme(unsigned int DomainSize)
        :
          Scheme<TSparseSpace, TDenseSpace>(),
          mRotationTool(DomainSize,DomainSize+1,IS_STRUCTURE,0.0) // Second argument is number of matrix rows per node: monolithic elements have velocity and pressure dofs.
	{
	}


        /** Destructor.
         */
        virtual ~PFEM2MonolithicSlipScheme() {}


        /*@} */
        /**@name Operators
         */
        /*@{ */

        /**
                Performing the update of the solution.
         */
        //***************************************************************************

        virtual void Update(ModelPart& r_model_part,
                            DofsArrayType& rDofSet,
                            TSystemMatrixType& A,
                            TSystemVectorType& Dv,
                            TSystemVectorType& b) override
        {
            KRATOS_TRY;

            mRotationTool.RotateVelocities(r_model_part);

            BasicUpdateOperations(r_model_part, rDofSet, A, Dv, b);

            mRotationTool.RecoverVelocities(r_model_part);

            KRATOS_CATCH("")
        }

        //***************************************************************************

        virtual void BasicUpdateOperations(ModelPart& rModelPart,
                                           DofsArrayType& rDofSet,
                                           TSystemMatrixType& A,
                                           TSystemVectorType& Dv,
                                           TSystemVectorType& b)
        {
            KRATOS_TRY

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector DofSetPartition;
            OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

            //update of velocity (by DOF)
            #pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                typename DofsArrayType::iterator DofSetBegin = rDofSet.begin() + DofSetPartition[k];
                typename DofsArrayType::iterator DofSetEnd = rDofSet.begin() + DofSetPartition[k + 1];

                for (typename DofsArrayType::iterator itDof = DofSetBegin; itDof != DofSetEnd; itDof++) {
                    if (itDof->IsFree()) {
                        itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dv, itDof->EquationId());
                    }
                }

            }

            KRATOS_CATCH("")
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
            //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
            //KRATOS_WATCH(LHS_Contribution);
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
        virtual void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
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
            //mRotationTool.Rotate(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());
            //mRotationTool.ApplySlipCondition(LHS_Contribution,RHS_Contribution,rCurrentCondition->GetGeometry());

            KRATOS_CATCH("")
        }

        virtual void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
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

        virtual void InitializeSolutionStep(ModelPart& r_model_part,
                                            TSystemMatrixType& A,
                                            TSystemVectorType& Dx,
                                            TSystemVectorType& b) override
        {
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

            double DeltaTime = CurrentProcessInfo[DELTA_TIME];

            if (DeltaTime == 0)
                KRATOS_THROW_ERROR(std::logic_error, "detected delta_time = 0 ... check if the time step is created correctly for the current model part", "");

        }
        //*************************************************************************************
        //*************************************************************************************

        virtual void InitializeNonLinIteration(ModelPart& r_model_part,
                                               TSystemMatrixType& A,
                                               TSystemVectorType& Dx,
                                               TSystemVectorType& b) override
        {
            KRATOS_TRY

            KRATOS_CATCH("")
        }

        virtual void FinalizeNonLinIteration(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
            /*
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            //if orthogonal subscales are computed

            if (CurrentProcessInfo[OSS_SWITCH] == 1.0) {
                if (rModelPart.GetCommunicator().MyPID() == 0)
                    std::cout << "Computing OSS projections" << std::endl;
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

                // Correction for periodic conditions
                this->PeriodicConditionProjectionCorrection(rModelPart);

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
	    */
        }

        void FinalizeSolutionStep(ModelPart &rModelPart, TSystemMatrixType &A, TSystemVectorType &Dx, TSystemVectorType &b) override
        {
	    /*
            Element::EquationIdVectorType EquationId;
            LocalSystemVectorType RHS_Contribution;
            LocalSystemMatrixType LHS_Contribution;
            ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

            for (ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
            {
                itNode->FastGetSolutionStepValue(REACTION_X,0) = 0.0;
                itNode->FastGetSolutionStepValue(REACTION_Y,0) = 0.0;
                itNode->FastGetSolutionStepValue(REACTION_Z,0) = 0.0;
            }

            for (ModelPart::ElementsContainerType::ptr_iterator itElem = rModelPart.Elements().ptr_begin(); itElem != rModelPart.Elements().ptr_end(); ++itElem)
            {
                (*itElem)->InitializeNonLinearIteration(CurrentProcessInfo);
                //KRATOS_WATCH(LHS_Contribution);
                //basic operations for the element considered
                (*itElem)->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

                (*itElem)->EquationIdVector(EquationId, CurrentProcessInfo);



                GeometryType& rGeom = (*itElem)->GetGeometry();
                unsigned int NumNodes = rGeom.PointsNumber();
                unsigned int Dimension = rGeom.WorkingSpaceDimension();
                unsigned int index = 0;

                for (unsigned int i = 0; i < NumNodes; i++)
                {
                    rGeom[i].FastGetSolutionStepValue(REACTION_X,0) -= RHS_Contribution[index++];
                    rGeom[i].FastGetSolutionStepValue(REACTION_Y,0) -= RHS_Contribution[index++];
                    if (Dimension == 3) rGeom[i].FastGetSolutionStepValue(REACTION_Z,0) -= RHS_Contribution[index++];
                    index++; // skip pressure dof
                }
            }

            rModelPart.GetCommunicator().AssembleCurrentData(REACTION);
            */
            // Base scheme calls FinalizeSolutionStep method of elements and conditions
            Scheme<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(rModelPart, A, Dx, b);
        }

        //************************************************************************************************
        //************************************************************************************************

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

        /*@} */
        /**@name Protected Operators*/
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
