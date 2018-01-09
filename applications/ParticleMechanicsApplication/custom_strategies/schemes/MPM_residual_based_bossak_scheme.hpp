//
//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:                    ilaria $
//   Date:                $Date:                   July 2015 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME )
#define  KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
//#include "solid_mechanics_application.h"
#include "particle_mechanics_application.h"
namespace Kratos
{
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
/*@} */

template<class TSparseSpace,  class TDenseSpace >
class MPMResidualBasedBossakScheme: public Scheme<TSparseSpace,TDenseSpace>
{
protected:

    struct GeneralAlphaMethod
    {

        double f;  //alpha Hilbert
        double m;  //alpha Bosssak

    };

    struct NewmarkMethod
    {

        double beta;
        double gamma;

        //system constants
        double c0;
        double c1;
        double c2;
        double c3;
        double c4;
        double c5;
        double c6;

        //static-dynamic parameter
        double static_dynamic;

    };


    struct  GeneralMatrices
    {

        std::vector< Matrix > M;     //first derivative matrix  (usually mass matrix)
        std::vector< Matrix > D;     //second derivative matrix (usually damping matrix)

    };

    struct GeneralVectors
    {

        std::vector< Vector > v;    //velocity
        std::vector< Vector > a;    //acceleration
        std::vector< Vector > ap;   //previous acceleration

    };



public:


    /**@name Type Definitions */

    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( MPMResidualBasedBossakScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef typename BaseType::Pointer                     BaseTypePointer;

    /*@} */

    /**
     * Constructor.
     * The bossak method
     */
    MPMResidualBasedBossakScheme(ModelPart& grid_model_part, double rAlpham=0,double rDynamic=1)
        :Scheme<TSparseSpace,TDenseSpace>(),mr_grid_model_part(grid_model_part)
    {
        //For pure Newmark Scheme
        mAlpha.f= 0;
        mAlpha.m= rAlpham;

        mNewmark.beta= 0.25;//(1.0+mAlpha.f-mAlpha.m)*(1.0+mAlpha.f-mAlpha.m)*0.25;
        mNewmark.gamma= 0.5;//0.5+mAlpha.f-mAlpha.m;

        mNewmark.static_dynamic= rDynamic;

        //std::cout << " MECHANICAL SCHEME: The Bossak Time Integration Scheme [alpha_m= "<<mAlpha.m<<" beta= "<<mNewmark.beta<<" gamma= "<<mNewmark.gamma<<"]"<<std::endl;


        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(NumThreads);
        mMatrix.D.resize(NumThreads);

        mVector.v.resize(NumThreads);
        mVector.a.resize(NumThreads);
        mVector.ap.resize(NumThreads);


    }


    /** Copy Constructor.
     */
    MPMResidualBasedBossakScheme(MPMResidualBasedBossakScheme& rOther)
        :BaseType(rOther)
        ,mAlpha(rOther.mAlpha)
        ,mNewmark(rOther.mNewmark)
        ,mMatrix(rOther.mMatrix)
        ,mVector(rOther.mVector)
        ,mr_grid_model_part(rOther.mr_grid_model_part)
    {
    }


    /** Destructor.
     */
    virtual ~MPMResidualBasedBossakScheme
    () {}

    /*@} */
    /**@name Operators
     */
    /*@{ */


    /**
     * Clone
     */
    virtual BaseTypePointer Clone()
    {
        return BaseTypePointer( new MPMResidualBasedBossakScheme(*this) );
    }



    //***************************************************************************
    //***************************************************************************

    /**
     * Performing the update of the solution
     * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param r_model_part
     * @param rDofSet set of all primary variables
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY

        //std::cout << " Update " << std::endl;
        //update of displacement (by DOF)
        for (typename DofsArrayType::iterator i_dof = rDofSet.begin(); i_dof != rDofSet.end(); ++i_dof)
        {
            if (i_dof->IsFree() )
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }


		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(r_model_part.Nodes().size()); ++iter)
		{
			
			auto i = r_model_part.NodesBegin() + iter;
			array_1d<double, 3 > & DeltaDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);


            array_1d<double, 3 > & CurrentVelocity      = (i)->FastGetSolutionStepValue(VELOCITY, 0);
            array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);

            array_1d<double, 3 > & CurrentAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 0);
            array_1d<double, 3 > & PreviousAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION, 1);





            UpdateVelocityUpdate     (CurrentVelocity, DeltaDisplacement, PreviousVelocity);

            UpdateAcceleration (CurrentAcceleration, DeltaDisplacement, PreviousVelocity, PreviousAcceleration);
		}



        

        KRATOS_CATCH( "" )
    }


    //***************************************************************************
    //***************************************************************************

    //predicts the solution for the current step:
    // x = 0

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        //std::cout << " Prediction " << std::endl;
        //array_1d<double, 3 > DeltaDisplacement;

		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(r_model_part.Nodes().size()); ++iter)
		{
			
			auto i = r_model_part.NodesBegin() + iter;
			array_1d<double, 3 > & PreviousVelocity     = (i)->FastGetSolutionStepValue(VELOCITY, 1);
            array_1d<double, 3 > & PreviousDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT, 1);
            array_1d<double, 3 > & CurrentDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
            //array_1d<double, 3 > & ImposedDisplacement  = (i)->FastGetSolutionStepValue(IMPOSED_DISPLACEMENT);
            array_1d<double, 3 > & PreviousAcceleration  = (i)->FastGetSolutionStepValue(ACCELERATION, 1);

            


            if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
            {
                CurrentDisplacement[0] = 0.0;
            }
            else
            {
                CurrentDisplacement[0]  = PreviousDisplacement[0];// + ImposedDisplacement[0];
            }

            if (i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false)
            {
                CurrentDisplacement[1] = 0.0;
            }
            else
            {
                CurrentDisplacement[1]  = PreviousDisplacement[1];// + ImposedDisplacement[1];
            }


            if (i->HasDofFor(DISPLACEMENT_Z))
            {
                if (i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false)
                {
                    CurrentDisplacement[2] = 0.0;
                }
                else
                {
                    CurrentDisplacement[2]  = PreviousDisplacement[2];// + ImposedDisplacement[2];
                }
            }

            //(i)->FastGetSolutionStepValue(DISPLACEMENT_AUX) = CurrentDisplacement;


            if (i->HasDofFor(PRESSURE))
            {
                double& PreviousPressure    = (i)->FastGetSolutionStepValue(PRESSURE, 1);
                double& CurrentPressure     = (i)->FastGetSolutionStepValue(PRESSURE);

                if ((i->pGetDof(PRESSURE))->IsFixed() == false)
                    CurrentPressure = PreviousPressure;
                //CurrentPressure = 0.0;

            }



            //updating time derivatives

            array_1d<double, 3 > & CurrentVelocity       = (i)->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3 > & CurrentAcceleration   = (i)->FastGetSolutionStepValue(ACCELERATION);



            UpdateVelocityPredict     (CurrentVelocity, CurrentDisplacement, PreviousVelocity, PreviousAcceleration);

            UpdateAcceleration (CurrentAcceleration, CurrentDisplacement, PreviousVelocity, PreviousAcceleration);
			
		}


        

    }

    //***************************************************************************
    //***************************************************************************

    /**
    this is the place to initialize the elements.
    This is intended to be called just once when the strategy is initialized
     */
    void InitializeElements(ModelPart& rModelPart)
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin() + ElementPartition[k];
            ElementsArrayType::iterator ElemEnd = rModelPart.Elements().begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)
            {
                itElem->Initialize(); //function to initialize the element
            }

        }

        this->mElementsAreInitialized = true;



        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /**
    this is the place to initialize the conditions.
    This is intended to be called just once when the strategy is initialized
    */
    void InitializeConditions(ModelPart& rModelPart)
    {
        KRATOS_TRY

        if(this->mElementsAreInitialized==false)
            KRATOS_THROW_ERROR( std::logic_error, "Before initilizing Conditions, initialize Elements FIRST", "" )

            int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin() + ConditionPartition[k];
            ConditionsArrayType::iterator CondEnd = rModelPart.Conditions().begin() + ConditionPartition[k + 1];

            for (ConditionsArrayType::iterator itCond = CondBegin; itCond != CondEnd; itCond++)
            {
                itCond->Initialize(); //function to initialize the condition
            }

        }

        this->mConditionsAreInitialized = true;
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /**
     * initializes time step solution
     * only for reasons if the time step solution is restarted
     * @param r_model_part
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        ProcessInfo CurrentProcessInfo= r_model_part.GetProcessInfo();

        //LOOP OVER THE GRID NODES PERFORMED FOR CLEAR ALL NODAL INFORMATION


		#pragma omp parallel for
		for(int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
		{
			
			auto i = mr_grid_model_part.NodesBegin() + iter;
			if( (i)->SolutionStepsDataHas(NODAL_MOMENTUM) && (i)->SolutionStepsDataHas(NODAL_MASS) && (i)->SolutionStepsDataHas(NODAL_INERTIA))//&& (i)->SolutionStepsDataHas(NODAL_INTERNAL_FORCE) )
            {

                array_1d<double, 3 > & NodalMomentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                array_1d<double, 3 > & NodalInertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                double & NodalMass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                //double & NodalMPressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                double & NodalPressure = (i)->FastGetSolutionStepValue(PRESSURE,1);

                double & NodalDensity = (i)->FastGetSolutionStepValue(DENSITY);
                double & NodalAuxR = (i)->FastGetSolutionStepValue(AUX_R);
                array_1d<double, 3 > & NodalAuxRVel = (i)->FastGetSolutionStepValue(AUX_R_VEL);
                array_1d<double, 3 > & NodalAuxRAcc = (i)->FastGetSolutionStepValue(AUX_R_ACC);

                NodalMomentum.clear();
                NodalInertia.clear();
                NodalMass= 0.0;
                //NodalMPressure = 0.0;
                NodalPressure = 0.0;

                NodalDensity = 0.0;
                NodalAuxR = 0.0;
                NodalAuxRVel.clear();
                NodalAuxRAcc.clear();
                //std::cout<< "NodalDensity "<< (i)->FastGetSolutionStepValue(DENSITY)<<std::endl;
            }

            if((i)->SolutionStepsDataHas(DISPLACEMENT) && (i)->SolutionStepsDataHas(VELOCITY) && (i)->SolutionStepsDataHas(ACCELERATION) )
            {
                array_1d<double, 3 > & NodalDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);
                double & NodalPressure = (i)->FastGetSolutionStepValue(PRESSURE);
                array_1d<double, 3 > & NodalVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & NodalAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                array_1d<double, 3 > & DeltaNodalVelocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY);
                array_1d<double, 3 > & DeltaNodalAcceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION);
                //array_1d<double, 3 > & NodalStress = (i)->FastGetSolutionStepValue(NODAL_STRESSES);

                NodalDisplacement.clear();
                NodalPressure = 0.0;
                NodalVelocity.clear();
                NodalAcceleration.clear();
                DeltaNodalVelocity.clear();
                DeltaNodalAcceleration.clear();
                //NodalStress.clear();
            }
		}


        





        //double RatioNormVel = 1.0;
        //double RatioNormAcc = 1.0;
        //double RatioNormPres = 1.0;
        double NormVel = 1.0;
        double NormAcc = 1.0;
        double NormPres = 1.0;
        double NormDeltaVel = 1.0;
        double NormDeltaAcc = 1.0;
        double NormDeltaPres = 1.0;
        //double TolVel = 5.0e-4;
        //double TolAcc = 5.0e-4;
        //double TolPres = 1.0e-4;
        int ItNum = 1;
        //if (CurrentProcessInfo.Has(PRESSURE)==false)
        //{
            //RatioNormPres = 1.0e-5;

        //}
        //for (unsigned int i = 0; i<20; i++)//
        //while(RatioNormVel > TolVel || RatioNormAcc > TolAcc)
        //while (RatioNormVel > TolVel || RatioNormPres > TolPres)// && ItNum <1000)
        while (ItNum <2)
        {
            //std::cout<< "ItNum "<<ItNum<<std::endl;

//**************************************************************************************
            
            
            #pragma omp parallel for
			for( int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
			{
			
			auto i = mr_grid_model_part.NodesBegin() + iter;
            if( (i)->SolutionStepsDataHas(NODAL_MOMENTUM) && (i)->SolutionStepsDataHas(NODAL_MASS) && (i)->SolutionStepsDataHas(NODAL_INERTIA))//&& (i)->SolutionStepsDataHas(NODAL_INTERNAL_FORCE) )
                {

                    array_1d<double, 3 > & NodalMomentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    array_1d<double, 3 > & NodalInertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                    double & NodalMPressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);
                    array_1d<double, 3 > & DeltaNodalVelocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY,1);
                    array_1d<double, 3 > & DeltaNodalAcceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION,1);

                    double & NodalMass = (i)->FastGetSolutionStepValue(NODAL_MASS);
                    NodalMomentum.clear();
                    NodalInertia.clear();
                    NodalMPressure = 0.0;
                    DeltaNodalVelocity.clear();
                    DeltaNodalAcceleration.clear();

                    NodalMass = 0.0;
                }
			}
            
            
            //IterativeExtrapolation evaluate again the global nodal nodal momentum inertia in function of a auxiliary function
            //std::cout<<"BEFORE CALLING FOR THE INITIALIZE SOLUTION STEP OF THE ELEMENT AFTER THE FIRST TIME"<<std::endl;
            NormVel = 0.0;
            NormAcc = 0.0;
            NormPres = 0.0;
            NormDeltaVel = 0.0;
            NormDeltaAcc = 0.0;
            NormDeltaPres = 0.0;
            //TolVel = 1;
            //TolAcc = 1;
            int nodes_counter = 0;
            Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(r_model_part,A,Dx,b);


			#pragma omp parallel for
			for(int iter = 0; iter < static_cast<int>(mr_grid_model_part.Nodes().size()); ++iter)
			{
			
			auto i = mr_grid_model_part.NodesBegin() + iter;
			double & NodalMass     = (i)->FastGetSolutionStepValue(NODAL_MASS);
                
                double DeltaNodalPressure = 0.0;


                if (NodalMass > 1.0e-16 )//> 1.0e-18)
                {
                    

                    array_1d<double, 3 > & DeltaNodalVelocity = (i)->FastGetSolutionStepValue(AUX_VELOCITY,1);
                    array_1d<double, 3 > & DeltaNodalAcceleration = (i)->FastGetSolutionStepValue(AUX_ACCELERATION,1);

                    array_1d<double, 3 > & NodalMomentum     = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    array_1d<double, 3 > & NodalInertia    = (i)->FastGetSolutionStepValue(NODAL_INERTIA);
                    double & NodalMPressure = (i)->FastGetSolutionStepValue(NODAL_MPRESSURE);

                    array_1d<double, 3 > & NodalVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
                    array_1d<double, 3 > & NodalAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
                    double & NodalPressure = (i)->FastGetSolutionStepValue(PRESSURE,1);

                    
                    if (i->HasDofFor(PRESSURE))
                    {
                        DeltaNodalPressure = NodalMPressure/NodalMass;

                    }

                    if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == false)
                    {
                        DeltaNodalVelocity[0] = NodalMomentum[0]/NodalMass;
                        DeltaNodalAcceleration[0] = NodalInertia[0]/NodalMass;
                    }
                    else
                    {
                        DeltaNodalVelocity[0] = 0.0;
                        DeltaNodalAcceleration[0] = 0.0;
                        //DeltaNodalAcceleration[0] = NodalInertia[0]/NodalMass;

                    }
                    if ((i->pGetDof(DISPLACEMENT_Y))->IsFixed() == false)
                    {


                        DeltaNodalAcceleration[1] = NodalInertia[1]/NodalMass;
                        DeltaNodalVelocity[1] = NodalMomentum[1]/NodalMass;


                        
                    }
                    else
                    {
                        DeltaNodalVelocity[1] = 0.0;
                        DeltaNodalAcceleration[1] = 0.0;
                        //DeltaNodalVelocity[1] = NodalMomentum[1]/NodalMass;

                    }
                    if (i->HasDofFor(DISPLACEMENT_Z))
                    {
                        if ((i->pGetDof(DISPLACEMENT_Z))->IsFixed() == false)
                        {
                            DeltaNodalVelocity[2] = NodalMomentum[2]/NodalMass;
                            DeltaNodalAcceleration[2] = NodalInertia[2]/NodalMass;
                        }
                        else
                        {
                            DeltaNodalVelocity[2] = 0.0;
                            DeltaNodalAcceleration[2] = 0.0;
                            //DeltaNodalAcceleration[2] = NodalInertia[2]/NodalMass;

                        }
                    }
//************************************************************************************************************************************************************
                 
                    NodalVelocity += DeltaNodalVelocity;
                    NodalAcceleration += DeltaNodalAcceleration;

                    NodalPressure += DeltaNodalPressure;

                    


                    NormDeltaVel += (DeltaNodalVelocity[0]*DeltaNodalVelocity[0]+DeltaNodalVelocity[1]*DeltaNodalVelocity[1]+DeltaNodalVelocity[2]*DeltaNodalVelocity[2]);
                    NormDeltaAcc += (DeltaNodalAcceleration[0]*DeltaNodalAcceleration[0]+DeltaNodalAcceleration[1]*DeltaNodalAcceleration[1]+DeltaNodalAcceleration[2]*DeltaNodalAcceleration[2]);
                    NormDeltaPres += (DeltaNodalPressure * DeltaNodalPressure);

                    NormVel += (NodalVelocity[0]*NodalVelocity[0]+NodalVelocity[1]*NodalVelocity[1]+NodalVelocity[2]*NodalVelocity[2]);
                    NormAcc += (NodalAcceleration[0]*NodalAcceleration[0]+NodalAcceleration[1]*NodalAcceleration[1]+NodalAcceleration[2]*NodalAcceleration[2]);
                    NormPres += (NodalPressure * NodalPressure);
                    
                    ++nodes_counter;
                }
			
			
			}


            
            

            NormVel = sqrt(NormVel);
            NormAcc = sqrt(NormAcc);
            NormPres = sqrt(NormPres);

            NormDeltaVel = sqrt(NormDeltaVel);
            NormDeltaAcc = sqrt(NormDeltaAcc);
            NormDeltaPres = sqrt(NormDeltaPres);
            
            ++ItNum;
        }



        double DeltaTime = CurrentProcessInfo[DELTA_TIME];

        if (DeltaTime == 0)
            KRATOS_THROW_ERROR( std::logic_error, "detected delta_time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )


            //initializing Newmark constants
            mNewmark.c0 = ( 1.0 / (mNewmark.beta * DeltaTime * DeltaTime) );
        mNewmark.c1 = ( mNewmark.gamma / (mNewmark.beta * DeltaTime) );
        mNewmark.c2 = ( 1.0 / (mNewmark.beta * DeltaTime) );
        mNewmark.c3 = ( 0.5 / (mNewmark.beta) - 1.0 );
        mNewmark.c4 = ( (mNewmark.gamma / mNewmark.beta) - 1.0  );
        mNewmark.c5 = ( DeltaTime * 0.5 * ( ( mNewmark.gamma / mNewmark.beta ) - 2 ) );


        //std::cout<<" Newmark Variables "<<mNewmark.c0<<" "<<mNewmark.c1<<" "<<mNewmark.c2<<" "<<mNewmark.c3<<" "<<mNewmark.c4<<" "<<mNewmark.c5<<std::endl;

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************
    /**
    Function called once at the beginning of each solution step in the initialize solution step.
    The basic operations to be carried in there are the following:
    - evaluation of the auxiliar variables in terms of velocity and acceleration on the particles
    - evaluate again the global values of nodal momentum and global inertia on the nodes of the connectivities
     */
    void IterativeExtrapolation(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        //iterative extrapolation
        //reset all nodal momentum and nodal inertia values
        //for (ModelPart::NodeIterator i = mr_grid_model_part.NodesBegin();
        //i != mr_grid_model_part.NodesEnd(); ++i)
        //{
        //if( (i)->SolutionStepsDataHas(NODAL_MOMENTUM) && (i)->SolutionStepsDataHas(NODAL_MASS) && (i)->SolutionStepsDataHas(NODAL_INERTIA))//&& (i)->SolutionStepsDataHas(NODAL_INTERNAL_FORCE) )
        //{

        //array_1d<double, 3 > & NodalMomentum = (i)->FastGetSolutionStepValue(NODAL_MOMENTUM);
        //array_1d<double, 3 > & NodalInertia = (i)->FastGetSolutionStepValue(NODAL_INERTIA);

        ////double & NodalMass = (i)->FastGetSolutionStepValue(NODAL_MASS);
        //NodalMomentum.clear();
        //NodalInertia.clear();
        //}
        //}


        //ElementsArrayType& pElements = rModelPart.Elements();
        //ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        //{
        //((UpdatedLagrangian*)it) -> IterativeExtrapolation(CurrentProcessInfo);
        //}

        //ConditionsArrayType& pConditions = rModelPart.Conditions();
        //for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        //{
        //UpdatedLagrangian * myCond = (*it);
        //myCond -> IterativeExtrapolation(CurrentProcessInfo);
        //}
        KRATOS_CATCH("")
    }

    //***************************************************************************
    //***************************************************************************


    /**
    function called once at the end of a solution step, after convergence is reached if
    an iterative process is needed
     */
    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY
        //finalizes solution step for all of the elements

        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        
//*********************************************************************************************************************************************
        for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
                i != rModelPart.NodesEnd(); ++i)
        {
            array_1d<double, 3 > & NodalVelocity = (i)->FastGetSolutionStepValue(VELOCITY,0);
            array_1d<double, 3 > & NodalAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,0);
            
            //NODAL VELOCITY AND NODAL ACCELERATION ARE SET TO ZERO FOR FIXED NODE, BEFORE MAPPING FROM NODES TO PARTICLES
            //THE NODAL INFORMATION (I am not sure..)
            if ((i->pGetDof(DISPLACEMENT_X))->IsFixed() == true)
            {
                NodalVelocity[0] = 0.0;
                NodalAcceleration[0] = 0.0;

            }
            if ((i->pGetDof(DISPLACEMENT_Y))->IsFixed() == true)
            {
                NodalVelocity[1] = 0.0;
                NodalAcceleration[1] = 0.0;

            }
        }
//*********************************************************************************************************************************************
        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            ElementsArrayType::iterator ElementsEnd = rElements.begin() + ElementPartition[k + 1];
            
            for (ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {

                itElem->FinalizeSolutionStep(CurrentProcessInfo);

                
            }
            
        }


        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            ConditionsArrayType::iterator ConditionsEnd = rConditions.begin() + ConditionPartition[k + 1];

            for (ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                itCond->FinalizeSolutionStep(CurrentProcessInfo);
            }
        }
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinIteration(ModelPart& r_model_part,
                                   TSystemMatrixType& A,
                                   TSystemVectorType& Dx,
                                   TSystemVectorType& b)
    {
        KRATOS_TRY
        ElementsArrayType& pElements = r_model_part.Elements();
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (ElementsArrayType::iterator it = pElements.begin(); it != pElements.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }

        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (ConditionsArrayType::iterator it = pConditions.begin(); it != pConditions.end(); ++it)
        {
            (it) -> InitializeNonLinearIteration(CurrentProcessInfo);
        }
        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
                                      ProcessInfo& CurrentProcessInfo)
    {
        (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }


    //***************************************************************************
    //***************************************************************************

    void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
                                      ProcessInfo& CurrentProcessInfo)
    {
        (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
        for (ModelPart::NodeIterator i = mr_grid_model_part.NodesBegin();
                i != mr_grid_model_part.NodesEnd(); ++i)
        {
            if( (i)->SolutionStepsDataHas(EXTERNAL_FORCE) && (i)->SolutionStepsDataHas(INTERNAL_FORCE) )
            {
                array_1d<double, 3 > & ExternalForce = (i)->FastGetSolutionStepValue(EXTERNAL_FORCE);
                array_1d<double, 3 > & InternalForce = (i)->FastGetSolutionStepValue(INTERNAL_FORCE);
                ExternalForce.clear();
                InternalForce.clear();
            }
        }
    }

    //***************************************************************************
    //***************************************************************************

    /** this function is designed to be called in the builder and solver to introduce*/

    void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread],CurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        }
        //if ((rCurrentElement) ->Id() == 365)
        //{

        //std::cout<<" final LHS_Contribution in integratios scheme "<<LHS_Contribution<<std::endl;
        //std::cout<<" final RHS_Contribution "<<RHS_Contribution<<std::endl;
        //}

        //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {

        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);
        }

        (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }

        KRATOS_CATCH( "" )

    }

    //***************************************************************************
    //***************************************************************************

    /** functions totally analogous to the precedent but applied to
          the "condition" objects
    */
    void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {


        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current element
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);


        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToLHS  (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);
        }

        //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);


        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        int thread = OpenMPUtils::ThisThread();

        //Initializing the non linear iteration for the current condition
        //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

        //basic operations for the element considered
        (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        if(mNewmark.static_dynamic !=0)
        {

            (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

            (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

        }

        (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

        //adding the dynamic contributions (static is already included)

        if(mNewmark.static_dynamic !=0)
        {

            AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        }

        KRATOS_CATCH( "" )
    }

    //***************************************************************************
    //***************************************************************************

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    void GetElementalDofList(
        Element::Pointer rCurrentElement,
        Element::DofsVectorType& ElementalDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /** Function that returns the list of Degrees of freedom to be
    assembled in the system for a Given Element
     */
    void GetConditionDofList(
        Condition::Pointer rCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
    }

    //***************************************************************************
    //***************************************************************************

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY

        int err = Scheme<TSparseSpace, TDenseSpace>::Check(r_model_part);
        if(err!=0) return err;

        //check for variables keys
        //verify that the variables are correctly initialized
        if(DISPLACEMENT.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )
            if(VELOCITY.Key() == 0)
                KRATOS_THROW_ERROR( std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered", "" )
                if(ACCELERATION.Key() == 0)
                    KRATOS_THROW_ERROR( std::invalid_argument,"ACCELERATION has Key zero! (check if the application is correctly registered", "" )

                    //check that variables are correctly allocated
                    for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                            it!=r_model_part.NodesEnd(); it++)
                    {
                        if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
                            KRATOS_THROW_ERROR( std::logic_error, "DISPLACEMENT variable is not allocated for node ", it->Id() )
                            if (it->SolutionStepsDataHas(VELOCITY) == false)
                                KRATOS_THROW_ERROR( std::logic_error, "VELOCITY variable is not allocated for node ", it->Id() )
                                if (it->SolutionStepsDataHas(ACCELERATION) == false)
                                    KRATOS_THROW_ERROR( std::logic_error, "ACCELERATION variable is not allocated for node ", it->Id() )
                                }

        //check that dofs exist
        for(ModelPart::NodesContainerType::iterator it=r_model_part.NodesBegin();
                it!=r_model_part.NodesEnd(); it++)
        {
            if(it->HasDofFor(DISPLACEMENT_X) == false)
                KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_X dof on node ",it->Id() )
                if(it->HasDofFor(DISPLACEMENT_Y) == false)
                    KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Y dof on node ",it->Id() )
                    if(it->HasDofFor(DISPLACEMENT_Z) == false)
                        KRATOS_THROW_ERROR( std::invalid_argument,"missing DISPLACEMENT_Z dof on node ",it->Id() )
                    }


        //check for admissible value of the AlphaBossak
        if(mAlpha.m > 0.0 || mAlpha.m < -0.3)
            KRATOS_THROW_ERROR( std::logic_error,"Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is ", mAlpha.m )

            //check for minimum value of the buffer index
            //verify buffer size
            if (r_model_part.GetBufferSize() < 2)
                KRATOS_THROW_ERROR( std::logic_error, "insufficient buffer size. Buffer size should be greater than 2. Current size is", r_model_part.GetBufferSize() )


                return 0;
        KRATOS_CATCH( "" )
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

protected:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */

    GeneralAlphaMethod  mAlpha;
    NewmarkMethod       mNewmark;

    GeneralMatrices     mMatrix;
    GeneralVectors      mVector;

    ModelPart& mr_grid_model_part;

    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocityPredict(array_1d<double, 3 > & CurrentVelocity,
                                      const array_1d<double, 3 > & DeltaDisplacement,
                                      const array_1d<double, 3 > & PreviousVelocity,
                                      const array_1d<double, 3 > & PreviousAcceleration)
    {

        noalias(CurrentVelocity) =  (mNewmark.c1 * DeltaDisplacement - mNewmark.c4 * PreviousVelocity
                                     - mNewmark.c5 * PreviousAcceleration) * mNewmark.static_dynamic;

    }
    inline void UpdateVelocityUpdate(array_1d<double, 3 > & CurrentVelocity,
                                     const array_1d<double, 3 > & DeltaDisplacement,
                                     const array_1d<double, 3 > & PreviousVelocity)
    {

        noalias(CurrentVelocity) =  (mNewmark.c1 * DeltaDisplacement - mNewmark.c4 * PreviousVelocity) * mNewmark.static_dynamic;

    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    inline void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
                                   const array_1d<double, 3 > & DeltaDisplacement,
                                   const array_1d<double, 3 > & PreviousVelocity,
                                   const array_1d<double, 3 > & PreviousAcceleration)
    {

        noalias(CurrentAcceleration) =  (mNewmark.c0 * DeltaDisplacement - mNewmark.c2 * PreviousVelocity
                                         -  mNewmark.c3 * PreviousAcceleration) * mNewmark.static_dynamic;


    }


    //Elements:
    //****************************************************************************

    /**
    Atangent = M*c0 + D*c1 + K

     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {

        // adding mass contribution to the dynamic stiffness
        if (M.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += M * (1-mAlpha.m) * mNewmark.c0 * mNewmark.static_dynamic;

        }

        //adding  damping contribution
        if (D.size1() != 0) // if M matrix declared
        {
            noalias(LHS_Contribution) += D * (1-mAlpha.f) * mNewmark.c1 * mNewmark.static_dynamic;

        }

    }

    //Elements:
    //****************************************************************************

    /**
    bdyn = b - M*a - D*v

     */
    void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m) * mNewmark.static_dynamic ;

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread] * mNewmark.static_dynamic;



            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);


        }

        //adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mNewmark.static_dynamic ;

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }


    }


    //Conditions:
    //****************************************************************************

    void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        //adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha.m) * mNewmark.static_dynamic;

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha.m * mVector.ap[thread] * mNewmark.static_dynamic;

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        //adding damping contribution
        //damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            (mVector.v[thread]) *= mNewmark.static_dynamic ;

            noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
        }

    }

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
private:
    /**@name Static Member Variables */
    /*@{ */
    /*@} */
    /**@name Member Variables */
    /*@{ */
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
    /**@name Unaccessible methods */
    /*@{ */
}; /* Class MPMResidualBasedBossakScheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_MPM_RESIDUAL_BASED_BOSSAK_SCHEME defined */


