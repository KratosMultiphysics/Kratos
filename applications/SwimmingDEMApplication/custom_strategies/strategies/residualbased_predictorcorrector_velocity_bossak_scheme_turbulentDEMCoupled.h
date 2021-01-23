//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//


#if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_DEM_COUPLED_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_TURBULENT_DEM_COUPLED_SCHEME

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "containers/array_1d.h"
#include "utilities/parallel_utilities.h"
#include "utilities/dof_updater.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "processes/process.h"

// Application includes
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
        :ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize),
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
        :ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, DomainSize, rPeriodicIdVar),
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
        :ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, rSlipFlag),
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
        :ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>(NewAlphaBossak, MoveMeshStrategy, DomainSize, pTurbulenceModel),
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
    ~ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentDEMCoupled() override {}


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

    void CalculateFluidFraction(
        ModelPart& r_model_part,
        ProcessInfo& r_current_process_info)
    {

        double delta_time = r_current_process_info[DELTA_TIME];
        double delta_time_inv = 1.0 / delta_time;

        //for(int k = 0; k<static_cast<int>(r_model_part.Nodes().size()); k++)
        block_for_each(r_model_part.Nodes(), [&](Node<3>& rNode)
        {
        const double fluid_fraction = rNode.FastGetSolutionStepValue(FLUID_FRACTION);
        const double fluid_fraction_old = rNode.FastGetSolutionStepValue(FLUID_FRACTION_OLD);
        const double fluid_fraction_rate = delta_time_inv * (fluid_fraction - fluid_fraction_old);
        rNode.FastGetSolutionStepValue(FLUID_FRACTION_RATE) = fluid_fraction_rate;
        });
    }

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        ProcessInfo current_process_info = r_model_part.GetProcessInfo();

        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);
        this->CalculateFluidFraction(r_model_part, current_process_info);

    }


    //*************************************************************************************
    //*************************************************************************************

    void FinalizeSolutionStep(
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        //for(int k = 0; k<static_cast<int>(rModelPart.Nodes().size()); k++)
        block_for_each(r_model_part.Nodes(), [&](Node<3>& rNode)
        {
            rNode.FastGetSolutionStepValue(FLUID_FRACTION_OLD) = rNode.FastGetSolutionStepValue(FLUID_FRACTION);
        });
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<TSparseSpace, TDenseSpace>::FinalizeSolutionStep(r_model_part, A, Dx, b);
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

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME_DEM_COUPLED defined */