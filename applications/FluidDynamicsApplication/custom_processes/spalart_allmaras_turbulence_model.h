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
//                   Riccardo Rossi
//


#if !defined(KRATOS_SPALART_ALLMARAS_TURBULENCE_H_INCLUDED )
#define  KRATOS_SPALART_ALLMARAS_TURBULENCE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "processes/process.h"
#include "includes/cfd_variables.h"
#include "solving_strategies/strategies/solving_strategy.h"
//#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
// #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incremental_aitken_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"

// Application includes
#include "custom_utilities/periodic_condition_utilities.h"
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

/// An impelementation of the Spalart-Allmaras turbulence model for incompressible flows.
/** Detail class definition.
 */
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class SpalartAllmarasTurbulenceModel : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SpalartAllmarasTurbulenceModel
    KRATOS_CLASS_POINTER_DEFINITION(SpalartAllmarasTurbulenceModel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for the Spalart-Allmaras turbulence model.
    /**
      * @param rModelPart ModelPart for the flow problem
      * @param pLinearSolver Pointer to the linear solver to use in the solution of the viscosity transport problem
      * @param DomainSize Spatial dimension of the problem (2 or 3)
      * @param NonLinearTol Relative tolerance for the turbulent viscosity transport problem (convergence is checked using the norm of the residual)
      * @param MaxIter Maximum number of iterations for the solution of the viscosity transport problem
      * @param ReformDofSet True if the degrees of freedom change during the problem (for example due to remeshing) false otherwise
      * @param TimeOrder Order for time integration (1 - Backward Euler will be used, 2 - BDF2 method)
      */
    SpalartAllmarasTurbulenceModel(
        ModelPart& rModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        unsigned int DomainSize,
        double NonLinearTol,
        unsigned int MaxIter,
        bool ReformDofSet,
        unsigned int TimeOrder)
        : mr_model_part(rModelPart),
        mrspalart_model_part(rModelPart.GetOwnerModel().CreateModelPart("SpalartModelPart")),
        mdomain_size(DomainSize),
        mtol(NonLinearTol),
        mmax_it(MaxIter),
        mtime_order(TimeOrder),
        madapt_for_fractional_step(false)
    {
        //************************************************************************************************
        //check that the variables needed are in the model part
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(DISTANCE)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", DISTANCE);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", VELOCITY);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(MOLECULAR_VISCOSITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", MOLECULAR_VISCOSITY);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(TURBULENT_VISCOSITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", TURBULENT_VISCOSITY);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", MESH_VELOCITY);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(VISCOSITY)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", VISCOSITY);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_AREA)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", NODAL_AREA);
        if (!(rModelPart.NodesBegin()->SolutionStepsDataHas(TEMP_CONV_PROJ)))
            KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", TEMP_CONV_PROJ);

        if (mr_model_part.GetBufferSize() < 3)
            KRATOS_THROW_ERROR(std::logic_error, "insufficient buffer size for BDF2, currently buffer size is ", mr_model_part.GetBufferSize());

        //************************************************************************************************
        //construct a new auxiliary model part
        mrspalart_model_part.GetNodalSolutionStepVariablesList() = mr_model_part.GetNodalSolutionStepVariablesList();
        mrspalart_model_part.SetBufferSize(3);
        mrspalart_model_part.Nodes() = mr_model_part.Nodes();
        mrspalart_model_part.SetProcessInfo(mr_model_part.pGetProcessInfo());
        mrspalart_model_part.SetProperties(mr_model_part.pProperties());

        std::string ElementName;
        if (DomainSize == 2)
            ElementName = std::string("SpalartAllmaras2D");
        else
            ElementName = std::string("SpalartAllmaras3D");

        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        //generating the elements
        for (ModelPart::ElementsContainerType::iterator iii = mr_model_part.ElementsBegin(); iii != mr_model_part.ElementsEnd(); iii++)
        {
            Properties::Pointer properties = iii->pGetProperties();
            Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
            mrspalart_model_part.Elements().push_back(p_element);
        }

        // pointer types for the solution strategy construcion
        typedef typename Scheme< TSparseSpace, TDenseSpace >::Pointer SchemePointerType;
        typedef typename ConvergenceCriteria< TSparseSpace, TDenseSpace >::Pointer ConvergenceCriteriaPointerType;
        typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
        typedef typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer StrategyPointerType;

        // Solution scheme: Aitken iterations
        const double DefaultAitkenOmega = 1.0;
        SchemePointerType pScheme = SchemePointerType( new ResidualBasedIncrementalAitkenStaticScheme< TSparseSpace, TDenseSpace > (DefaultAitkenOmega) );
//         SchemePointerType pScheme = SchemePointerType( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > () );

        // Convergence criteria
        const double NearlyZero = 1.0e-20;
        ConvergenceCriteriaPointerType pConvCriteria = ConvergenceCriteriaPointerType( new ResidualCriteria<TSparseSpace,TDenseSpace>(NonLinearTol,NearlyZero) );

        // Builder and solver
        BuilderSolverTypePointer pBuildAndSolver = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> > (pLinearSolver, TURBULENT_VISCOSITY));

        // Strategy
        bool CalculateReactions = false;
        bool MoveMesh = false;

        mpSolutionStrategy = StrategyPointerType( new ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(mrspalart_model_part,pScheme,pLinearSolver,pConvCriteria,pBuildAndSolver,MaxIter,CalculateReactions,ReformDofSet,MoveMesh));
        mpSolutionStrategy->SetEchoLevel(0);
        mpSolutionStrategy->Check();
    }

    /// Destructor.

    ~SpalartAllmarasTurbulenceModel() override
    {
        Model& r_model = mrspalart_model_part.GetOwnerModel();
        r_model.DeleteModelPart("SpalartModelPart");
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /// Solve an iteration of the turbulent viscosity
    void Execute() override
    {
        KRATOS_TRY

        if(madapt_for_fractional_step == true)
        {
            if (!(mrspalart_model_part.NodesBegin()->SolutionStepsDataHas(FRACT_VEL)))
                KRATOS_THROW_ERROR(std::logic_error, "Variable is not in the model part:", FRACT_VEL);

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mrspalart_model_part.Nodes().size()); i++)
            {
                ModelPart::NodesContainerType::iterator it = mrspalart_model_part.NodesBegin() + i;
                it->FastGetSolutionStepValue(VELOCITY) = it->FastGetSolutionStepValue(FRACT_VEL);
            }
        }

        AuxSolve();

        //update viscosity on the nodes
        for (ModelPart::NodeIterator i = mrspalart_model_part.NodesBegin();
                i != mrspalart_model_part.NodesEnd(); ++i)
        {
            double molecular_viscosity = i->FastGetSolutionStepValue(MOLECULAR_VISCOSITY);
            double turbulent_viscosity = i->FastGetSolutionStepValue(TURBULENT_VISCOSITY);

            if(turbulent_viscosity < 0)
            {
                i->FastGetSolutionStepValue(TURBULENT_VISCOSITY) = 1e-9;
                i->FastGetSolutionStepValue(VISCOSITY) = molecular_viscosity;
            }
            else
            {
                const double cv1_3 = 7.1*7.1*7.1;

                double xi = turbulent_viscosity / molecular_viscosity;
                double xi_3 = xi*xi*xi;
                double fv1 = xi_3 / (xi_3 + cv1_3);

                double viscosity = fv1 * turbulent_viscosity + molecular_viscosity;
                i->FastGetSolutionStepValue(VISCOSITY) = viscosity;
            }
        }

        KRATOS_CATCH("");
    }

    void SetMaxIterations(unsigned int max_it)
    {
        KRATOS_TRY

        mmax_it = max_it;

        KRATOS_CATCH("");
    }

    void AdaptForFractionalStep()
    {
        KRATOS_TRY

        madapt_for_fractional_step = true;

        KRATOS_CATCH("");
    }

    void ActivateDES(double CDES)
    {
        KRATOS_TRY;

        mrspalart_model_part.GetProcessInfo()[C_DES] = CDES;
/*
        //update viscosity on the nodes
        for (ModelPart::NodeIterator i = mrspalart_model_part.NodesBegin();
                i != mrspalart_model_part.NodesEnd(); ++i)
        {
            double distance = i->FastGetSolutionStepValue(DISTANCE);
            const array_1d<double,3>& xc = i->Coordinates();

            double h_max = 0.0;

            //compute nodal h (by max edge size)
            WeakPointerVector<Node<3> >& neigbours = i->GetValue(NEIGHBOUR_NODES);
            for(WeakPointerVector<Node<3> >::iterator ineighb=neigbours.begin(); ineighb!=neigbours.end(); ineighb++)
            {
                array_1d<double,3> aux = ineighb->Coordinates();
                aux -=  xc;

                double h = norm_2(aux);

                if(h > h_max) h_max=h;
            }

            if(h_max == 0.0)
                KRATOS_THROW_ERROR(std::logic_error,"unexpected isolated node. Wrong node has Id ",i->Id());

            if(distance > h_max*CDES)
                i->FastGetSolutionStepValue(DISTANCE) = h_max*CDES;

        }*/

        KRATOS_CATCH("");
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
        buffer << "SpalartAllmarasTurbulenceModel";
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SpalartAllmarasTurbulenceModel";
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
    }


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

    ModelPart& mr_model_part;
    ModelPart& mrspalart_model_part;
    unsigned int mdomain_size;
    double mtol;
    unsigned int mmax_it;
    unsigned int mtime_order;
    bool madapt_for_fractional_step;
    typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mpSolutionStrategy;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    /// Protected constructor, initializing only the references (for derived classes)
    SpalartAllmarasTurbulenceModel(ModelPart& rModelPart)
        :
        Process(),
        mr_model_part(rModelPart),
        mrspalart_model_part(rModelPart.GetOwnerModel().CreateModelPart("SpalartModelPart"))
    {}

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //*********************************************************************************
    //**********************************************************************

    /*double*/
    void AuxSolve()
    {
        KRATOS_TRY

        //calculate the BDF coefficients
        ProcessInfo& rCurrentProcessInfo = mrspalart_model_part.GetProcessInfo();
        double Dt = rCurrentProcessInfo[DELTA_TIME];

        if (mtime_order == 2)
        {
            double dt_old = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

            double rho = dt_old / Dt;
            double coeff = 1.0 / (Dt * rho * rho + Dt * rho);

            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(3);
            BDFcoeffs[0] = coeff * (rho * rho + 2.0 * rho); //coefficient for step n+1
            BDFcoeffs[1] = -coeff * (rho * rho + 2.0 * rho + 1.0); //coefficient for step n
            BDFcoeffs[2] = coeff;
        }
        else
        {
            Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
            BDFcoeffs.resize(2);
            BDFcoeffs[0] = 1.0 / Dt; //coefficient for step n+1
            BDFcoeffs[1] = -1.0 / Dt; //coefficient for step n
        }



//            unsigned int iter = 0;
//            double ratio;
//            bool is_converged = false;
//            double dT_norm = 0.0;
//            double T_norm = 0.0;

        int current_fract_step = rCurrentProcessInfo[FRACTIONAL_STEP];
        rCurrentProcessInfo[FRACTIONAL_STEP] = 2;

        CalculateProjection();

        rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
        mpSolutionStrategy->Solve();

        rCurrentProcessInfo[FRACTIONAL_STEP] = current_fract_step;

//            while (iter++ < mmax_it && is_converged == false)
//            {
//                rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
//                dT_norm = mpSolutionStrategy->Solve();
//                T_norm = CalculateVarNorm();
//                CalculateProjection();
////                KRATOS_WATCH(dT_norm)
////                KRATOS_WATCH(T_norm)

//                ratio = 1.00;
//                if (T_norm != 0.00)
//                    ratio = dT_norm / T_norm;
//                else
//                {
//                    std::cout << "Nu norm = " << T_norm << " dNu_norm = " << dT_norm << std::endl;
//                }

//                if (dT_norm < 1e-11)
//                    ratio = 0; //converged


//                if (ratio < mtol)
//                    is_converged = true;

//                std::cout << "   SA iter = " << iter << " ratio = " << ratio << std::endl;

//            }

//            return dT_norm;
        KRATOS_CATCH("")
    }


    //******************************************************************************************************
    //******************************************************************************************************
    ///calculation of temperature norm

    double CalculateVarNorm()
    {
        KRATOS_TRY;

        double norm = 0.00;



        for (ModelPart::NodeIterator i = mrspalart_model_part.NodesBegin();
                i != mrspalart_model_part.NodesEnd(); ++i)
        {
            norm += pow(i->FastGetSolutionStepValue(TURBULENT_VISCOSITY), 2);
        }

        return sqrt(norm);

        KRATOS_CATCH("")
    }

    ///calculation of projection

    void CalculateProjection()
    {
        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo = mrspalart_model_part.GetProcessInfo();

        //first of all set to zero the nodal variables to be updated nodally
        for (ModelPart::NodeIterator i = mrspalart_model_part.NodesBegin();
                i != mrspalart_model_part.NodesEnd(); ++i)
        {
            (i)->FastGetSolutionStepValue(TEMP_CONV_PROJ) = 0.00;
            (i)->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
        }

        //add the elemental contributions for the calculation of the velocity
        //and the determination of the nodal area


        for (ModelPart::ElementIterator i = mrspalart_model_part.ElementsBegin();
                i != mrspalart_model_part.ElementsEnd(); ++i)
        {
            (i)->InitializeSolutionStep(rCurrentProcessInfo);
        }

        Communicator& rComm = mrspalart_model_part.GetCommunicator();

        rComm.AssembleCurrentData(NODAL_AREA);
        rComm.AssembleCurrentData(TEMP_CONV_PROJ);

        // Obtain nodal projection of the residual
        for (ModelPart::NodeIterator i = mrspalart_model_part.NodesBegin();
                i != mrspalart_model_part.NodesEnd(); ++i)
        {
            const double NodalArea = i->FastGetSolutionStepValue(NODAL_AREA);
            if(NodalArea > 0.0)
            {
                double& rConvProj = i->FastGetSolutionStepValue(TEMP_CONV_PROJ);
                rConvProj /= NodalArea;
            }
        }

        KRATOS_CATCH("")
    }


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

    SpalartAllmarasTurbulenceModel & operator=(SpalartAllmarasTurbulenceModel const& rOther)
    {
        return *this;
    }

    /// Copy constructor.

    SpalartAllmarasTurbulenceModel(SpalartAllmarasTurbulenceModel const& rOther)
        : mr_model_part(rOther.mr_model_part), mdomain_size(rOther.mdomain_size)
    {
    }


    ///@}

}; // Class SpalartAllmarasTurbulenceModel

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
inline std::istream & operator >>(std::istream& rIStream,
                                  SpalartAllmarasTurbulenceModel<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    return rIStream;
}

/// output stream function

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const SpalartAllmarasTurbulenceModel<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_SPALART_ALLMARAS_TURBULENCE_H_INCLUDED  defined


