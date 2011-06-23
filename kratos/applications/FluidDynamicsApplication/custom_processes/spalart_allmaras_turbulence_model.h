//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SPALART_ALLMARAS_H_INCLUDED )
#define  KRATOS_SPALART_ALLMARAS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "spaces/ublas_space.h"

namespace Kratos
{
    ///@addtogroup ApplicationNameApplication
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

    /// Short class definition.

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

        /// Default constructor.

        SpalartAllmarasTurbulenceModel(
                ModelPart& rmodel_part,
                typename TLinearSolver::Pointer pNewLinearSolver,
                unsigned int domain_size,
                double non_linear_tol,
                unsigned int max_it,
                bool reform_dofset,
                unsigned int time_order)
        : mr_model_part(rmodel_part), mdomain_size(domain_size), mtol(non_linear_tol), mmax_it(max_it), mtime_order(time_order),madapt_for_fractional_step(false)
        {
            //************************************************************************************************
            //check that the variables needed are in the model part DISTANCE, MOLECULAR_VISCOSITY, TURBLUENT_VISCOSITY, VELOCITY, MESH_VELOCITY, VISCOSITY,NODALA_AREA,TEMP_CONV_PROJ
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(DISTANCE)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", DISTANCE);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(VELOCITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", VELOCITY);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(MOLECULAR_VISCOSITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", MOLECULAR_VISCOSITY);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(TURBULENT_VISCOSITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", TURBULENT_VISCOSITY);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", MESH_VELOCITY);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(VISCOSITY)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", VISCOSITY);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", NODAL_AREA);
            if (!(rmodel_part.NodesBegin()->SolutionStepsDataHas(TEMP_CONV_PROJ)))
                KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", TEMP_CONV_PROJ);

            //************************************************************************************************
            //construct a new auxiliary model part
            mspalart_model_part.SetBufferSize(mr_model_part.GetBufferSize());
            mspalart_model_part.Nodes() = mr_model_part.Nodes();
            mspalart_model_part.SetProcessInfo(mr_model_part.pGetProcessInfo());
            mspalart_model_part.SetProperties(mr_model_part.pProperties());
            //            KRATOS_WATCH( mspalart_model_part.pGetProcessInfo());
            //            KRATOS_WATCH( mr_model_part.pGetProcessInfo()      );

            std::string ElementName;
            if (domain_size == 2)
                ElementName = std::string("SpalartAllmaras2D");
            else
                KRATOS_ERROR(std::logic_error, "Spalart Allmaras not yet implemented in 3D", "")

            const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

            //generating the elements
            for (ModelPart::ElementsContainerType::iterator iii = mr_model_part.ElementsBegin(); iii != mr_model_part.ElementsEnd(); iii++)
            {
                Properties::Pointer properties = iii->pGetProperties();
                Element::Pointer p_element = rReferenceElement.Create(iii->Id(), iii->GetGeometry(), properties);
                mspalart_model_part.Elements().push_back(p_element);
            }

            //************************************************************************************************
            //construct strategy

//            ProcessInfo& rCurrentProcessInfo = mspalart_model_part.GetProcessInfo();

            //initializing fractional velocity solution step
            typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
            typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                    (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

            bool CalculateReactions = false;
            bool CalculateNormDxFlag = true;

            typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;

            BuilderSolverTypePointer componentwise_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> > (pNewLinearSolver, TURBULENT_VISCOSITY));
            mpstep1 = typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (mspalart_model_part, pscheme, pNewLinearSolver, componentwise_build, CalculateReactions, reform_dofset, CalculateNormDxFlag));
            mpstep1->SetEchoLevel(0);
        }

        void Execute()
        {
            KRATOS_TRY

            //            KRATOS_WATCH( mspalart_model_part.pGetProcessInfo());
            //            KRATOS_WATCH( mr_model_part.pGetProcessInfo()      );
            //            KRATOS_WATCH( *(mspalart_model_part.pGetProcessInfo()) );
            //            KRATOS_WATCH( *(mr_model_part.pGetProcessInfo() )     );

            if(madapt_for_fractional_step == true)
            {
                if (!(mspalart_model_part.NodesBegin()->SolutionStepsDataHas(FRACT_VEL)))
                    KRATOS_ERROR(std::logic_error, "Variable is not in the model part:", FRACT_VEL);

                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(mspalart_model_part.Nodes().size()); i++)
                {
                    ModelPart::NodesContainerType::iterator it = mspalart_model_part.NodesBegin() + i;
                    it->FastGetSolutionStepValue(VELOCITY) = it->FastGetSolutionStepValue(FRACT_VEL);
                }
            }

            AuxSolve();

            //update viscosity on the nodes
            for (ModelPart::NodeIterator i = mspalart_model_part.NodesBegin();
                    i != mspalart_model_part.NodesEnd(); ++i)
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
                    const double cv1 = 7.1;

                    double xi = turbulent_viscosity / molecular_viscosity;
                    double fv1 = (xi * xi * xi) / (xi * xi * xi + cv1 * cv1 * cv1);

                    double viscosity = fv1 * turbulent_viscosity + molecular_viscosity;
    //                KRATOS_WATCH(viscosity)
    //                KRATOS_WATCH(turbulent_viscosity)
                    i->FastGetSolutionStepValue(VISCOSITY) = viscosity;
                }
            }

            KRATOS_CATCH("");
        }

        /// Destructor.

        virtual ~SpalartAllmarasTurbulenceModel()
        {
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
            KRATOS_TRY

            //update viscosity on the nodes
            for (ModelPart::NodeIterator i = mspalart_model_part.NodesBegin();
                    i != mspalart_model_part.NodesEnd(); ++i)
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
                    KRATOS_ERROR(std::logic_error,"unexpected isolated node. Wrong node has Id ",i->Id());

                if(distance > h_max*CDES)
                    i->FastGetSolutionStepValue(DISTANCE) = h_max;

            }

            KRATOS_CATCH("");
        }


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{


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
            buffer << "SpalartAllmarasTurbulenceModel";
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "SpalartAllmarasTurbulenceModel";
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
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


        ///@}

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{
        ModelPart& mr_model_part;
        ModelPart mspalart_model_part;
        unsigned int mdomain_size;
        double mtol;
        unsigned int mmax_it;
        bool madapt_for_fractional_step;
        unsigned int mtime_order;
        typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mpstep1;

        ///@}
        ///@name Private Operators
        ///@{

        ///@}
        ///@name Private Operations
        ///@{

        //*********************************************************************************
        //**********************************************************************

        double AuxSolve()
        {
            KRATOS_TRY

            //calculate the BDF coefficients
            ProcessInfo& rCurrentProcessInfo = mspalart_model_part.GetProcessInfo();
            double Dt = rCurrentProcessInfo[DELTA_TIME];

            if (mtime_order == 2)
            {
                if (mspalart_model_part.GetBufferSize() < 3)
                    KRATOS_ERROR(std::logic_error, "insufficient buffer size for BDF2, currently buffer size is ", mspalart_model_part.GetBufferSize())

                    double dt_old = rCurrentProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME];

                double rho = dt_old / Dt;
                double coeff = 1.0 / (Dt * rho * rho + Dt * rho);

                rCurrentProcessInfo[BDF_COEFFICIENTS].resize(3);
                Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
                BDFcoeffs[0] = coeff * (rho * rho + 2.0 * rho); //coefficient for step n+1
                BDFcoeffs[1] = -coeff * (rho * rho + 2.0 * rho + 1.0); //coefficient for step n
                BDFcoeffs[2] = coeff;
            } else
            {
                rCurrentProcessInfo[BDF_COEFFICIENTS].resize(2);
                Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
                BDFcoeffs[0] = 1.0 / Dt; //coefficient for step n+1
                BDFcoeffs[1] = -1.0 / Dt; //coefficient for step n
            }



            unsigned int iter = 0;
            double ratio;
            bool is_converged = false;
            double dT_norm = 0.0;
            double T_norm = 0.0;

            while (iter++ < mmax_it && is_converged == false)
            {
                rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
                dT_norm = mpstep1->Solve();
                T_norm = CalculateVarNorm();
                CalculateProjection();
//                KRATOS_WATCH(dT_norm)
//                KRATOS_WATCH(T_norm)

                ratio = 1.00;
                if (T_norm != 0.00)
                    ratio = dT_norm / T_norm;
                else
                {
                    std::cout << "Nu norm = " << T_norm << " dNu_norm = " << dT_norm << std::endl;
                }

                if (dT_norm < 1e-11)
                    ratio = 0; //converged


                if (ratio < mtol)
                    is_converged = true;

                std::cout << "   SA iter = " << iter << " ratio = " << ratio << std::endl;

            }

            return dT_norm;
            KRATOS_CATCH("")
        }


        //******************************************************************************************************
        //******************************************************************************************************
        ///calculation of temperature norm

        double CalculateVarNorm()
        {
            KRATOS_TRY;

            double norm = 0.00;

            ProcessInfo& rCurrentProcessInfo = mspalart_model_part.GetProcessInfo();



            for (ModelPart::NodeIterator i = mspalart_model_part.NodesBegin();
                    i != mspalart_model_part.NodesEnd(); ++i)
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

            ProcessInfo& rCurrentProcessInfo = mspalart_model_part.GetProcessInfo();

            //first of all set to zero the nodal variables to be updated nodally
            for (ModelPart::NodeIterator i = mspalart_model_part.NodesBegin();
                    i != mspalart_model_part.NodesEnd(); ++i)
            {
                if ((i->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
                {
                    (i)->FastGetSolutionStepValue(TEMP_CONV_PROJ) = 0.00;
                    (i)->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
                } else
                {
                    (i)->FastGetSolutionStepValue(NODAL_AREA) = 1.00;
                }
            }

            //add the elemental contributions for the calculation of the velocity
            //and the determination of the nodal area
            rCurrentProcessInfo[FRACTIONAL_STEP] = 2;
            for (ModelPart::ElementIterator i = mspalart_model_part.ElementsBegin();
                    i != mspalart_model_part.ElementsEnd(); ++i)
            {
                (i)->InitializeSolutionStep(rCurrentProcessInfo);
            }

            //solve nodally for the velocity
            for (ModelPart::NodeIterator i = mspalart_model_part.NodesBegin();
                    i != mspalart_model_part.NodesEnd(); ++i)
            {
                double& conv_proj = (i)->FastGetSolutionStepValue(TEMP_CONV_PROJ);
                double temp = 1.00 / (i)->FastGetSolutionStepValue(NODAL_AREA);
                conv_proj *= temp;
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

#endif // KRATOS_SPALART_ALLMARAS_H_INCLUDED  defined


