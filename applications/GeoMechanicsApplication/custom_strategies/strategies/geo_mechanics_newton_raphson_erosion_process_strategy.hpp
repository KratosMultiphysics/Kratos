// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall,
//                   Aron Noordam
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

// Application includes
#include "boost/range/adaptor/filtered.hpp"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "geo_mechanics_application_variables.h"

#include <chrono>
#include <ctime>
#include <iostream>
#include <tuple>

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicsNewtonRaphsonErosionProcessStrategy
    : public GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsNewtonRaphsonErosionProcessStrategy);

    using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;
    using filtered_elements =
        typename boost::range_detail::filtered_range<std::function<bool(Element*)>, std::vector<Element*>>;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;

    GeoMechanicsNewtonRaphsonErosionProcessStrategy(ModelPart&                    model_part,
                                                    typename TSchemeType::Pointer pScheme,
                                                    typename TLinearSolver::Pointer pNewLinearSolver,
                                                    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                                    typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                                    Parameters& rParameters,
                                                    int         MaxIterations          = 30,
                                                    bool        CalculateReactions     = false,
                                                    bool        ReformDofSetAtEachStep = false,
                                                    bool        MoveMeshFlag           = false)
        : GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              model_part,
              pScheme,
              pNewLinearSolver,
              pNewConvergenceCriteria,
              pNewBuilderAndSolver,
              rParameters,
              MaxIterations,
              CalculateReactions,
              ReformDofSetAtEachStep,
              MoveMeshFlag)
    {
        rank              = model_part.GetCommunicator().MyPID();
        mPipingIterations = rParameters["max_piping_iterations"].GetInt();
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
            << "Max Piping Iterations: " << mPipingIterations << std::endl;

        bool grow = true;

        // get piping elements
        std::vector<Element*> PipeElements = GetPipingElements();
        unsigned int          n_el         = PipeElements.size(); // number of piping elements

        // get initially open pipe elements
        unsigned int openPipeElements = this->InitialiseNumActivePipeElements(PipeElements);

        if (PipeElements.empty()) {
            KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                << "No Pipe Elements -> Finalizing Solution " << std::endl;
            GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::FinalizeSolutionStep();
            return;
        }
        // calculate max pipe height and pipe increment
        double amax = CalculateMaxPipeHeight(PipeElements);

        // continue this loop, while the pipe is growing in length
        while (grow && (openPipeElements < n_el)) {
            // todo: JDN (20220817) : grow not used.
            // bool Equilibrium = false;

            // get tip element and activate
            Element* tip_element = PipeElements.at(openPipeElements);
            openPipeElements += 1;
            tip_element->SetValue(PIPE_ACTIVE, true);

            // Get all open pipe elements
            std::function<bool(Element*)> filter = [](Element* i) {
                return i->Has(PIPE_ACTIVE) && i->GetValue(PIPE_ACTIVE);
            };
            filtered_elements OpenPipeElements = PipeElements | boost::adaptors::filtered(filter);
            KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                << "Number of Open Pipe Elements: " << boost::size(OpenPipeElements) << std::endl;

            // non-lin picard iteration, for deepening the pipe
            // Todo JDN (20220817):: Deal with Equilibrium redundancy
            // Equilibrium = check_pipe_equilibrium(OpenPipeElements, amax, mPipingIterations);
            check_pipe_equilibrium(OpenPipeElements, amax, mPipingIterations);

            // check if pipe should grow in length
            std::tie(grow, openPipeElements) =
                check_status_tip_element(openPipeElements, n_el, amax, PipeElements);

            // if n open elements is lower than total pipe elements, save pipe height current
            // growing iteration or reset to previous iteration in case the pipe should not grow.
            if (openPipeElements < n_el) {
                save_or_reset_pipe_heights(OpenPipeElements, grow);
            }
            // recalculate groundwater flow
            bool converged = Recalculate();

            // error check
            KRATOS_ERROR_IF_NOT(converged) << "Groundwater flow calculation failed to converge." << std::endl;
        }

        GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::FinalizeSolutionStep();
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();
        this->GetBuilderAndSolver()->Check(BaseType::GetModelPart());
        this->GetScheme()->Check(BaseType::GetModelPart());
        return 0;

        KRATOS_CATCH("")
    }

    //-----------------------------Get Piping Elements--------------------------------------

    std::vector<Element*> GetPipingElements()
    {
        ModelPart&            CurrentModelPart = this->GetModelPart();
        std::vector<Element*> PipeElements;
        double                PipeElementStartX;

        for (Element& element : CurrentModelPart.Elements()) {
            if (element.GetProperties().Has(PIPE_START_ELEMENT)) {
                PipeElements.push_back(&element);

                long unsigned int startElement = element.GetProperties()[PIPE_START_ELEMENT];

                KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                    << element.Id() << " " << startElement << std::endl;

                if (element.Id() == startElement) {
                    PipeElementStartX = element.GetGeometry().GetPoint(0)[0];
                }
            }
        }

        if (PipeElements.empty()) {
            return PipeElements;
        }

        // Get Maximum X Value in Pipe
        auto rightPipe = std::max_element(PipeElements.begin(), PipeElements.end(),
                                          [](const Element* a, const Element* b) {
            return a->GetGeometry().GetPoint(0)[0] < b->GetGeometry().GetPoint(0)[0];
        });

        // Get Minimum X Value in Pipe
        auto leftPipe = std::min_element(PipeElements.begin(), PipeElements.end(),
                                         [](const Element* a, const Element* b) {
            return a->GetGeometry().GetPoint(0)[0] < b->GetGeometry().GetPoint(0)[0];
        });

        double minX = leftPipe[0]->GetGeometry().GetPoint(0)[0];
        double maxX = rightPipe[0]->GetGeometry().GetPoint(0)[0];

        KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
            << minX << " " << maxX << " " << PipeElementStartX << std::endl;

        if ((minX != PipeElementStartX) && (maxX != PipeElementStartX)) {
            KRATOS_ERROR << "Unable to determine pipe direction (multiple directions possible) -  "
                            "Check PIPE_START_ELEMENT"
                         << std::endl;
        }

        if (minX == PipeElementStartX) {
            // Pipe Left -> Right
            sort(PipeElements.begin(), PipeElements.end(), [](const Element* lhs, const Element* rhs) {
                return lhs->GetGeometry().GetPoint(0)[0] < rhs->GetGeometry().GetPoint(0)[0];
            });
        } else {
            // Pipe Right -> Left
            sort(PipeElements.begin(), PipeElements.end(), [](const Element* lhs, const Element* rhs) {
                return lhs->GetGeometry().GetPoint(0)[0] > rhs->GetGeometry().GetPoint(0)[0];
            });
        }

        KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
            << "Number of Pipe Elements: " << PipeElements.size() << std::endl;
        for (const Element* pipeElement : PipeElements) {
            KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                << "PipeElementIDs (in order): " << pipeElement->Id() << std::endl;
        }

        return PipeElements;
    }

private:
    unsigned int mPipingIterations; /// This is used to calculate the pipingLength
    int          rank;
    double       small_pipe_height    = 1e-10;
    double       pipe_height_accuracy = small_pipe_height * 10;

    /// <summary>
    /// Initialises the number of open pipe elements. This value can be greater than 0 in a multi
    /// staged analysis.
    /// </summary>
    /// <param name="PipeElements"></param>
    /// <returns></returns>
    int InitialiseNumActivePipeElements(std::vector<Element*> PipeElements)
    {
        int nOpenElements = 0;

        for (Element* pipe_element : PipeElements) {
            if (pipe_element->GetValue(PIPE_ACTIVE)) {
                nOpenElements += 1;
            }
        }
        return nOpenElements;
    }

    /// <summary>
    /// Calculates the pipe particle diameter according to the modified Sellmeijer rule if selected.
    /// Else the pipe particle diameter is equal to the d70.
    /// </summary>
    /// <param name="Prop"></param>
    /// <returns></returns>
    double CalculateParticleDiameter(const PropertiesType& Prop)
    {
        double diameter;

        if (Prop[PIPE_MODIFIED_D]) diameter = 2.08e-4 * pow((Prop[PIPE_D_70] / 2.08e-4), 0.4);
        else diameter = Prop[PIPE_D_70];
        return diameter;
    }

    /// <summary>
    /// Calculates the maximum pipe height which the algorithm will allow
    /// </summary>
    /// <param name="pipe_elements"> vector of all pipe elements</param>
    /// <returns></returns>
    double CalculateMaxPipeHeight(std::vector<Element*> pipe_elements)
    {
        double max_diameter  = 0;
        double height_factor = 100;

        // loop over all elements
        for (Element* pipe_element : pipe_elements) {
            // calculate pipe particle diameter of pipe element
            PropertiesType prop              = pipe_element->GetProperties();
            double         particle_diameter = this->CalculateParticleDiameter(prop);

            // get maximum pipe particle diameter of all pipe elements
            if (particle_diameter > max_diameter) {
                max_diameter = particle_diameter;
            }
        }

        // max pipe height is maximum pipe particle diameter * a constant
        return max_diameter * height_factor;
    }

    /// <summary>
    /// Calculates the pipe height increment.
    /// </summary>
    /// <param name="max_pipe_height"> maximum allowed pipe height</param>
    /// <param name="n_steps"> number of non lineair piping iteration steps</param>
    /// <returns></returns>
    double CalculatePipeHeightIncrement(double max_pipe_height, const unsigned int n_steps)
    {
        return max_pipe_height / (n_steps - 1);
    }

    bool Recalculate()
    {
        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0 && rank == 0)
            << "Recalculating" << std::endl;
        // KRATOS_INFO_IF("PipingLoop") << "Recalculating" << std::endl;
        // ModelPart& CurrentModelPart = this->GetModelPart();
        // this->Clear();

        // Reset displacements to the initial (Assumes Water Pressure is the convergence criteria)
        /* block_for_each(CurrentModelPart.Nodes(), [&](Node& rNode) {
             auto dold = rNode.GetSolutionStepValue(WATER_PRESSURE, 1);
             rNode.GetSolutionStepValue(WATER_PRESSURE, 0) = dold;
             });*/

        GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::InitializeSolutionStep();
        GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Predict();
        return GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::SolveSolutionStep();
    }

    bool check_pipe_equilibrium(filtered_elements open_pipe_elements, double amax, unsigned int mPipingIterations)
    {
        bool         equilibrium = false;
        bool         converged   = true;
        unsigned int PipeIter    = 0;
        // todo: JDN (20220817) : grow not used.
        // bool grow = true;

        // calculate max pipe height and pipe increment
        double da = CalculatePipeHeightIncrement(amax, mPipingIterations);

        while (PipeIter < mPipingIterations && !equilibrium && converged) {
            // set equilibrium on true
            equilibrium = true;

            // perform a flow calculation and stop growing if the calculation doesn't converge
            converged = Recalculate();

            // todo: JDN (20220817) : grow not used.
            // if (!converged)
            //{
            //    grow = false;
            //}

            if (converged) {
                // Update depth of open piping Elements
                equilibrium = true;
                for (auto OpenPipeElement : open_pipe_elements) {
                    auto pElement = static_cast<SteadyStatePwPipingElement<2, 4>*>(OpenPipeElement);

                    // get open pipe element geometry and properties
                    auto& Geom = OpenPipeElement->GetGeometry();
                    auto& prop = OpenPipeElement->GetProperties();

                    // calculate equilibrium pipe height and get current pipe height
                    double eq_height = pElement->CalculateEquilibriumPipeHeight(
                        prop, Geom, OpenPipeElement->GetValue(PIPE_ELEMENT_LENGTH));
                    double current_height = OpenPipeElement->GetValue(PIPE_HEIGHT);

                    // set erosion on true if current pipe height is greater than the equilibrium height
                    if (current_height > eq_height) {
                        OpenPipeElement->SetValue(PIPE_EROSION, true);
                    }

                    // check this if statement, I don't understand the check for pipe erosion
                    if (((!OpenPipeElement->GetValue(PIPE_EROSION) || (current_height > eq_height)) &&
                         current_height < amax)) {
                        OpenPipeElement->SetValue(PIPE_HEIGHT, OpenPipeElement->GetValue(PIPE_HEIGHT) + da);
                        equilibrium = false;
                    }

                    // check if equilibrium height and current pipe heights are diverging, stop
                    // picard iterations if this is the case and set pipe height on zero
                    if (!OpenPipeElement->GetValue(PIPE_EROSION) && (PipeIter > 1) &&
                        ((eq_height - current_height) > OpenPipeElement->GetValue(DIFF_PIPE_HEIGHT))) {
                        equilibrium = true;
                        OpenPipeElement->SetValue(PIPE_HEIGHT, small_pipe_height);
                    }

                    // calculate difference between equilibrium height and current pipe height
                    OpenPipeElement->SetValue(DIFF_PIPE_HEIGHT, eq_height - current_height);
                }
                // increment piping iteration number
                PipeIter += 1;
            }
        }
        return equilibrium;
    }

    /// <summary>
    /// Checks the status of the tip element and checks if it should continue or stop growing
    /// </summary>
    /// <param name="n_open_elements"> number of open pipe elements</param>
    /// <param name="n_elements"> number of pipe elements</param>
    /// <param name="max_pipe_height"> maximum allowed pipe height</param>
    /// <param name="PipeElements"> vector of all pipe elements</param>
    /// <returns>tuple of grow bool and number of open pipe elements</returns>
    std::tuple<bool, int> check_status_tip_element(unsigned int          n_open_elements,
                                                   unsigned int          n_elements,
                                                   double                max_pipe_height,
                                                   std::vector<Element*> PipeElements)
    {
        bool grow = true;
        // check status of tip element, stop growing if pipe_height is zero or greater than maximum
        // pipe height or if all elements are open
        if (n_open_elements < n_elements) {
            Element* tip_element = PipeElements.at(n_open_elements - 1);
            double   pipe_height = tip_element->GetValue(PIPE_HEIGHT);

            if ((pipe_height > max_pipe_height + std::numeric_limits<double>::epsilon()) ||
                (pipe_height < pipe_height_accuracy)) {
                // stable element found; pipe length does not increase during current time step
                grow = false;
                tip_element->SetValue(PIPE_EROSION, false);
                tip_element->SetValue(PIPE_ACTIVE, false);
                n_open_elements -= 1;

                KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                    << "Number of Open Pipe Elements: " << n_open_elements << std::endl;
            }
        } else {
            grow = false;
        }
        return std::make_tuple(grow, n_open_elements);
    }

    /// <summary>
    /// Saves pipe heights if pipe grows, else reset pipe heights to previous grow step.
    /// </summary>
    /// <param name="open_pipe_elements"> open pipe elements</param>
    /// <param name="grow"> boolean to check if pipe grows</param>
    /// <returns></returns>
    void save_or_reset_pipe_heights(filtered_elements open_pipe_elements, bool grow)
    {
        for (Element* OpenPipeElement : open_pipe_elements) {
            if (grow) {
                OpenPipeElement->SetValue(PREV_PIPE_HEIGHT, OpenPipeElement->GetValue(PIPE_HEIGHT));
            } else {
                OpenPipeElement->SetValue(PIPE_HEIGHT, OpenPipeElement->GetValue(PREV_PIPE_HEIGHT));
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; // Class GeoMechanicsNewtonRaphsonStrategy

} // namespace Kratos