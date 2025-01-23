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
#include "custom_elements/geo_steady_state_Pw_piping_element.h"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"

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
    using PropertiesType           = Properties;
    using NodeType                 = Node;
    using GeometryType             = Geometry<NodeType>;
    using SizeType                 = std::size_t;

    GeoMechanicsNewtonRaphsonErosionProcessStrategy(ModelPart&                    model_part,
                                                    typename TSchemeType::Pointer pScheme,
                                                    typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                                    typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                                    Parameters& rParameters,
                                                    int         MaxIterations          = 30,
                                                    bool        CalculateReactions     = false,
                                                    bool        ReformDofSetAtEachStep = false,
                                                    bool        MoveMeshFlag           = false)
        : GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
              model_part, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, rParameters, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        rank              = model_part.GetCommunicator().MyPID();
        mPipingIterations = rParameters["max_piping_iterations"].GetInt();
    }

    template <typename PipingElementType>
    std::optional<std::vector<PipingElementType*>> TryDownCastToPipingElement(const std::vector<Element*>& rPipeElements)
    {
        std::vector<PipingElementType*> result{rPipeElements.size()};
        std::transform(rPipeElements.begin(), rPipeElements.end(), result.begin(),
                       [](auto p_element) { return dynamic_cast<PipingElementType*>(p_element); });

        const auto number_of_piping_elements = static_cast<std::size_t>(std::count_if(
            result.begin(), result.end(), [](auto p_element) { return p_element != nullptr; }));
        if (number_of_piping_elements == 0) return std::nullopt;

        KRATOS_ERROR_IF(number_of_piping_elements != rPipeElements.size())
            << "Unexpected number of piping elements of type " << rPipeElements[0]->Info() << ": expected "
            << rPipeElements.size() << " but got " << number_of_piping_elements << std::endl;

        return result;
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
            << "Max Piping Iterations: " << mPipingIterations << std::endl;

        const auto piping_elements = GetPipingElements();
        if (piping_elements.empty()) {
            KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                << "No Pipe Elements -> Finalizing Solution " << std::endl;
            this->BaseClassFinalizeSolutionStep();
            return;
        }

        if (const auto piping_interface_elements =
                TryDownCastToPipingElement<SteadyStatePwPipingElement<2, 4>>(piping_elements)) {
            this->DetermineOpenPipingElements(piping_interface_elements.value());
            this->BaseClassFinalizeSolutionStep();
            return;
        }

        if (const auto piping_2D_line_elements =
                TryDownCastToPipingElement<GeoSteadyStatePwPipingElement<2, 2>>(piping_elements)) {
            this->DetermineOpenPipingElements(piping_2D_line_elements.value());
            this->BaseClassFinalizeSolutionStep();
            return;
        }

        if (const auto piping_3D_line_elements =
                TryDownCastToPipingElement<GeoSteadyStatePwPipingElement<3, 2>>(piping_elements)) {
            this->DetermineOpenPipingElements(piping_3D_line_elements.value());
            this->BaseClassFinalizeSolutionStep();
            return;
        }

        KRATOS_ERROR << "Unexpected type of piping element\n";
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

    std::vector<Element*> GetPipingElements()
    {
        ModelPart&            r_current_model_part = this->GetModelPart();
        std::vector<Element*> pipe_elements;
        double                pipe_element_start_x;

        for (auto& r_element : r_current_model_part.Elements()) {
            if (r_element.GetProperties().Has(PIPE_START_ELEMENT)) {
                pipe_elements.push_back(&r_element);

                long unsigned int start_element = r_element.GetProperties()[PIPE_START_ELEMENT];

                KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                    << r_element.Id() << " " << start_element << std::endl;

                if (r_element.Id() == start_element) {
                    pipe_element_start_x = r_element.GetGeometry().GetPoint(0)[0];
                }
            }
        }

        if (pipe_elements.empty()) {
            return pipe_elements;
        }

        // Get Maximum X Value in Pipe
        auto right_pipe = std::max_element(pipe_elements.begin(), pipe_elements.end(),
                                           [](const Element* a, const Element* b) {
            return a->GetGeometry().GetPoint(0)[0] < b->GetGeometry().GetPoint(0)[0];
        });

        // Get Minimum X Value in Pipe
        auto left_pipe = std::min_element(pipe_elements.begin(), pipe_elements.end(),
                                          [](const Element* a, const Element* b) {
            return a->GetGeometry().GetPoint(0)[0] < b->GetGeometry().GetPoint(0)[0];
        });

        const auto min_x = left_pipe[0]->GetGeometry().GetPoint(0)[0];
        const auto max_x = right_pipe[0]->GetGeometry().GetPoint(0)[0];

        KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
            << min_x << " " << max_x << " " << pipe_element_start_x << std::endl;

        if ((min_x != pipe_element_start_x) && (max_x != pipe_element_start_x)) {
            KRATOS_ERROR << "Unable to determine pipe direction (multiple directions possible) -  "
                            "Check PIPE_START_ELEMENT"
                         << std::endl;
        }

        if (min_x == pipe_element_start_x) {
            // Pipe Left -> Right
            std::sort(pipe_elements.begin(), pipe_elements.end(), [](const Element* lhs, const Element* rhs) {
                return lhs->GetGeometry().GetPoint(0)[0] < rhs->GetGeometry().GetPoint(0)[0];
            });
        } else {
            // Pipe Right -> Left
            std::sort(pipe_elements.begin(), pipe_elements.end(), [](const Element* lhs, const Element* rhs) {
                return lhs->GetGeometry().GetPoint(0)[0] > rhs->GetGeometry().GetPoint(0)[0];
            });
        }

        KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
            << "Number of Pipe Elements: " << pipe_elements.size() << std::endl;
        for (const Element* pipeElement : pipe_elements) {
            KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                << "PipeElementIDs (in order): " << pipeElement->Id() << std::endl;
        }

        return pipe_elements;
    }

private:
    unsigned int mPipingIterations; /// This is used to calculate the pipingLength
    int          rank;
    double       mSmallPipeHeight    = 1e-10;
    double       mPipeHeightAccuracy = mSmallPipeHeight * 10;

    /// <summary>
    /// Initialises the number of open pipe elements. This value can be greater than 0 in a multi
    /// staged analysis.
    /// </summary>
    /// <param name="rPipeElements"></param>
    /// <returns></returns>
    template <typename FilteredElementsType>
    SizeType GetNumberOfActivePipeElements(const FilteredElementsType& rPipeElements)
    {
        auto is_pipe_active = [](auto p_element) { return p_element->GetValue(PIPE_ACTIVE); };
        return std::count_if(std::begin(rPipeElements), std::end(rPipeElements), is_pipe_active);
    }

    /// <summary>
    /// Calculates the maximum pipe height which the algorithm will allow
    /// </summary>
    /// <param name="rPipeElements"> vector of all pipe elements</param>
    /// <returns></returns>
    template <typename FilteredElementsType>
    double CalculateMaxPipeHeight(const FilteredElementsType& rPipeElements)
    {
        // get maximum pipe particle diameter over all pipe elements
        double max_diameter = 0;
        for (auto pipe_element : rPipeElements) {
            max_diameter = std::max(max_diameter, GeoTransportEquationUtilities::CalculateParticleDiameter(
                                                      pipe_element->GetProperties()));
        }

        // max pipe height is maximum pipe particle diameter * a constant
        const double height_factor = 100.;
        return max_diameter * height_factor;
    }

    /// <summary>
    /// Calculates the pipe height increment.
    /// </summary>
    /// <param name="MaxPipeHeight"> maximum allowed pipe height</param>
    /// <param name="NumberOfSteps"> number of non lineair piping iteration steps</param>
    /// <returns></returns>
    [[nodiscard]] double CalculatePipeHeightIncrement(double MaxPipeHeight, SizeType NumberOfSteps) const
    {
        return MaxPipeHeight / (static_cast<double>(NumberOfSteps) - 1.);
    }

    virtual bool Recalculate()
    {
        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0 && rank == 0)
            << "Recalculating" << std::endl;

        GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::InitializeSolutionStep();
        GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Predict();
        return GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::SolveSolutionStep();
    }

    virtual void BaseClassFinalizeSolutionStep()
    {
        // to override in a unit test
        GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::FinalizeSolutionStep();
    }

    template <typename FilteredElementsType>
    bool CheckPipeEquilibrium(const FilteredElementsType& rOpenPipeElements, double MaxPipeHeight, unsigned int MaxNumberOfPipingIterations)
    {
        bool         equilibrium      = false;
        bool         converged        = true;
        unsigned int piping_iteration = 0;

        // calculate max pipe height and pipe increment
        const auto pipe_height_increment =
            CalculatePipeHeightIncrement(MaxPipeHeight, MaxNumberOfPipingIterations);

        while (piping_iteration < MaxNumberOfPipingIterations && !equilibrium && converged) {
            equilibrium = true;

            // perform a flow calculation and stop growing if the calculation doesn't converge
            converged = this->Recalculate();

            if (!converged) continue;

            // Update depth of open piping Elements
            equilibrium = true;
            for (auto p_open_pipe_element : rOpenPipeElements) {
                // get open pipe element geometry and properties
                auto& r_geom = p_open_pipe_element->GetGeometry();
                auto& r_prop = p_open_pipe_element->GetProperties();

                // calculate equilibrium pipe height and get current pipe height
                double eq_height = p_open_pipe_element->CalculateEquilibriumPipeHeight(
                    r_prop, r_geom, p_open_pipe_element->GetValue(PIPE_ELEMENT_LENGTH));
                const auto current_height = p_open_pipe_element->GetValue(PIPE_HEIGHT);

                // set erosion on true if current pipe height is greater than the equilibrium height
                p_open_pipe_element->SetValue(PIPE_EROSION, p_open_pipe_element->GetValue(PIPE_EROSION) ||
                                                                current_height > eq_height);
                // check this if statement, I don't understand the check for pipe erosion
                if ((!p_open_pipe_element->GetValue(PIPE_EROSION) || current_height > eq_height) &&
                    current_height < MaxPipeHeight) {
                    p_open_pipe_element->SetValue(PIPE_HEIGHT, current_height + pipe_height_increment);
                    equilibrium = false;
                }

                // check if equilibrium height and current pipe heights are diverging, stop
                // Picard iterations if this is the case and set pipe height on zero
                if (!p_open_pipe_element->GetValue(PIPE_EROSION) && piping_iteration > 1 &&
                    eq_height - current_height > p_open_pipe_element->GetValue(DIFF_PIPE_HEIGHT)) {
                    p_open_pipe_element->SetValue(PIPE_HEIGHT, mSmallPipeHeight);
                    equilibrium = true;
                }
                // calculate difference between equilibrium height and current pipe height
                p_open_pipe_element->SetValue(DIFF_PIPE_HEIGHT, eq_height - current_height);
            }
            ++piping_iteration;
        }

        return equilibrium;
    }

    /// <summary>
    /// Checks the status of the tip element and checks if it should continue or stop growing
    /// </summary>
    /// <param name="NumberOfOpenPipeElements"> number of open pipe elements</param>
    /// <param name="NumberOfPipeELements"> number of pipe elements</param>
    /// <param name="MaxPipeHeight"> maximum allowed pipe height</param>
    /// <param name="PipeElements"> vector of all pipe elements</param>
    /// <returns>tuple of grow bool and number of open pipe elements</returns>
    template <typename FilteredElementsType>
    std::tuple<bool, SizeType> CheckStatusTipElement(SizeType NumberOfOpenPipeElements,
                                                     SizeType NumberOfPipeELements,
                                                     double   MaxPipeHeight,
                                                     const FilteredElementsType& rPipeElements)
    {
        bool grow = true;
        // check status of tip element, stop growing if pipe_height is zero or greater than maximum
        // pipe height or if all elements are open
        if (NumberOfOpenPipeElements < NumberOfPipeELements) {
            auto       p_tip_element = rPipeElements.at(NumberOfOpenPipeElements - 1);
            const auto pipe_height   = p_tip_element->GetValue(PIPE_HEIGHT);

            if (pipe_height > MaxPipeHeight + std::numeric_limits<double>::epsilon() ||
                pipe_height < mPipeHeightAccuracy) {
                // stable element found; pipe length does not increase during current time step
                grow = false;
                p_tip_element->SetValue(PIPE_EROSION, false);
                p_tip_element->SetValue(PIPE_ACTIVE, false);
                --NumberOfOpenPipeElements;

                KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                    << "Number of Open Pipe Elements: " << NumberOfOpenPipeElements << std::endl;
            }
        } else {
            grow = false;
        }
        return std::make_tuple(grow, NumberOfOpenPipeElements);
    }

    /// <summary>
    /// Saves pipe heights if pipe grows, else reset pipe heights to previous grow step.
    /// </summary>
    /// <param name="rOpenPipeElements"> open pipe elements</param>
    /// <param name="Grow"> boolean to check if pipe grows</param>
    /// <returns></returns>
    template <typename FilteredElementsType>
    void SaveOrResetPipeHeights(const FilteredElementsType& rOpenPipeElements, bool Grow)
    {
        for (auto p_element : rOpenPipeElements) {
            if (Grow) {
                p_element->SetValue(PREV_PIPE_HEIGHT, p_element->GetValue(PIPE_HEIGHT));
            } else {
                p_element->SetValue(PIPE_HEIGHT, p_element->GetValue(PREV_PIPE_HEIGHT));
            }
        }
    }

    template <typename FilteredElementsType>
    void DetermineOpenPipingElements(const FilteredElementsType& rPipingElements)
    {
        bool       grow                      = true;
        const auto number_of_piping_elements = rPipingElements.size();
        auto number_of_open_piping_elements  = this->GetNumberOfActivePipeElements(rPipingElements);
        auto max_pipe_height                 = CalculateMaxPipeHeight(rPipingElements);
        while (grow && (number_of_open_piping_elements < number_of_piping_elements)) {
            // get tip element and activate
            auto p_tip_element = rPipingElements.at(number_of_open_piping_elements);
            ++number_of_open_piping_elements;
            p_tip_element->SetValue(PIPE_ACTIVE, true);

            // Get all open pipe elements
            auto filter = [](auto p_element) {
                return p_element->Has(PIPE_ACTIVE) && p_element->GetValue(PIPE_ACTIVE);
            };
            auto OpenPipeElements = rPipingElements | boost::adaptors::filtered(filter);
            KRATOS_INFO_IF("PipingLoop", this->GetEchoLevel() > 0 && rank == 0)
                << "Number of Open Pipe Elements: " << boost::size(OpenPipeElements) << std::endl;

            // nonlinear Picard iteration, for deepening the pipe
            // Todo JDN (20220817):: Deal with Equilibrium redundancy
            CheckPipeEquilibrium(OpenPipeElements, max_pipe_height, mPipingIterations);

            // check if pipe should grow in length
            std::tie(grow, number_of_open_piping_elements) = CheckStatusTipElement(
                number_of_open_piping_elements, number_of_piping_elements, max_pipe_height, rPipingElements);

            // if n open elements is lower than total pipe elements, save pipe height current
            // growing iteration or reset to previous iteration in case the pipe should not grow.
            if (number_of_open_piping_elements < number_of_piping_elements) {
                SaveOrResetPipeHeights(OpenPipeElements, grow);
            }
            // recalculate groundwater flow
            const auto converged = this->Recalculate();

            KRATOS_ERROR_IF_NOT(converged) << "Groundwater flow calculation failed to converge." << std::endl;
        }
    }
}; // Class GeoMechanicsNewtonRaphsonErosionProcessStrategy
} // namespace Kratos
