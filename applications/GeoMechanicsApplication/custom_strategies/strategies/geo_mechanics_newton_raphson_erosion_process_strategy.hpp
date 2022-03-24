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

#if !defined(KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_EROSION_PROCESS_STRATEGY)
#define KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_EROSION_PROCESS_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "geo_mechanics_application_variables.h"
#include "custom_elements/steady_state_Pw_piping_element.hpp"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "boost/range/adaptor/filtered.hpp"

#include <iostream>
#include <chrono>
#include <ctime> 

namespace Kratos
{

bool isOpen(Element* element) {
    if (element->Has(PIPE_ACTIVE))
    {
        return element->GetValue(PIPE_ACTIVE);
    }
    else
        return true;
}
	
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicsNewtonRaphsonErosionProcessStrategy :
    public GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsNewtonRaphsonErosionProcessStrategy);

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>              BaseType;
    typedef GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>  MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace>                 TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType                          TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType                                              TSchemeType;
    typedef typename BaseType::DofsArrayType                                          DofsArrayType;
    typedef typename BaseType::TSystemMatrixType                                  TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType                                  TSystemVectorType;


    typedef typename SteadyStatePwPipingElement<2, 4>                    SteadyStatePwPipingElement;
    typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
	
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ///Constructor
    GeoMechanicsNewtonRaphsonErosionProcessStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
		) : GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
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

    	mPipingIterations = rParameters["max_piping_iterations"].GetInt();
    	
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Destructor
    ~GeoMechanicsNewtonRaphsonErosionProcessStrategy() override {}
	

    void FinalizeSolutionStep() override
    {
        KRATOS_INFO("PipingLoop") << "Max Piping Iterations: " << mPipingIterations << std::endl;

        int openPipeElements = 0;
        bool grow = true;

        // get piping elements
        std::vector<Element*> PipeElements = GetPipingElements();
        unsigned int n_el = PipeElements.size(); // number of piping elements

        // calculate max pipe height and pipe increment
        double amax = CalculateMaxPipeHeight(PipeElements);
        double da = CalculatePipeHeightIncrement(amax, mPipingIterations);

        
        // continue this loop, while the pipe is growing in length
        while (grow)
        {
            bool Equilibrium = false;
            bool converged = true;
            int PipeIter = 0;

            // get tip element and open with 1 pipe height increment
            Element* tip_element = PipeElements.at(openPipeElements);
            openPipeElements += 1;
            tip_element->SetValue(PIPE_ACTIVE, true);
            //tip_element->SetValue(PIPE_HEIGHT, da);

            // non-lin picard iteration, for deepening the pipe
            while (PipeIter < mPipingIterations && !Equilibrium && converged)
            {

                // set equilibirum on true
                Equilibrium = true;

                // Get all open pipe elements
                auto OpenPipeElements = PipeElements | boost::adaptors::filtered(isOpen);
                KRATOS_INFO("PipingLoop") << "Number of Open Pipe Elements: " << boost::size(OpenPipeElements) << std::endl;

                // todo check if this part can be removed/refactored
                // Check if elements are in equilibrium
                for (auto OpenPipeElement : OpenPipeElements)
                {
                    SteadyStatePwPipingElement* pElement = static_cast<SteadyStatePwPipingElement*>(OpenPipeElement);
                    if (!pElement->InEquilibrium(OpenPipeElement->GetProperties(), OpenPipeElement->GetGeometry()))
                    {
                        Equilibrium = false;
                        break;
                    }
                }

                // perform a flow calculation and stop growing if the calculation doesnt converge

                // get underlying buffer
                //std::streambuf* orig_buf = std::cout.rdbuf();

                // set null
               // auto start = std::chrono::system_clock::now();
                //std::cout.rdbuf(NULL);
                //int t = 0;
               // while (t < 100)
                //{
                    converged = Recalculate();
                  //  t += 1;
                //}


                // restore buffer
               // std::cout.rdbuf(orig_buf);
                //auto end = std::chrono::system_clock::now();
                //std::chrono::duration<double> elapsed_seconds = end - start;
                //std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

                //converged = Recalculate();
                if (!converged)
                {
                    grow = false;
                }

                if (converged)
                {
                    // Update depth of open piping Elements 
                    Equilibrium = true;
                    for (auto OpenPipeElement : OpenPipeElements)
                    {
                        auto pElement = static_cast<SteadyStatePwPipingElement*>(OpenPipeElement);

                        // get open pipe element geometry and properties
                        auto& Geom = OpenPipeElement->GetGeometry();
                        auto& prop = OpenPipeElement->GetProperties();

                        // todo set this property in the input
                        prop.SetValue(PIPE_MODEL_FACTOR, 1);

                        // calculate equilibrium pipe height and get current pipe height
                        double eq_height = pElement->CalculateEquilibriumPipeHeight(prop, Geom, OpenPipeElement->GetValue(PIPE_ELEMENT_LENGTH));
                        double current_height = OpenPipeElement->GetValue(PIPE_HEIGHT);

                        // set erosion on true if current pipe height is greater than the equilibirum height
                        if (current_height > eq_height)
                        {
                            OpenPipeElement->SetValue(PIPE_EROSION, true);
                        }

                        
                        // check this if statement, I dont understand the check for pipe erosion
                        if (((!OpenPipeElement->GetValue(PIPE_EROSION) || (current_height > eq_height)) && current_height < amax))
                        {
                            OpenPipeElement->SetValue(PIPE_HEIGHT, OpenPipeElement->GetValue(PIPE_HEIGHT) + da);
                            Equilibrium = false;
                        }

                        // check if equilibrium height and current pipe heights are diverging, stop picard iterations if this is the case and set pipe height on zero
                        if (!OpenPipeElement->GetValue(PIPE_EROSION) && (PipeIter > 1) && ((eq_height - current_height) > OpenPipeElement->GetValue(DIFF_PIPE_HEIGHT)))
                        {
                            Equilibrium = true;
                            OpenPipeElement->SetValue(PIPE_HEIGHT, 1e-10);
                        }

                        // calculate difference between equilibrium height and current pipe height
                        OpenPipeElement->SetValue(DIFF_PIPE_HEIGHT, eq_height - current_height);

                    }
                    // increment piping iteration number
                    PipeIter += 1;
                }
            }

            // get open pipe elements
            auto OpenPipeElements = PipeElements | boost::adaptors::filtered(isOpen);

            // check status of tip element, stop growing if pipe_height is zero or greater than maximum pipe height or if all elements are open
            if (openPipeElements < n_el)
            {
                tip_element = PipeElements.at(openPipeElements-1);
                double pipe_height = tip_element->GetValue(PIPE_HEIGHT);
                if ((pipe_height > amax + 1e-10) || (pipe_height < 1e-9))
                {
                    // stable element found; pipe length does not increase during current time step
                    grow = false;
                    tip_element->SetValue(PIPE_EROSION, false);
                    tip_element->SetValue(PIPE_ACTIVE, false);
                    openPipeElements -= 1;
                    auto OpenPipeElements = PipeElements | boost::adaptors::filtered(isOpen);

                }
            }
            else
            {
                grow = false;
            }

            // save pipe height current growing iteration or reset to previous iteration in case the pipe should not grow.
            if (openPipeElements < n_el)
            {
                for (auto OpenPipeElement : OpenPipeElements)
                {
                    if (grow)
                    {
                        OpenPipeElement->SetValue(PREV_PIPE_HEIGHT, OpenPipeElement->GetValue(PIPE_HEIGHT));
                    }
                    else
                    {
                        OpenPipeElement->SetValue(PIPE_HEIGHT, OpenPipeElement->GetValue(PREV_PIPE_HEIGHT));
                    }
                }
            }

            // get underlying buffer
            std::streambuf* orig_buf = std::cout.rdbuf();

            // set null
            std::cout.rdbuf(NULL);
            converged = Recalculate();

            // restore buffer
            std::cout.rdbuf(orig_buf);
        }          

        GeoMechanicsNewtonRaphsonStrategy::FinalizeSolutionStep();
        
	}

	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    unsigned int mPipingIterations; /// This is used to calculate the pipingLength

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    double CalculateParticleDiameter(const PropertiesType& Prop)
    {
        double diameter;

        if (Prop[PIPE_MODIFIED_D])
            diameter = 2.08e-4 * pow((Prop[PIPE_D_70] / 2.08e-4), 0.4);
        else
            diameter = Prop[PIPE_D_70];
        return diameter;
    }

    double CalculateMaxPipeHeight(std::vector<Element*> pipe_elements)
    {
        double max_diameter = 0;
        double particle_diameter = 0;
        double height_factor = 100;

        for (Element* pipe_element : pipe_elements)
        {
            PropertiesType prop = pipe_element->GetProperties();
            particle_diameter = this->CalculateParticleDiameter(prop);
            if (particle_diameter > max_diameter)
            {
                max_diameter = particle_diameter;
            }
        }
        return max_diameter * height_factor;
    }

    double CalculatePipeHeightIncrement(double max_pipe_height, const unsigned int n_steps)
    {
        return max_pipe_height / (n_steps - 1);
    }


    bool Recalculate()
    {

        KRATOS_INFO("PipingLoop") << "Recalculating" << std::endl;
    	//ModelPart& CurrentModelPart = this->GetModelPart();
        //this->Clear();

        // Reset displacements to the initial (Assumes Water Pressure is the convergence criteria)
       /* block_for_each(CurrentModelPart.Nodes(), [&](Node<3>& rNode) {
            auto dold = rNode.GetSolutionStepValue(WATER_PRESSURE, 1);
            rNode.GetSolutionStepValue(WATER_PRESSURE, 0) = dold;
            });*/

    	GeoMechanicsNewtonRaphsonStrategy::InitializeSolutionStep();
        GeoMechanicsNewtonRaphsonStrategy::Predict();
        return GeoMechanicsNewtonRaphsonStrategy::SolveSolutionStep();

    }
	
    //-----------------------------Get Piping Elements--------------------------------------

    std::vector<Element*> GetPipingElements() {
        ModelPart& CurrentModelPart = this->GetModelPart();
        std::vector<Element*> PipeElements;
        double PipeElementStartX;
        ModelPart::ElementsContainerType::iterator element_it = CurrentModelPart.ElementsBegin();

        for (Element& element: CurrentModelPart.Elements())
        {
        	if (element.GetProperties().Has(PIPE_START_ELEMENT))
            {
               // ModelPart::ElementsContainerType::iterator element_it = rModelPart.ElementsBegin()
                element.SetValue(PIPE_ACTIVE, false);
                element.SetValue(PREV_PIPE_HEIGHT, 1e-10);
                element.SetValue(DIFF_PIPE_HEIGHT, 0);

                PipeElements.push_back(&element);
                
                auto startElement = element.GetProperties()[PIPE_START_ELEMENT];

        		KRATOS_INFO("PipingLoop") << element.Id() << " " << startElement << std::endl;

        		if (element.Id() == startElement)
                {
                    PipeElementStartX = element.GetGeometry().GetPoint(0)[0];
                }
            }
        }

        if (PipeElements.size() == 0)
        {
            return PipeElements;
        }

        // Get Maximum X Value in Pipe
        auto rightPipe = std::max_element(PipeElements.begin(), PipeElements.end(), [](const Element* a, const Element* b)
            {
                return a->GetGeometry().GetPoint(0)[0] < b->GetGeometry().GetPoint(0)[0];
            });

        // Get Minimum X Value in Pipe
        auto leftPipe = std::min_element(PipeElements.begin(), PipeElements.end(), [](const Element* a, const Element* b)
            {
                return a->GetGeometry().GetPoint(0)[0] < b->GetGeometry().GetPoint(0)[0];
            });

        double minX = leftPipe[0]->GetGeometry().GetPoint(0)[0];
        double maxX = rightPipe[0]->GetGeometry().GetPoint(0)[0];

        KRATOS_INFO("PipingLoop") << minX << " " << maxX << " " << PipeElementStartX << std::endl;

    	if ((minX != PipeElementStartX) && (maxX != PipeElementStartX))
        {
            KRATOS_ERROR << "Unable to determine pipe direction (multiple directions possible) -  Check PIPE_START_ELEMENT" << std::endl;
        }

        if (minX == PipeElementStartX)
        {
            // Pipe Left -> Right
            sort(PipeElements.begin(), PipeElements.end(), [](const Element* lhs, const Element* rhs) {
                return lhs->GetGeometry().GetPoint(0)[0] < rhs->GetGeometry().GetPoint(0)[0];
                });
        }
        else
        {
            // Pipe Right -> Left
            sort(PipeElements.begin(), PipeElements.end(), [](const Element* lhs, const Element* rhs) {
                return lhs->GetGeometry().GetPoint(0)[0] > rhs->GetGeometry().GetPoint(0)[0];
                });
        }

        KRATOS_INFO("PipingLoop") << "Number of Pipe Elements: " << PipeElements.size() << std::endl;
        for (const Element* pipeElement: PipeElements)
        {
            KRATOS_INFO("PipingLoop") << "PipeElementIDs (in order): " << pipeElement->Id() << std::endl;
        }

    	return PipeElements;
    }
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class GeoMechanicsNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_EROSION_PROCESS_STRATEGY  defined
