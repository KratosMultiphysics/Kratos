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

#include "boost/range/adaptor/filtered.hpp"

namespace Kratos
{

bool isOpen(Element element) {
    if (element.IsDefined(ACTIVE))
    {
        return element.Is(ACTIVE);
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

        double amax = 0.02; //todo calculate this value
        double da = 0.00001; // todo calculate this value

    	auto PipeElements = GetPipingElements();

        unsigned int n_el = PipeElements.size();
        //PipeElements.size()
        std::vector<double> prev_pipe_heights;
        std::vector<double> diff_pipe_heights;

        // initialise pipe heights in previous timestep as 0
        // todo, check how to handle multiple construction stages
        for (unsigned int i = 0; i < PipeElements.size(); ++i)
        {
            prev_pipe_heights.push_back(0);
            diff_pipe_heights.push_back(0);
        }


        bool grow = true;
        while (grow)
        {
            bool Equilibrium = false;
            int PipeIter = 0;

            // Open tip element of pipe (activate next pipe element)
            PipeElements.at(openPipeElements).Set(ACTIVE, true);
            openPipeElements += 1;

            // Loop over open pipe elements
            auto OpenPipeElements = PipeElements | boost::adaptors::filtered(isOpen);
            KRATOS_INFO("PipingLoop") << "Number of Open Pipe Elements: " << boost::size(OpenPipeElements) << std::endl;

            // Implement Piping Loop (non-lin picard iteration)
            while (PipeIter < mPipingIterations && !Equilibrium)
            {

                // Update the pipe_height by the pipe increment

                // Sellmeijer Piping Method 
                Equilibrium = true;


                Recalculate();

                // todo refactor
                // Check if elements are in equilibrium
                //for (auto OpenPipeElement : OpenPipeElements)
          //      {
          //          SteadyStatePwPipingElement& pElement = static_cast<SteadyStatePwPipingElement>(OpenPipeElement);
          //          if (!pElement.InEquilibrium(OpenPipeElement.GetProperties(), OpenPipeElement.GetGeometry()))
          //          {
          //              Equilibrium = false;
          //              break;
          //          }          
          //      }

           /*     if (!Equilibrium)
                {*/
                unsigned int pp = 0;
                // Update Piping Elements 
                for (auto OpenPipeElement : OpenPipeElements)
                {
                    auto& pElement = static_cast<SteadyStatePwPipingElement>(OpenPipeElement);

                    // Update here
                    auto& Geom = OpenPipeElement.GetGeometry();
                    auto& prop = OpenPipeElement.GetProperties();

                    // todo set this property in the input
                    prop.SetValue(PIPE_MODEL_FACTOR, 1);

                    // calculate eq height
                    double eq_height = pElement.CalculateEquilibriumPipeHeight(prop, Geom, OpenPipeElement.GetValue(PIPE_ELEMENT_LENGTH));
                    double current_height = OpenPipeElement.GetValue(PIPE_HEIGHT);

                    if (current_height > eq_height)
                    {
                        OpenPipeElement.SetValue(PIPE_EROSION, true);
                    }

                    // todo check max pipe height

                    // check this if statement, I dont understand the check for pipe erosion
                    if (((!OpenPipeElement.GetValue(PIPE_EROSION) || (current_height > eq_height)) && current_height < amax))
                    {
                        OpenPipeElement.SetValue(PIPE_HEIGHT, OpenPipeElement.GetValue(PIPE_HEIGHT) + da);
                        Equilibrium = false;
                    }

                    //check max pipe height
                    if ((current_height > eq_height) && (boost::size(OpenPipeElements) < PipeElements.size()) && (OpenPipeElement.GetValue(PIPE_HEIGHT) > amax))
                    {
                        KRATOS_ERROR << "max pipe height is too low" << std::endl;
                    }

                    // check divergence
                    if (!OpenPipeElement.GetValue(PIPE_EROSION) && (PipeIter > 1) && ((eq_height - current_height) > diff_pipe_heights[pp]))
                    {
                        Equilibrium = true;
                        //*converged = true;
                        OpenPipeElement.GetValue(PIPE_HEIGHT) = -1.0;
                    }

                    diff_pipe_heights[pp] = eq_height - current_height;

                    double test = OpenPipeElement.GetValue(PIPE_HEIGHT);
                    pp = pp + 1;

                }
                PipeIter += 1;
                //}
                //else
                //{
                //    // open new piping element
                //    openPipeElements += 1;
                //	PipeElements.at(openPipeElements).Set(ACTIVE, true);
                //}

            }

            // check status tip element
            bool grow = true;
            if (boost::size(OpenPipeElements) < PipeElements.size())
            {
                unsigned int tip = boost::size(OpenPipeElements) - 1;
                auto tip_element = PipeElements.at(tip);

                // no growth, close tip element
                if ((tip_element.GetValue(PIPE_HEIGHT) < DBL_EPSILON) || (tip_element.GetValue(PIPE_HEIGHT) > amax + DBL_EPSILON))
                {
                    grow = false;
                    tip_element.SetValue(PIPE_EROSION, false);
                    //tip_element.Set(ACTIVE, false);
                    openPipeElements -= 1;

                    // get open elements
                    //OpenPipeElements = PipeElements | boost::adaptors::filtered(isOpen);
                }
                else
                {
                    std::cout << "n_open: " << tip << "; height tip: " << tip_element.GetValue(PIPE_HEIGHT);
                }
            }
            if (boost::size(OpenPipeElements) < PipeElements.size())
            {
                unsigned int pp = 0;
                for (auto OpenPipeElement : OpenPipeElements)
                {
                    if (grow)
                    {
                        prev_pipe_heights[pp] = OpenPipeElement.GetValue(PIPE_HEIGHT);
                    }
                    else
                    {
                        OpenPipeElement.GetValue(PIPE_HEIGHT) = prev_pipe_heights[pp];
                    }
                    pp = pp + 1;
                }
            }

            Recalculate();
        }
        GeoMechanicsNewtonRaphsonStrategy::FinalizeSolutionStep();
        
	}

   

	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    unsigned int mPipingIterations; /// This is used to calculate the pipingLength

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    void Recalculate()
    {
        auto nodes = GeoMechanicsNewtonRaphsonStrategy::GetModelPart().Nodes();
        //for (auto node : GeoMechanicsNewtonRaphsonStrategy::GetModelPart().Nodes())
        //{
        //    double b = 1 + 1;
        //    //dold = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)
        //    //    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, dold)
        //}
        for (unsigned int i = 0; i < nodes.size(); ++i)
        {
            auto dold = nodes[i].GetSolutionStepValue(WATER_PRESSURE, 1);
            nodes[i].GetSolutionStepValue(WATER_PRESSURE, 0) = dold;
        }
        GeoMechanicsNewtonRaphsonStrategy::GetModelPart().GetProcessInfo()[TIME] = 0;
        GeoMechanicsNewtonRaphsonStrategy::GetModelPart().GetProcessInfo()[STEP] = 0;
        double time = GeoMechanicsNewtonRaphsonStrategy::GetModelPart().GetProcessInfo()[TIME];
        GeoMechanicsNewtonRaphsonStrategy::InitializeSolutionStep();
        GeoMechanicsNewtonRaphsonStrategy::SolveSolutionStep();
        //GeoMechanicsNewtonRaphsonStrategy::FinalizeSolutionStep();
        time = GeoMechanicsNewtonRaphsonStrategy::GetModelPart().GetProcessInfo()[TIME];
        double d_t = GeoMechanicsNewtonRaphsonStrategy::GetModelPart().GetProcessInfo()[DELTA_TIME];
        double a = 1 + 1;
        //self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
    }
	
    //-----------------------------Get Piping Elements--------------------------------------

    std::vector<Kratos::Element> GetPipingElements() {
        ModelPart& CurrentModelPart = this->GetModelPart();
        std::vector<Element> PipeElements;
        for (Element& element: CurrentModelPart.Elements())
        {
        	if (element.GetProperties().Has(PIPE_D_70))
            {
                element.Set(ACTIVE, false);
                PipeElements.push_back(element);
                
                
            }
        }

        // We assume that the piping elements are provided in order.
    
    	// todo JDN  - Make sure elements are sorted correctly

        KRATOS_INFO("PipingLoop") << "Number of Pipe Elements: " << PipeElements.size() << std::endl;
        return PipeElements;
    }
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class GeoMechanicsNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_EROSION_PROCESS_STRATEGY  defined
