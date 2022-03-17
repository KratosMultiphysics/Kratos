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

bool isOpen(Element* element) {
    if (element->IsDefined(ACTIVE))
    {
        return element->Is(ACTIVE);
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
    	
        int PipeIter = 0;
        int openPipeElements = 0;
        bool Equilibrium = false;

        double amax = 0.02; //todo calculate this value
        double da = 0.00001; // todo calculate this value

        std::vector<Element*> PipeElements = GetPipingElements();

        unsigned int n_el = PipeElements.size();
        //PipeElements.size()
        std::vector<double> prev_pipe_heights;


        for (unsigned int i = 0; i < PipeElements.size(); ++i)
        {
            prev_pipe_heights.push_back(0);
        }
        // Open tip element of pipe (activate next pipe element)
        Element* test_element = PipeElements.at(0);
        test_element->Set(ACTIVE, true);


    	// Implement Piping Loop (non-lin picard iteration)
    	while (PipeIter < mPipingIterations && !Equilibrium)
        {
            
    		// Update the pipe_height by the pipe increment
    		
    		// Sellmeijer Piping Method 
            Equilibrium = true;
    		
            // Loop over open pipe elements
            auto OpenPipeElements = PipeElements | boost::adaptors::filtered(isOpen);
            KRATOS_INFO("PipingLoop") << "Number of Open Pipe Elements: " << boost::size(OpenPipeElements) << std::endl;
			
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



    		if (!Equilibrium)
            {
                // Update Piping Elements 
            	for (auto OpenPipeElement : OpenPipeElements)
                {
                    auto pElement = static_cast<SteadyStatePwPipingElement*>(OpenPipeElement);

                    // Update here
                    auto& Geom = OpenPipeElement->GetGeometry();
                    auto& prop = OpenPipeElement->GetProperties();

                    // todo set this property in the input
                    prop.SetValue(PIPE_MODEL_FACTOR, 1);

                    // calculate eq height
                    double eq_height = pElement->CalculateEquilibriumPipeHeight(prop, Geom, OpenPipeElement->GetValue(PIPE_ELEMENT_LENGTH));
                    double current_height = OpenPipeElement->GetValue(PIPE_HEIGHT);

                    if (current_height > eq_height)
                    {
                        OpenPipeElement->SetValue(PIPE_EROSION, true);
                    }

                    // todo check max pipe height

                    // check this if statement, I dont understand the check for pipe erosion
                    if (((!OpenPipeElement->GetValue(PIPE_EROSION) ||(current_height > eq_height)) && current_height < amax))
                    {
                        OpenPipeElement->SetValue(PIPE_HEIGHT, OpenPipeElement->GetValue(PIPE_HEIGHT)+da);
                        Equilibrium = false;
                    }

                    //// check divergence
                    //if (!activeSecondary[pp] && (k > 1) && ((cc - aa) > d[pp]))
                    //{
                    //    equilibrium = true;
                    //    *converged = true;
                    //    a[pp] = primaryErosionSelected ? b[pp] : -1.0;
                    //}


                    double test = OpenPipeElement->GetValue(PIPE_HEIGHT);
                    /*prop[MINIMUM_JOINT_WIDTH] = prop[MINIMUM_JOINT_WIDTH + 0.2;*/
                    double a = 1 + 1;
                       
                }
            	bool converged = Recalculate();
                PipeIter += 1;  
            }
            else
            {
                // open new piping element
                openPipeElements += 1;
            	PipeElements.at(openPipeElements)->Set(ACTIVE, true);
            }
    		
        }

        GeoMechanicsNewtonRaphsonStrategy::FinalizeSolutionStep();
        
	}

   

	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    unsigned int mPipingIterations; /// This is used to calculate the pipingLength

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    bool Recalculate()
    {

        KRATOS_INFO("PipingLoop") << "Recalculating" << std::endl;
    	ModelPart& CurrentModelPart = this->GetModelPart();
        this->Clear();

        // Reset displacements to the initial (Assumes Water Pressure is the convergence criteria)
        block_for_each(CurrentModelPart.Nodes(), [&](Node<3>& rNode) {
            auto dold = rNode.GetSolutionStepValue(WATER_PRESSURE, 1);
            rNode.GetSolutionStepValue(WATER_PRESSURE, 0) = dold;
            });

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
                element.Set(ACTIVE, false);
                PipeElements.push_back(&element);
                
                auto startElement = element.GetProperties()[PIPE_START_ELEMENT];

                // Debug Statement
        		// startElement = 217;

        		KRATOS_INFO("PipingLoop") << element.Id() << " " << startElement << std::endl;

        		if (element.Id() == startElement)
                {
                    PipeElementStartX = element.GetGeometry().GetPoint(0)[0];
                }
            }
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
