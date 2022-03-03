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
    	
        int PipeIter = 0;
        int openPipeElements = 0;
        bool Equilibrium = false;


        

        double amax = 0.02; //todo calculate this value
        double da = 0.0001; // todo calculate this value

    	auto PipeElements = GetPipingElements();

        unsigned int n_el = PipeElements.size();
        //PipeElements.size()
        std::vector<double> prev_pipe_heights;


        for (unsigned int i = 0; i < PipeElements.size(); ++i)
        {
            prev_pipe_heights.push_back(0);
        }
        // Open tip element of pipe (activate next pipe element)
        PipeElements.at(0).Set(ACTIVE, true);


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
                SteadyStatePwPipingElement& pElement = static_cast<SteadyStatePwPipingElement>(OpenPipeElement);
                if (!pElement.InEquilibrium(OpenPipeElement.GetProperties(), OpenPipeElement.GetGeometry()))
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
                    if (((!OpenPipeElement.GetValue(PIPE_EROSION) ||(current_height > eq_height)) && current_height < amax))
                    {
                        OpenPipeElement.SetValue(PIPE_HEIGHT, OpenPipeElement.GetValue(PIPE_HEIGHT)+da);
                        Equilibrium = false;
                    }

                    //// check divergence
                    //if (!activeSecondary[pp] && (k > 1) && ((cc - aa) > d[pp]))
                    //{
                    //    equilibrium = true;
                    //    *converged = true;
                    //    a[pp] = primaryErosionSelected ? b[pp] : -1.0;
                    //}


                    double test = OpenPipeElement.GetValue(PIPE_HEIGHT);
                    /*prop[MINIMUM_JOINT_WIDTH] = prop[MINIMUM_JOINT_WIDTH + 0.2;*/
                    double a = 1 + 1;
                       
                }
            	Recalculate();
                PipeIter += 1;  
            }
            else
            {
                // open new piping element
                openPipeElements += 1;
            	PipeElements.at(openPipeElements).Set(ACTIVE, true);
            }
    		
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
        GeoMechanicsNewtonRaphsonStrategy::InitializeSolutionStep();
        GeoMechanicsNewtonRaphsonStrategy::SolveSolutionStep();
    }
	
    //-----------------------------Get Piping Elements--------------------------------------

    std::vector<Kratos::Element> GetPipingElements() {
        ModelPart& CurrentModelPart = this->GetModelPart();
        std::vector<Element> PipeElements;
        for (const Element element: CurrentModelPart.Elements())
        {
        	if (element.GetProperties().Has(PIPE_D_70))
            {
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
