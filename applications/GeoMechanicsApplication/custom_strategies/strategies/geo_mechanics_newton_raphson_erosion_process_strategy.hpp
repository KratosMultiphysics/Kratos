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
        bool Equilibrium = false;

    	auto PipeElements = GetPipingElements();

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
			
    		for (auto OpenPipeElement : OpenPipeElements)
            {
				// For all open pipe elements:
                //	if pipe_element_height > critical_height_element:
                //      Equilibrium = false // Loop again
                //      update the pipe_height by the pipe increment

            	KRATOS_INFO("PipingLoop") << "Pipe Element: " << OpenPipeElement.Id() << " = " << OpenPipeElement.Is(ACTIVE) << std::endl;
    			

            }
    			// JN TODO: For transient analysis reset state to start of solution step ???
                Recalculate();
                PipeIter += 1;
                Equilibrium = false;
                PipeElements.at(PipeIter).Set(ACTIVE, true);
            
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

        KRATOS_INFO("PipingLoop") << "Number of Pipe Elements: " << PipeElements.size() << std::endl;
        return PipeElements;
    }
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class GeoMechanicsNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_EROSION_PROCESS_STRATEGY  defined
