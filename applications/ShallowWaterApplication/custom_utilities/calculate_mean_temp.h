#if !defined(KRATOS_CALCUALTE_MEAN_TEMP_UTILITY_INCLUDED )
#define  KRATOS_CALCUALTE_MEAN_TEMP_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"


namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class CalculateMeanTemperature 
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(CalculateMeanTemperature);

		CalculateMeanTemperature(ModelPart& model_part)
			: mr_model_part(model_part)              //mr_model_part is saved as private variable (declared at the end of the file)  
		{
			KRATOS_TRY	
			std::cout << "Hello, I am the constructor of the Utility" << std::endl; 
			KRATOS_CATCH("")	
		}
		

		~CalculateMeanTemperature()
		{}

		
 		void Calculate()  
		{
			KRATOS_TRY
			
			double area;                    //we create the needed variables
			double sum_areas=0.0;
			double sum_temperatures=0.0;
			double one_third=1.0/3.0;

			//getting data for the given geometry
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); //looping the elements
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				Geometry<Node<3> >& geom = ielem->GetGeometry(); 
				area=CalculateArea(geom);                            //we call CalculateArea (private function)  
				sum_areas+=area;
				for (unsigned int k = 0; k < 3; k++)
				{
					sum_temperatures += geom[k].FastGetSolutionStepValue(TEMPERATURE)*one_third*area;
				}
			}
			const double mean_temperature = sum_temperatures / sum_areas;
			std::cout << "Finished, the mean temperature is" << mean_temperature << std::endl;   //we print the result  
			
			
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:
		double CalculateArea(Element::GeometryType& geom)
		{
			return 0.5 * ((geom[1].X() - geom[0].X())*(geom[2].Y() - geom[0].Y())- (geom[1].Y() - geom[0].Y())*(geom[2].X() - geom[0].X()));
		}
	
		ModelPart& mr_model_part;

	};

}  // namespace Kratos.

#endif // KRATOS_CALCUALTE_MEAN_TEMP_UTILITY_INCLUDED  defined

