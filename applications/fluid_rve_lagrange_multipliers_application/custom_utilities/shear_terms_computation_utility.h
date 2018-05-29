

#if !defined(KRATOS_SHEAR_TERM_COMPUTATION_UTILITY_INCLUDED )
#define KRATOS_SHEAR_TERM_COMPUTATION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "fluid_rve_lagrange_multipliers_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"


namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class ShearTermsComputationUtility2D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(ShearTermsComputationUtility2D);

		ShearTermsComputationUtility2D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 2d Utility" << std::endl;
			Check();		
			KRATOS_CATCH("")	
		}
		

		~ShearTermsComputationUtility2D()
		{}

		
		void ComputeTau(array_1d<double, 3 > & first_row,array_1d<double, 3 > & second_row)
		{
			KRATOS_TRY

			const int TDim = 2;

			boost::numeric::ublas::bounded_matrix<double, TDim, TDim > grad_U_UxU_total =  ZeroMatrix(TDim, TDim);


			boost::numeric::ublas::bounded_matrix<double,3*(TDim-1), TDim+1 > Ngauss; //GP, N_value(1-3)
			Ngauss(0,0)=2.0/3.0; Ngauss(0,1)=1.0/6.0; Ngauss(0,2)=1.0/6.0;
			Ngauss(1,0)=1.0/6.0; Ngauss(1,1)=2.0/3.0; Ngauss(1,2)=1.0/6.0;
			Ngauss(2,0)=1.0/6.0; Ngauss(2,1)=1.0/6.0; Ngauss(2,2)=2.0/3.0;

			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); 
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				boost::numeric::ublas::bounded_matrix<double,TDim+1, TDim > coords;
				boost::numeric::ublas::bounded_matrix<double,TDim+1, TDim > velocity;

				Geometry<Node<3> >& geom = ielem->GetGeometry();
				double viscosity = 0.0;
				for (unsigned int i = 0; i < 3; i++) //3 nodes
				{
					const array_1d<double, 3 > & xyz = geom[i].Coordinates();
					const array_1d<double, 3 > & vel_node = geom[i].FastGetSolutionStepValue(VELOCITY);

					//volumes[i] = 0.0;
					//KRATOS_WATCH(distances(i));
					for (unsigned int j = 0; j < TDim; j++)
					{
						coords(i, j) = xyz[j];
						velocity(i, j) = vel_node[j];
					}
					viscosity +=geom[i].GetSolutionStepValue(VISCOSITY)*1.0/3.0;
				}
				
				double Area;
				boost::numeric::ublas::bounded_matrix<double, (TDim+1), TDim > DN_DX;
				array_1d<double, (TDim+1) > N_singleGP;
				GeometryUtils::CalculateGeometryData(geom, DN_DX, N_singleGP, Area);
				
				boost::numeric::ublas::bounded_matrix<double, TDim, TDim > grad_U = prod ( trans(velocity), DN_DX ); //constant for all GPs
				boost::numeric::ublas::bounded_matrix<double, TDim, TDim > viscous_part  = 2.0 * Area * viscosity  * 0.5 * ( grad_U + trans(grad_U) )  ;

				
				boost::numeric::ublas::bounded_matrix<double, TDim, TDim > grad_U_UxU =  ZeroMatrix(TDim, TDim);
				
				//the term we must compute is grad_U UxX . Since UxX is quadratic, we must use 3 GP:
				for (unsigned int i = 0; i < 3; i++)
				{
					boost::numeric::ublas::bounded_matrix<double, 1 , TDim > X_GP = ZeroMatrix(1,TDim);
					boost::numeric::ublas::bounded_matrix<double, TDim , 1 > U_GP = ZeroMatrix(TDim,1);
					for (unsigned int j = 0; j < 3; j++) //3 nodes
					{
						X_GP(0,0) += coords(j, 0) * Ngauss(i,j);
						X_GP(0,1) += coords(j, 1) * Ngauss(i,j);
						
						U_GP(0,0) += velocity(j, 0) * Ngauss(i,j);
						U_GP(1,0) += velocity(j, 1) * Ngauss(i,j);
					}
					boost::numeric::ublas::bounded_matrix<double, TDim , TDim > UxX_GP = prod(U_GP,X_GP);

					grad_U_UxU += Area*1.0/3.0*prod(grad_U,UxX_GP) ;
				}
				
				grad_U_UxU_total += grad_U_UxU + viscous_part ;
			}
			
			//computing the trace:
			boost::numeric::ublas::bounded_matrix<double, TDim, TDim > volumetric_part =  1.0/2.0 * (grad_U_UxU_total(0,0) + grad_U_UxU_total(1,1) ) * IdentityMatrix(TDim,TDim) ;
			
			//substracting
			grad_U_UxU_total -= volumetric_part;
			
			
			first_row[0]= grad_U_UxU_total(0,0);       first_row[1]= grad_U_UxU_total(0,1);
			second_row[0]= grad_U_UxU_total(1,0);       second_row[1]= grad_U_UxU_total(1,1);

			
			//std::cout << "Finished adding lagrange multipliers for the velocity gradients (ELEMENT BASED)" << condition_number << std::endl;
			


			
			KRATOS_CATCH("")
		} 
		
	protected:

		void Check()
		{
			//if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NORMAL) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data","");
			//if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data","");
			//if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS variable on solution step data","");
		}

	private:
		ModelPart& mr_model_part;

	};
	
	
	
	

}  // namespace Kratos.

#endif // KRATOS_MEAN_VELOCITY_GRADIENTS_LAGRANGE_MULTIPLIER_CONDITION_UTILITY_INCLUDED


