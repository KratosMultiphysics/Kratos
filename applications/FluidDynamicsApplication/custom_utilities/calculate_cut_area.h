#if !defined(KRATOS_CALCULATE_CUT_AREA_H_INCLUDED)
#define  KRATOS_CALCULATE_CUT_AREA_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"


namespace Kratos
{
	class CalculateCutArea 
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(CalculateCutArea);

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

		CalculateCutArea(ModelPart& model_part)
			: mr_model_part(model_part)              //mr_model_part is saved as private variable (declared at the end of the file)  
		{
			KRATOS_TRY	
			KRATOS_CATCH("")	
		}
		

		~CalculateCutArea()
		{}

		
 		double Calculate()  
		{
			KRATOS_TRY
			
			double neg_vol = 0.0;              
			#pragma omp parallel for reduction(+: neg_vol)

			//getting data for the given geometry
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); //looping the elements
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				Matrix shape_functions;
                GeometryType::ShapeFunctionsGradientsType shape_derivatives;
                Geometry<Node<3>> & rGeom = ielem->GetGeometry(); 
                unsigned int pt_count_neg = 0;

                // instead of using data.isCut()
                for (unsigned int pt = 0; pt < rGeom.Points().size(); pt++){
                    if ( rGeom[pt].FastGetSolutionStepValue(DISTANCE) <= 0.0 ){
                        pt_count_neg++;
                    }
                }

                if ( pt_count_neg == rGeom.PointsNumber() ){
                    // all nodes are negative (pointer is necessary to maintain polymorphism of DomainSize())
                    neg_vol += ielem->pGetGeometry()->DomainSize();
                }
                else if ( 0 < pt_count_neg ){
                    // element is cut by the surface (splitting)
                    Kratos::unique_ptr<ModifiedShapeFunctions> p_modified_sh_func = nullptr;
                    Vector w_gauss_neg_side(3, 0.0);

                    Vector Distance( rGeom.PointsNumber(), 0.0 );
                    for (unsigned int i = 0; i < rGeom.PointsNumber(); i++){
                        // Control mechanism to avoid 0.0 ( is necessary because "distance_modification" possibly not yet executed )
                        double& r_dist = rGeom[i].FastGetSolutionStepValue(DISTANCE);
                        if (std::abs(r_dist) < 1.0e-12) {
                            const double aux_dist = 1.0e-6* rGeom[i].GetValue(NODAL_H);
                            if (r_dist > 0.0) {
                                #pragma omp critical
                                r_dist = aux_dist;
                            } else {
                                #pragma omp critical
                                r_dist = -aux_dist;
                            }
                        }
                        Distance[i] = rGeom[i].FastGetSolutionStepValue(DISTANCE);
                    }


                    if ( rGeom.PointsNumber() == 3 ){ p_modified_sh_func = Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(ielem->pGetGeometry(), Distance); }
                    else if ( rGeom.PointsNumber() == 4 ){ p_modified_sh_func = Kratos::make_unique<Tetrahedra3D4ModifiedShapeFunctions>(ielem->pGetGeometry(), Distance); }
                    else { KRATOS_ERROR << "The process can not be applied on this kind of element" << std::endl; }

                    // Call the negative side modified shape functions calculator
                    // Object p_modified_sh_func has full knowledge of slit geometry
                    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                            shape_functions,                    // N
                            shape_derivatives,                  // DN
                            w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                            GeometryData::GI_GAUSS_1);          // first order Gauss integration

                    for ( unsigned int i = 0; i < w_gauss_neg_side.size(); i++){
                        neg_vol += w_gauss_neg_side[i];
                    }
                }
            }
            return neg_vol;
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:

        ModelPart& mr_model_part;

	};
}  // namespace Kratos.

#endif // KRATOS_CALCULATE_CUT_AREA_H_INCLUDED  defined