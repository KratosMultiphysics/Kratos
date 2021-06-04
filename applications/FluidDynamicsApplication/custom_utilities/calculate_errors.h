#if !defined(KRATOS_CALCULATE_NORML_H_INCLUDED)
#define  KRATOS_CALCULATE_NORML_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <tuple>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{
	class CalculateNormL
	{
	public:

		KRATOS_CLASS_POINTER_DEFINITION(CalculateNormL);

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

		CalculateNormL(ModelPart& model_part)
			: mr_model_part(model_part)              //mr_model_part is saved as private variable (declared at the end of the file)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}


		~CalculateNormL()
		{}


 		double Calculate(Vector exact_pressure)
		{
			KRATOS_TRY

			//getting data for the _given geometry
            double err_L2 {0}, norm_L2 {0};
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); //looping the elements
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				Matrix shape_functions;
                GeometryType::ShapeFunctionsGradientsType shape_derivatives;
                const GeometryType & rGeom = ielem->GetGeometry();
                GeometryType::Pointer p_geom = ielem->pGetGeometry();
                unsigned int pt_count_neg = 0;
                double errL2_e {0}, normL2_e {0};

                const unsigned int NumNodes = 4;

                Vector distances_vector( NumNodes, 0.0 );
                Vector pressure_vector( NumNodes, 0.0 );
                Vector density_vector( NumNodes, 0.0 );
                std::vector<array_1d<double, NumNodes>> velocity_vector;

                for (unsigned int i = 0; i < NumNodes; i++){
                    distances_vector[i] = rGeom[i].FastGetSolutionStepValue(DISTANCE);
                    pressure_vector[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE);
                    density_vector[i] = rGeom[i].FastGetSolutionStepValue(DENSITY);
                    velocity_vector.push_back(rGeom[i].FastGetSolutionStepValue(VELOCITY));
                    if ( distances_vector[i] <= 0.0 ){
                        pt_count_neg++;
                    }
                }

                if ( pt_count_neg > 0 && pt_count_neg <  NumNodes){
                    // element is cut by the surface (splitting)

                    // Construct the modified shape fucntions utility
                    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
                    if (NumNodes == 4) {
                        p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, distances_vector);
                    } else if (NumNodes == 3) {
                        p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, distances_vector);
                    } else { KRATOS_ERROR << "The process can not be applied on this kind of element" << std::endl; }

                    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients, int_positive_side_sh_func_gradients;
                    Vector positive_side_weights, int_positive_side_weights;
                    Matrix positive_side_sh_func, int_positive_side_sh_func;

                    // Call the interface outwards normal area vector calculator
			        std::vector<Vector> positive_side_area_normals;

			        p_modified_sh_func->ComputePositiveSideInterfaceAreaNormals(
			            positive_side_area_normals,
			            GeometryData::GI_GAUSS_2);

                    p_modified_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                        int_positive_side_sh_func,
                        int_positive_side_sh_func_gradients,
                        int_positive_side_weights,
                        GeometryData::GI_GAUSS_2);

                    unsigned int NGauss = int_positive_side_sh_func.size1();

                    //ENERGY OF THE POSITIVE SIDE
                    p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                        positive_side_sh_func,
                        positive_side_sh_func_gradients,
                        positive_side_weights,
                        GeometryData::GI_GAUSS_2);
                    NGauss = positive_side_sh_func.size1();

                    for (unsigned int g = 0; g < NGauss; g++)
                    {
                        array_1d<double,NumNodes> v_gauss_p = ZeroVector(NumNodes); //(0*rGeom[0].FastGetSolutionStepValue(VELOCITY);
                        array_1d<double,NumNodes> p_grad_p = 0*rGeom[0].FastGetSolutionStepValue(VELOCITY);
                        double cexact_gauss {0}, cnum_gauss {0};

                        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                        {
                            //div_vp_p       += positive_side_sh_func_gradients[g](iNode,0)*velocity_vector.at(iNode)(0)+ positive_side_sh_func_gradients[g](iNode,1)*velocity_vector.at(iNode)(1);
                            //v_gauss_p      += positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(VELOCITY);
                            cexact_gauss  += positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
                            cnum_gauss    += positive_side_sh_func(g, iNode)*exact_pressure[rGeom[iNode].Id()];
                        }

                        errL2_e  += (cnum_gauss - cexact_gauss)*(cnum_gauss - cexact_gauss)*positive_side_weights(g); 
                        normL2_e += (cexact_gauss)*(cexact_gauss)*positive_side_weights(g); 
                    }
                } else {
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                    const unsigned int NumGauss = IntegrationPoints.size();

                    Matrix N = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                    Vector DetJ;
                    GeometryType::ShapeFunctionsGradientsType DN_DX;
                    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,DetJ,GeometryData::GI_GAUSS_2);

                    for (unsigned int g = 0; g < NumGauss; g++)
                    {
                        double Weight = DetJ(g) * IntegrationPoints[g].Weight(); //"Jacobian" is 2.0*A for triangles
                        double exact_gauss {0}, num_gauss {0};

                        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                        {
                            num_gauss   += N(iNode,g)*rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
                            exact_gauss += N(iNode,g)*exact_pressure[rGeom[iNode].Id()];
                        }

                        errL2_e  += (num_gauss - exact_gauss)*(num_gauss - exact_gauss)*Weight; 
                        normL2_e += (exact_gauss)*(exact_gauss)*Weight; 
                    }
                }

                err_L2 += errL2_e;
                norm_L2 += normL2_e;
            }
            err_L2 = std::pow(err_L2, 0.5)/std::pow(err_L2, 0.5);
            return  err_L2;
			KRATOS_CATCH("")
		}


	protected:


	private:

        ModelPart& mr_model_part;

	};
}  // namespace Kratos.

#endif // KRATOS_CALCULATE_NORMAL_VECTOR_H_INCLUDED  defined
