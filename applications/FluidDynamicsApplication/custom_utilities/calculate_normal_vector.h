#if !defined(KRATOS_CALCULATE_NORMAL_VECTOR_H_INCLUDED)
#define  KRATOS_CALCULATE_NORMAL_VECTOR_H_INCLUDED


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
	class CalculateNormalVector
	{
	public:

		KRATOS_CLASS_POINTER_DEFINITION(CalculateNormalVector);

        typedef Node<3> NodeType;
        typedef Geometry<NodeType> GeometryType;

		CalculateNormalVector(ModelPart& model_part)
			: mr_model_part(model_part)              //mr_model_part is saved as private variable (declared at the end of the file)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}


		~CalculateNormalVector()
		{}


 		std::tuple<double, double> Calculate()
		{
			KRATOS_TRY

            double work_rate {0}, total_energy{0};
            const double A {1e5}, B {3e8}, gamma {7.15}, rho_0 {880};

			//getting data for the given geometry
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); //looping the elements
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				Matrix shape_functions;
                GeometryType::ShapeFunctionsGradientsType shape_derivatives;
                const GeometryType & rGeom = ielem->GetGeometry();
                GeometryType::Pointer p_geom = ielem->pGetGeometry();
                unsigned int pt_count_neg = 0;
                double elemental_energy {0};

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

                    for (unsigned int g = 0; g < NGauss; g++)
                    {
                        array_1d<double,NumNodes> v_gauss_i = 0*rGeom[0].FastGetSolutionStepValue(VELOCITY);
                        double p_gauss_i {0}, proj_i {0};
                        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                        {
                            v_gauss_i   += int_positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(VELOCITY);
                            p_gauss_i   += int_positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
                        }
                        proj_i = MathUtils<double>::Dot(v_gauss_i, positive_side_area_normals.at(g)) / norm_2(positive_side_area_normals.at(g));
                        work_rate += proj_i*p_gauss_i*int_positive_side_weights(g);
                    }

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
                        double p_gauss_p {0}, div_vp_p {0}, rho_gauss_p {0};

                        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                        {
                            div_vp_p       += positive_side_sh_func_gradients[g](iNode,0)*velocity_vector.at(iNode)(0)+ positive_side_sh_func_gradients[g](iNode,1)*velocity_vector.at(iNode)(1);
                            v_gauss_p      += positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(VELOCITY);
                            p_gauss_p      += positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
                            rho_gauss_p    += positive_side_sh_func(g, iNode)*rGeom[iNode].FastGetSolutionStepValue(DENSITY);
                        }

                        //double e = B*std::pow(rho_gauss_p/rho_0,gamma)/(gamma - 1) + B - A;
                        elemental_energy += positive_side_weights(g)*(0.5*MathUtils<double>::Dot(v_gauss_p, v_gauss_p)*rho_gauss_p + p_gauss_p);
                        //elemental_energy += positive_side_weights(g)*(MathUtils<double>::Dot(v_gauss_p, p_grad_p) + p_gauss_p*div_vp_p);
                    }
                } else if (pt_count_neg == 0){
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                    const unsigned int NumGauss = IntegrationPoints.size();

                    Matrix N = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
                    Vector DetJ;
                    GeometryType::ShapeFunctionsGradientsType DN_DX;
                    rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,DetJ,GeometryData::GI_GAUSS_2);

                    for (unsigned int g = 0; g < NumGauss; g++)
                    {
                        double Weight = DetJ(g) * IntegrationPoints[g].Weight(); //"Jacobian" is 2.0*A for triangles
                        array_1d<double,NumNodes> vgauss = 0*rGeom[0].FastGetSolutionStepValue(VELOCITY);
                        array_1d<double,NumNodes> p_grad = 0*rGeom[0].FastGetSolutionStepValue(VELOCITY);
                        double pgauss {0}, div_vel {0}, rho_gauss {0};

                        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                        {
                            div_vel   += DN_DX[g](iNode,0)*velocity_vector.at(iNode)(0)+ DN_DX[g](iNode,1)*velocity_vector.at(iNode)(1);
                            vgauss    += N(iNode,g)*rGeom[iNode].FastGetSolutionStepValue(VELOCITY);
                            pgauss    += N(iNode,g)*rGeom[iNode].FastGetSolutionStepValue(PRESSURE);
                            rho_gauss    += N(iNode,g)*rGeom[iNode].FastGetSolutionStepValue(DENSITY);
                        }

                        //double e = B*std::pow(rho_gauss/rho_0,gamma)/(gamma - 1) + B - A;
                        elemental_energy += Weight*(0.5*MathUtils<double>::Dot(vgauss, vgauss)*rho_gauss + pgauss);
                        // elemental_energy += Weight*(MathUtils<double>::Dot(vgauss, p_grad) + pgauss*div_vel);
                    }
                }

                total_energy += elemental_energy;
            }
            return  std::make_tuple(total_energy, work_rate);
			KRATOS_CATCH("")
		}


	protected:


	private:

        ModelPart& mr_model_part;

	};
}  // namespace Kratos.

#endif // KRATOS_CALCULATE_NORMAL_VECTOR_H_INCLUDED  defined
