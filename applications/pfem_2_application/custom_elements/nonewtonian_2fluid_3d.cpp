//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/nonewtonian_2fluid_3d.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"



namespace Kratos
{



	//with matrixtype, constant coefficient
	void NoNewtonianMonolithicPFEM23D::AddViscousTerm(MatrixType& rDampMatrix,
                         const BoundedMatrix<double, 4, 3 >& rShapeDeriv,
                         double& Viscosity, const double Area)
	{
		double theta = 0.0;
		double Cohesion = 0.0;

        double base_viscosity = Viscosity;

		BoundedMatrix<double, 12,6 > B_matrix = ZeroMatrix(12,6);
		for (unsigned int i=0; i!=(4); i++) //i node
		{
			for (unsigned int j=0; j!=(3); j++) //x,y,z
				B_matrix(i*(3)+j,j)=rShapeDeriv(i,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(i*(3)+0,3)=rShapeDeriv(i,1);
			B_matrix(i*(3)+1,3)=rShapeDeriv(i,0);

			B_matrix(i*(3)+1,4)=rShapeDeriv(i,2);
			B_matrix(i*(3)+2,4)=rShapeDeriv(i,1);

			B_matrix(i*(3)+0,5)=rShapeDeriv(i,2);
			B_matrix(i*(3)+2,5)=rShapeDeriv(i,0);

		}

		int counter=0;
		BoundedMatrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);

		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		double pressure = 0.0;

		double negative_nodes=0.0;
		//double has_fixed_vel=false; //if some nodes are fixed, then we are on a boundary and then we can set different material properties (for example to simulate a lower basal friction angle)
		for (unsigned int i=0; i!=4; i++) //i node
		{
			if (this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0)
			{
				pressure += this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
				theta += this->GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE);
                Cohesion += this->GetGeometry()[i].FastGetSolutionStepValue(YIELD_STRESS);
				negative_nodes +=1.0;
			}
			//if (this->GetGeometry()[i].IsFixed(VELOCITY_X)==true)
			//    has_fixed_vel = true;
		}


		pressure /=negative_nodes;
		theta /=negative_nodes;
        Cohesion /=negative_nodes;

		if (pressure<0.0)
			pressure=0.0;
		double YieldStress = Cohesion + tan(theta) * pressure;
		//if (theta>basal_theta && has_fixed_vel)
		//        YieldStress = tan(basal_theta) * pressure;
		if (negative_nodes>0.5)
		{
			Viscosity = this->EffectiveViscosity(base_viscosity,YieldStress,rShapeDeriv);
			Viscosity*=negative_nodes/4.0; //number of nodes:
		}


		C_matrix *= Viscosity*Area;

		BoundedMatrix<double, 6 , 12  > temp_matrix = prod(C_matrix,trans(B_matrix));
		BoundedMatrix<double, 12 , 12  > viscosity_matrix = prod(B_matrix, temp_matrix );
		for (unsigned int i=0; i!=4; i++) //i node
		{
			for (unsigned int j=0; j!=4; j++) //j neighbour
			{
				for (unsigned int k=0; k!=3; k++) //xyz
					for (unsigned int l=0; l!=3; l++) //xyz
						rDampMatrix(i*4+k,j*4+l)+=viscosity_matrix(i*3+k,j*3+l);
			}
		}

		//OutputMatrix += viscosity_matrix;
	}


	double NoNewtonianMonolithicPFEM23D::EffectiveViscosity(double DynamicViscosity,
									  double YieldStress,
                                      const BoundedMatrix<double, 3+1, 3> &rDN_DX)
    {

        // Read the viscosity for the fluidified phase from the nodes
        // In Kratos, the viscosity is assumed to be given in kinematic units (m^2/s)
        double GammaDot = this->EquivalentStrainRate(rDN_DX);
        double m = 1.0e5;
        double OutputDynamicViscosity=DynamicViscosity;
        if (GammaDot > 1e-12) // Normal behaviour
        {
            double Regularization = 1.0 - std::exp(-m*GammaDot);
            OutputDynamicViscosity += Regularization * YieldStress / GammaDot;
        }
        else // fallback to avoid division by zero
        {
            // In this case dynamic viscosity goes to infinity,
            // understand the following as a large number times yield stress
            OutputDynamicViscosity += m*YieldStress;
        }
        if (OutputDynamicViscosity<0.0)
			KRATOS_WATCH(OutputDynamicViscosity);

        //this->GetValue(TEMPERATURE) = GammaDot;
        return OutputDynamicViscosity;
    }

	double NoNewtonianMonolithicPFEM23D::EquivalentStrainRate(const BoundedMatrix<double, 3+1, 3> &rDN_DX) // TDim+1,TDim
	{
		const int TDim=3;
		const GeometryType& rGeom = this->GetGeometry();
		const unsigned int NumNodes = rGeom.PointsNumber();
		// Calculate Symetric gradient
		BoundedMatrix<double,TDim,TDim> S = ZeroMatrix(TDim,TDim);
		for (unsigned int n = 0; n < NumNodes; ++n)
		{
			const array_1d<double,3>& rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY);
			for (unsigned int i = 0; i < TDim; ++i)
				for (unsigned int j = 0; j < TDim; ++j)
					S(i,j) += 0.5 * ( rDN_DX(n,j) * rVel[i] + rDN_DX(n,i) * rVel[j] );
		}
		// Norm of symetric gradient
		double NormS = 0.0;
		for (unsigned int i = 0; i < TDim; ++i)
			for (unsigned int j = 0; j < TDim; ++j)
				NormS += S(i,j) * S(i,j);
		return std::sqrt(2.0*NormS);
	}



} // Namespace Kratos
