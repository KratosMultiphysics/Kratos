/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_SMAGORINSKY_INCLUDED )
#define  KRATOS_SMAGORINSKY_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
 	

	class SmagorinskyTurbulentModel
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		//
		template<unsigned int TDim>
		void CalculateTurbulentViscosity(ModelPart& ThisModelPart, double MolecularViscosity, double Cs)
		{			
			KRATOS_TRY;


			
			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
			array_1d<double,TDim+1> N;
			boost::numeric::ublas::bounded_matrix<double,TDim,TDim> dv_dx;
			boost::numeric::ublas::bounded_matrix<double,TDim,TDim> S;
			const unsigned int nnodes = TDim+1;
			//unsigned int dim = TDim;
			double lumping_coeff = 1.0/double(TDim+1);

			//compute the projection (first step
			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			{
				it->FastGetSolutionStepValue(VISCOSITY) = 0.0;
			}				



			//compute the projection (first step
			for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin(); 
				i!=ThisModelPart.ElementsEnd(); i++)
			{	
				//calculating shape functions values
				Geometry< Node<3> >& geom = i->GetGeometry();
				double Volume;
				
				GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

				//constructing the velocity gradient 
				array_1d<double,3>& v = geom[0].FastGetSolutionStepValue(VELOCITY);


				//calculate grad v
				for(unsigned int i = 0; i<TDim; i++)
					for(unsigned int j = 0; j<TDim; j++)
						dv_dx(i,j) = DN_DX(0,j) * v[i];
				for(unsigned int I = 1; I<nnodes; I++)
				{
					array_1d<double,3>& v = geom[I].FastGetSolutionStepValue(VELOCITY);

					for(unsigned int i = 0; i<TDim; i++)
						for(unsigned int j = 0; j<TDim; j++)
							dv_dx(i,j) += DN_DX(I,j) * v[i];
				}

				//calculate S y normS

				for(unsigned int i = 0; i<TDim; i++)
					for(unsigned int j = 0; j<TDim; j++)
						S(i,j) = 0.5 * (dv_dx(i,j) + dv_dx(j,i));

				double normS = 0.0;
				for(unsigned int i = 0; i<TDim; i++)
					for(unsigned int j = 0; j<TDim; j++)
						normS += 2 * S(i,j) * S(i,j);
				normS = sqrt(normS);


				//calculate eddy viscosity (of the element)
				double delta = pow(Volume, 1.0/3.0);
				double elemental_viscosity = pow((Cs * delta), 2);
				elemental_viscosity *= normS;

				
					
				//adding the force to the nodes
				for(unsigned int I = 0; I<nnodes; I++)
				{
					geom[I].FastGetSolutionStepValue(VISCOSITY) += (lumping_coeff * Volume ) * elemental_viscosity;
				}

			}

			//divide by the nodal area
			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			{
				it->FastGetSolutionStepValue(VISCOSITY) *= it->FastGetSolutionStepValue(DENSITY)/it->FastGetSolutionStepValue(NODAL_MASS);
				it->FastGetSolutionStepValue(VISCOSITY) += MolecularViscosity;
			}

			KRATOS_CATCH("")
		}



	private:


	};

}  // namespace Kratos.

#endif // KRATOS_SMAGORINSKY_INCLUDED  defined 


