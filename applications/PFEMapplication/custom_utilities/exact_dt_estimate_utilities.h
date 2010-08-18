/*
==============================================================================
KratosPFEMApplication 
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
//   Last Modified by:    $Author: kazem $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_EXACT_DT_ESTIMATE_UTILITIES_INCLUDED )
#define  KRATOS_EXACT_DT_ESTIMATE_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	///@} 
	///@name Type Definitions
	///@{ 


	///@} 
	///@name  Enum's
	///@{

	///@}
	///@name  Functions 
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
		calculate the nodal H for all the nodes depending on the min distance
		of the neighbouring nodes.

		lonely nodes are given the average value of the H
	*/

	class ExactDtEstimateUtilities 
		
	{
	public:
		///@name Type Definitions
		///@{


		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.


		///@}
		///@name Operations
		///@{

		double CubicExactDt(double dt_max, ModelPart& model_part, int dimension)
		{
			KRATOS_TRY

			double deltatime = dt_max;
			double machin_zero = 1.0e-10;

		  if(dimension == 2)
		    {
			array_1d<double,3> dxprim12, dx32, dxprim32, dv12, dv32, dvprim12;
			double B, A, C;

			for(ModelPart::ElementsContainerType::iterator i = model_part.ElementsBegin();
				i!=model_part.ElementsEnd(); i++)
			{
				//calculating velocity and displacement vectors
				Geometry< Node<3> >& geom = i->GetGeometry();
                                noalias(dv32) = geom[2].FastGetSolutionStepValue(VELOCITY);
                                dv32 -= geom[1].FastGetSolutionStepValue(VELOCITY);

                                noalias(dv12) = geom[0].FastGetSolutionStepValue(VELOCITY);
                                dv12 -= geom[1].FastGetSolutionStepValue(VELOCITY);

                                dvprim12[0] = dv12[1];
                                dvprim12[1] = -1*dv12[0];

                                dxprim12[0] = geom[0].Y() - geom[1].Y();
                                dxprim12[1] = geom[1].X() - geom[0].X();

                                dx32[0] = geom[2].X() - geom[1].X();
                                dx32[1] = geom[2].Y() - geom[1].Y();

                                dxprim32[0] = geom[1].Y() - geom[2].Y();
                                dxprim32[1] = geom[2].X() - geom[1].X();

                                //Calculate A,B and C to calculate zero´s of dt
                                B = inner_prod(dv32,dxprim12);
                                B += inner_prod(dv12,dxprim32);

                                A = inner_prod(dv32,dvprim12);

                                C = inner_prod(dx32, dxprim12);

                               double zerodt = dt_max;
                               double Zdt1;
                               double Zdt2;

                                if(A==0.0)
                                  {
                                    if(B != 0.0) zerodt = -1.0*C/B;

                                  }

                                else
                                    {
                                      double delta = 0.0;
                                      delta = B*B - 4.0*A*C;
                                      if(delta > 0.0)
                                      {
                                           Zdt1 = (-1*B + sqrt(delta))/(2.0*A);
                                           Zdt2 = (-1*B - sqrt(delta))/(2.0*A);

                                          //zerodt = -1.0;
                                          if(Zdt1 > machin_zero) zerodt = Zdt1;
                                          if(zerodt>Zdt2 && Zdt2 > machin_zero) zerodt = Zdt2;
                                      }

                                    }

                                //negative race means evert dt is acceptable
                                  if(zerodt < machin_zero) zerodt = dt_max;

                                //global check
                                if(zerodt < deltatime)
                                {
                                    deltatime = zerodt;
                                    KRATOS_WATCH("(((((((((((((((INSIDE ESTIMATE)))))))))))))))");
                                }
                        }
			
		    }

	      if(dimension == 3)
		    {
			for(ModelPart::ElementsContainerType::iterator i = model_part.ElementsBegin();
				i!=model_part.ElementsEnd(); i++)
			{	
		          Geometry< Node<3> >& geom = i->GetGeometry();
			  int str = 0;
			    for(int ii=0; ii < geom.size(); ++ii)
				if(geom[ii].FastGetSolutionStepValue(IS_STRUCTURE) == 1.0)
				      str++;
			      

			  if(geom.size() == 4 && str==0)
			    {
			      array_1d<double,3> DX1,DX2,DX3,V1,V2,V3;

                                DX1[0] = geom[1].X() - geom[0].X();
                                DX1[1] = geom[1].Y() - geom[0].Y();
                                DX1[2] = geom[1].Z() - geom[0].Z();

                                DX2[0] = geom[2].X() - geom[0].X();
                                DX2[1] = geom[2].Y() - geom[0].Y();
                                DX2[2] = geom[2].Z() - geom[0].Z();

                                DX3[0] = geom[3].X() - geom[0].X();
                                DX3[1] = geom[3].Y() - geom[0].Y();
                                DX3[2] = geom[3].Z() - geom[0].Z();

                                noalias(V1) = geom[1].FastGetSolutionStepValue(VELOCITY);
                                V1 -= geom[0].FastGetSolutionStepValue(VELOCITY);

                                noalias(V2) = geom[2].FastGetSolutionStepValue(VELOCITY);
                                V2 -= geom[0].FastGetSolutionStepValue(VELOCITY);

                                noalias(V3) = geom[3].FastGetSolutionStepValue(VELOCITY);
                                V3 -= geom[0].FastGetSolutionStepValue(VELOCITY);

				//calculate coeficients of cubic equation using tensors
				//C3*t³ + C2*t² + C1*t + C0 = 0.0
				array_1d<double,4>  Coeficients = ZeroVector(4);
				array_1d<double,3>  roots = ZeroVector(3);
				int num_Rroots;
				double zerodt = dt_max;

				int Kronecker_delta = 1.0;
	
				for( int i = 0; i < 3; ++i )
				  for( int j = 0; j < 3; ++j )
			            for( int k = 0; k < 3; ++k )
				      if( i!=j  &&  j!=k  &&  i!=k )
					{
					  Kronecker_delta = (j-i) * (k-i) * (k-j);
					  Kronecker_delta /= (abs(j-i)*abs(k-i)*abs(k-j));

					  Coeficients[0] += DX3[i]*DX1[j]*DX2[k]*Kronecker_delta;

					  Coeficients[1] += (DX3[i]*(DX1[j]*V2[k] + DX2[k]*V1[j]) + DX1[j]*DX2[k]*V3[i])*Kronecker_delta;

					  Coeficients[2] += (V3[i]*(DX1[j]*V2[k] + DX2[k]*V1[j]) + DX3[i]*V1[j]*V2[k])*Kronecker_delta;

					  Coeficients[3] += (V3[i]*V1[j]*V2[k])*Kronecker_delta;
					}
				  
				  CubicRoots(Coeficients, roots, num_Rroots);

				  //roots are in ascending order
				  if(num_Rroots == 1 && roots[0] > machin_zero )
				      zerodt = roots[0];

				  if(num_Rroots == 2)
				      {
					if(roots[0] > machin_zero)
					    zerodt=roots[0];
					else if(roots[1] > machin_zero)
					    zerodt=roots[1];
				      }

				  if(num_Rroots == 3)
				      {
					if(roots[0] > machin_zero)
					    zerodt=roots[0];
					else if(roots[1] > machin_zero)
					    zerodt=roots[1];
					else if(roots[2] > machin_zero)
					    zerodt=roots[2];
				      } 
// KRATOS_WATCH(Coeficients);	
// KRATOS_WATCH(roots)
// KRATOS_WATCH(num_Rroots);		    
// KRATOS_WATCH(zerodt);
			  if(zerodt < deltatime )
			    {

				//deltatime = zerodt;
				KRATOS_WATCH("(((((((((((((((INSIDE ESTIMATE)))))))))))))))");
				KRATOS_WATCH(roots)
				KRATOS_WATCH("Coordinates");
				    array_1d<double,3> first_out_pre;
				    MathUtils<double>::CrossProduct(first_out_pre,DX1 , DX2);
				    KRATOS_WATCH("pre volume");
				    KRATOS_WATCH(first_out_pre);
				    KRATOS_WATCH(DX1);
				    KRATOS_WATCH(DX2);
				    KRATOS_WATCH(DX3);

				    KRATOS_WATCH(inner_prod(DX3,first_out_pre));
				    DX1 += V1*deltatime;
				    DX2 += V2*deltatime;
				    DX3 += V3*deltatime;
				    array_1d<double,3> first_out;
				    MathUtils<double>::CrossProduct(first_out,DX1 , DX2);

				    KRATOS_WATCH("Volume calculated");
				    KRATOS_WATCH(inner_prod(DX3,first_out));

			    }
			  }//end of if 4 noded elements loop
			}//end of element loop
			
		    }//end of 3d loop

		  return deltatime;
			KRATOS_CATCH("")
		}



		///@}
		///@name Access
		///@{ 


		///@}
		///@name Inquiry
		///@{


		///@}      
		///@name Input and output
		///@{

		/// Turn back information as a string.
		virtual std::string Info() const
		{
			return "ExactDtEstimateUtilities";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "ExactDtEstimateUtilities";
		}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const
		{
		}


		///@}      
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 


		///@} 
		///@name Protected  Access 
		///@{ 


		///@}      
		///@name Protected Inquiry 
		///@{ 


		///@}    
		///@name Protected LifeCycle 
		///@{ 


		///@}

	private:
		///@name Static Member Variables 
		///@{ 


		///@} 
		///@name Member Variables 
	



	     //***********************************************************************
	    //***********************************************************************
	    void CubicRoots(array_1d<double,4> A, array_1d<double,3>&  roots, int& num_Rroots)
		      {
			    KRATOS_TRY
			    const double PI = 3.1415926535897932;
			    const double THIRD = 1./3.;
			    double W, P, Q, DIS, PHI;
			    array_1d<double,3> U;

			    //define cubic root as statement function
			    // In C, the function is defined outside of the cubic fct

			    // ====determine the degree of the polynomial ====

			    if (A[3] != 0.0)
			    {
				    //cubic problem
				    W = A[2]/A[3]*THIRD;
				    P = pow((A[1]/A[3]*THIRD - pow(W,2)),3);
				    Q = -.5*(2.0*pow(W,3)-(A[1]*W-A[0])/A[3] );
				    DIS = pow(Q,2)+P;
				    if ( DIS < 0.0 )
				    {
					    //three real solutions!
					    //Confine the argument of ACOS to the interval [-1;1]!
					    PHI = acos(fmin(1.0,fmax(-1.0,Q/sqrt(-P))));
					    P=2.0*pow((-P),(5.e-1*THIRD));
					    for (int i=0;i<3;i++)	U[i] = P*cos((PHI+2*(double(i))*PI)*THIRD)-W;
					    roots[0] = fmin(U[0], fmin(U[1], U[2]));
					    roots[1] = fmax(fmin(U[0], U[1]),fmax( fmin(U[0], U[2]), fmin(U[1], U[2])));
					    roots[2] = fmax(U[0], fmax(U[1], U[2]));
					    num_Rroots = 3;
				    }
				    else
				    {
					    // only one real solution!
					    DIS = sqrt(DIS);
					    roots[0] = CBRT(Q+DIS)+CBRT(Q-DIS)-W;
					    num_Rroots=1;
				    }
			    }
			    else if (A[2] != 0.0)
			    {
				    // quadratic problem
				    P = 0.5*A[1]/A[2];
				    DIS = pow(P,2)-A[0]/A[2];
				    if (DIS > 0.0)
				    {
					    // 2 real solutions
					    roots[0] = -P - sqrt(DIS);
					    roots[1] = -P + sqrt(DIS);
					    num_Rroots=2;
				    }
				    else
				    {
					    // no real solution
					    num_Rroots=0;
				    }
			    }
			    else if (A[1] != 0.0)
			    {
				    //linear equation
				    roots[0] =A[0]/A[1];
				    num_Rroots=1;
			    }
			    else
			    {
				    //no equation
				    num_Rroots=0;
			    }
		    /*
		      *     ==== perform one step of a newton iteration in order to minimize
		      *          round-off errors ====
		      */
			    for (int i=0;i<num_Rroots;i++)
			    {
				    roots[i] = roots[i] - (A[0]+roots[i]*(A[1]+roots[i]*(A[2]+roots[i]*A[3])))/(A[1]+roots[i]*(2.0*A[2]+roots[i]*3.0*A[3]));
			    //	printf("\n X inside cubic %.15e\n", X[i]);
			    }

			    KRATOS_CATCH("")
		  }
		//***********************************************************************
		//***********************************************************************

		  int signR(double Z)
		  {
			  int ret;
			  if (Z > 0.0)	ret = 1;
			  if (Z < 0.0)	ret = -1;
			  if (Z == 0.0)	ret =0;

			  return ret;
		  }
		//***********************************************************************
		//***********************************************************************
		  double CBRT(double Z)
		  {
			  double ret;
			  const double THIRD = 1./3.;
			  //define cubic root as statement function
			  //SIGN has different meanings in both C and Fortran
			  // Was unable to use the sign command of C, so wrote my own
			  // that why a new variable needs to be introduced that keeps track of the sign of
			  // SIGN is supposed to return a 1, -1 or 0 depending on what the sign of the argument is
			  ret = fabs(pow(fabs(Z),THIRD)) * signR(Z);
			  return ret;
		  }
		///@} 
		///@name Private Operators
		///@{ 

		///@} 
		///@name Private Operations
		///@{ 


		///@} 
		///@name Private  Access 
		///@{ 


		///@}    
		///@name Private Inquiry 
		///@{ 


		///@}    
		///@name Un accessible methods 
		///@{ 

		/// Assignment operator.
		ExactDtEstimateUtilities& operator=(ExactDtEstimateUtilities const& rOther);

		/// Copy constructor.
		//ExactDtEstimateUtilities(ExactDtEstimateUtilities const& rOther);


		///@}    

	}; // Class ExactDtEstimateUtilities 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		ExactDtEstimateUtilities& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const ExactDtEstimateUtilities& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_H_PROCESS_INCLUDED  defined 


