// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

#if !defined(KRATOS_TOPOLOGY_UPDATING_UTILITIES_H_INCLUDED)
#define  KRATOS_TOPOLOGY_UPDATING_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <algorithm>

// External includes
#include <pybind11/pybind11.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"

#include "mma_solver.cpp"


// Application includes
#include "topology_optimization_application.h"


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

/// Solution utility that updates response values for next iteration.
/** Detail class definition.

 */

class TopologyUpdatingUtilities
{
public:

	///@name Type Definitions
	///@{

	/// Pointer definition of TopologyUpdatingUtilities
	KRATOS_CLASS_POINTER_DEFINITION(TopologyUpdatingUtilities);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	TopologyUpdatingUtilities( ModelPart& model_part )
	: mrModelPart(model_part)
	{
	}

	/// Destructor.
	virtual ~TopologyUpdatingUtilities()
	{
	}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	// ---------------------------------------------------------------------------------------------------------------------------------------------
	// --------------------------------- UPDATE DENSITIES  -----------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------------------------------------------------------------

	/// Finds the value of the X_PHYS (density) and updates it into the optimization problem
	void UpdateDensitiesUsingOCMethod( char update_type[], double volfrac, double greyscale , double OptItr , double qmax)
	{
		KRATOS_TRY;

		if ( strcmp( update_type , "oc_algorithm" ) == 0 ){
			clock_t begin = clock();
			std::cout << "  Optimality Criterion Method (OC) chosen to solve the optimization problem" << std::endl;

			// Check if Grey Scale Filter should be used
			double q = 1;
			if (greyscale == 1)
			{
				if (OptItr < 15)
					q = 1;
				else
					q = std::min(qmax, 1.01*q);

				std::cout << "  Grey Scale Filter activated, q = " << q << std::endl;
			}
			else
			{
				std::cout << "  Grey Scale Filter deactivated, q = " << q << std::endl;
			}

			// Update Densities procedure
			double l1     = 0.0;
			double l2     = 1000000000.0;
			double move   = 0.2;
			double sum_X_Phys;
			int nele;
			double x_new = 0.0;
			double lmid = 0.0;

			// Bisection algorithm to find Lagrange Multiplier so that volume constraint is satisfied (lmid)
			while ((l2-l1)/(l1+l2) > 0.001)
			{
				lmid = 0.5*(l2+l1);
				sum_X_Phys = 0.0;
				nele = 0;
				x_new = 0.0;

				for( ModelPart::ElementIterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
				{
					double x_old = element_i->GetValue(X_PHYS_OLD);
					int solid_void = element_i->GetValue(SOLID_VOID);
					double dcdx  = element_i->GetValue(DCDX);
					double dvdx  = element_i->GetValue(DVDX);

					// Update Density
					// When q = 1, Grey Scale Filter is not activated, i.e., the results are in the classical OC update method
					switch(solid_void)
					{
					// NORMAL elements
					case 0:
					{
						x_new = std::max(0.0, std::max(x_old - move, std::min(1.0, pow(std::min(x_old + move, x_old * sqrt(-dcdx/dvdx/lmid)),q))));
						break;
					}
					// ACTIVE elements (solid elements)
					case 1:
					{
						x_new = 1;
						break;
					}
					// PASSIVE elements (void elements)
					case 2:
					{
						x_new = 0;
						break;
					}
					default:
					{
						// If no element identification was found
						std::cout << "This value for SOLID_VOID does not exist."<< std::endl;
					}
					}

					// Update of the calculated X_PHYS for the next iteration
					element_i->SetValue(X_PHYS, x_new);

					// Updating additional quantities to determine the correct Lagrange Multiplier (lmid)
					sum_X_Phys = sum_X_Phys + x_new;
					nele = nele + 1;
				}

				if( sum_X_Phys > (volfrac*nele))
					l1=lmid;
				else
					l2=lmid;
			}

			// Printing of results
			clock_t end = clock();
			std::cout << "  Updating of values performed               [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;
		} else {
			KRATOS_ERROR << "No valid optimization_algorithm selected for the simulation. Selected one: " << update_type << std::endl;
		}

		KRATOS_CATCH("");

	}



	void UpdateDensitiesUsingMMAMethod( char update_type[], double volfrac,  double OptItr )
	{
		KRATOS_TRY;
		if ( strcmp( update_type , "MMA_algorithm" ) == 0 )
		{
			clock_t begin = clock();
			std::cout << "  Method of Moving Asymthodes (MMA) chosen to solve the optimization problem" << std::endl;


			// Create object of updating function
			int nn = mrModelPart.NumberOfElements();
			int mm = 1; /// constraints


			int iteration = 0; 

			Vector xold;
			xold.resize(mrModelPart.NumberOfElements());

			double *x = new double[nn];
			double *df = new double[nn];
			double *g = new double[mm];
			double *dg = new double[nn*mm];
			double *xmin = new double[nn];
			double *xmax = new double[nn];
			double Vol_sum =0;
			double vol_frac_iteration = 0;
			


			for( ModelPart::ElementsContainerType::iterator element_i = mrModelPart.ElementsBegin(); element_i!= mrModelPart.ElementsEnd(); element_i++ )
			{	
				
			
				double xval = element_i->GetValue(X_PHYS);
			/* 	double gx = volfrac; */
				double dfdx = element_i->GetValue(DCDX);
				double dgdx = element_i->GetValue(DVDX);
				double Xmin = 0;
				double Xmax= 1;
				Vol_sum += xval;

				
				x[iteration]= xval;
				df[iteration]= dfdx;
				dg[iteration] = dgdx;
				xmax[iteration] = Xmax;
				xmin[iteration] = Xmin; 
				iteration ++;
				

				

			}
			g[0] = 0;
			vol_frac_iteration = Vol_sum/nn;
			std::cout << "  Updating of values performed " <<  vol_frac_iteration << std::endl;
			g[0] = vol_frac_iteration  - volfrac ;
			std::cout << "  						g0 " <<  g[0] << std::endl;
			


 			// Initialize MMA
			MMASolver *mma = new MMASolver(nn,mm);

			double ch = 1;
			int itr = 0;

			while (ch > 0.002 && itr < 100)
			{	
				
				itr ++;
				double Xminn = 0;
				double Xmaxx= 1;
				double movlim = 0.1;
				double vol_summ = 0;

				for (int it = 0; it < nn; it++) 
				{
				xmax[it] = std::min(Xmaxx, x[it] + movlim);
				xmin[it] = std::max(Xminn, x[it] - movlim); 
				}

				mma->Update(x,df,g,dg,xmin,xmax);



				
				ch = 0.0;
				g[0]= 0;
				for (int itt=0;itt < nn;itt++) 
				{
				ch = std::max(ch,std::abs(x[itt]-xold[itt]));
			 	vol_summ += x[itt]; 
				xold[itt] = x[itt];
				}
				double vol_it_loop = vol_summ/nn;
				g[0] = vol_it_loop - vol_frac_iteration; 
			}




			int jiter = 0;

			for(ModelPart::ElementsContainerType::iterator elem_i = mrModelPart.ElementsBegin();
					elem_i!=mrModelPart.ElementsEnd(); elem_i++)
				elem_i->SetValue(X_PHYS, x[jiter++]);
				



			// Printing of results
			clock_t end = clock();
			std::cout << "  Updating of values performed               [ spent time =  " << double(end - begin) / CLOCKS_PER_SEC << " ] " << std::endl;

 /* 			delete[] x;
			delete[] df;
			delete[] g;
			delete[] dg;
			delete[] xmin;
			delete[] xmax;  */
			delete mma; 
		}
		else 
		{

			KRATOS_ERROR << "No valid optimization_algorithm selected for the simulation. Selected one: " << update_type << std::endl;
		}

		KRATOS_CATCH("");

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
		return "TopologyUpdatingUtilities";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "TopologyUpdatingUtilities";
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
	///@{

	ModelPart& mrModelPart;

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
	//TopologyUpdatingUtilities& operator=(TopologyUpdatingUtilities const& rOther);

	/// Copy constructor.
	//TopologyUpdatingUtilities(TopologyUpdatingUtilities const& rOther);


	///@}

}; // Class TopologyUpdatingUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_TOPOLOGY_UPDATING_UTILITIES_H_INCLUDED */
