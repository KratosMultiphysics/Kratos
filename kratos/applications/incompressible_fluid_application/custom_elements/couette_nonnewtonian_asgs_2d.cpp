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
//   Last modified by:    $Author: antonia $
//   Date:                $Date: 2010-10-26 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//
 
//#define GRADPN_FORM
//#define STOKES

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/couette_nonnewtonian_asgs_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	CouetteNonNewtonianASGS2D::CouetteNonNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: NoNewtonianASGS2D(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
		
	}

	//************************************************************************************
	//************************************************************************************
	CouetteNonNewtonianASGS2D::CouetteNonNewtonianASGS2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: NoNewtonianASGS2D(NewId, pGeometry, pProperties)
	{
			
	}

	Element::Pointer CouetteNonNewtonianASGS2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		
		KRATOS_TRY
		return Element::Pointer(new CouetteNonNewtonianASGS2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
		KRATOS_CATCH("");
	}

	CouetteNonNewtonianASGS2D::~CouetteNonNewtonianASGS2D()
	{
	}



  //************************************************************************************
    //************************************************************************************

    void CouetteNonNewtonianASGS2D::CalculateApparentViscosity(double & app_mu, double & app_mu_derivative,
	    array_1d<double, 3 >&  grad_sym_vel, double & gamma_dot,
            const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,
            const double & mu) {
        KRATOS_TRY
        app_mu = 0.0;
// 	KRATOS_WATCH("COUETTE NON NEWTONIAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
	unsigned int nodes_number = 3;
	double yield = 0.0;
	double gamma_dot_inv;
        double mcoef = 300;
	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
	      yield +=  GetGeometry()[ii].FastGetSolutionStepValue(YIELD_STRESS);
	}
	yield /= nodes_number;
	double aux_1;
	CalculateGradSymVel(grad_sym_vel, gamma_dot, B);
	
        if (gamma_dot > 1e-10) {
            aux_1 = 1.0 - exp(-(mcoef * gamma_dot));
            app_mu = mu + (yield / gamma_dot) * aux_1;
// 			gamma_dot_inv = 1.0/gamma_dot;
            if (app_mu < mu) {
                KRATOS_ERROR(std::logic_error, "!!!!!!!!!!!  APPARENT VISCOSITY < VISCOSITY !!!!!!!!", this->Id());
            }
        } else {
            app_mu = mu + yield*mcoef ;
// 			gamma_dot_inv = 0.0;
        }
        
        if (gamma_dot <= 1e-10) gamma_dot_inv=1e10;
        else  gamma_dot_inv= 1.0/gamma_dot;
        
        app_mu_derivative = yield * gamma_dot_inv*(- gamma_dot_inv + exp(-(mcoef * gamma_dot))*(gamma_dot_inv + mcoef));
	
        KRATOS_CATCH("")
    }
	

	  

//     void CouetteNonNewtonianASGS2D::CalculateApparentViscosityStbl(double & app_mu, double & app_mu_derivative,
// 	    array_1d<double, 3 >&  grad_sym_vel, double & gamma_dot,
//             const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & B,
//             const double & mu) {
//         KRATOS_TRY
//         app_mu = 0.0;
// // 	KRATOS_WATCH("COUETTE NON NEWTONIAAAAAAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
// 	unsigned int nodes_number = 3;
// 	double yield = 0.0;
// // 	double gamma_dot_inv;
// 	for (unsigned int ii = 0; ii < nodes_number; ++ii) {
// 	      yield +=  GetGeometry()[ii].FastGetSolutionStepValue(YIELD_STRESS);
// 	}
// 	yield /= nodes_number;
// // 	double aux_1;
// 	CalculateGradSymVel(grad_sym_vel, gamma_dot, B);
// 	
// 
// ////////BILINEAR FORM
// 	double mu_s = 1e7;
// 	double gamma_dot_lim = yield/mu_s;
// 	
// 	if(gamma_dot >= gamma_dot_lim){
// 	  if (gamma_dot_lim <= 1e-16){
// 	    app_mu = mu ; //newtonian fluid, no rigid behavior. yield ~ 0;
// 	  }
// 	  else{
// 	    app_mu = (mu*(gamma_dot-gamma_dot_lim) + yield)/gamma_dot ;
// 	  }
// 	}
// 	else{
// 	  app_mu = mu_s ;
// 	}
//	
// 	
//        KRATOS_CATCH("")
//    }


} // Namespace Kratos


