/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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

/* *********************************************************
*
*   Last Modified by:    $Author: Nelson Lafontaine $
*   Date:                $Date: 26-06-2009 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

#include "custom_utilities/tensor_utils.h"
#include "includes/ublas_interface.h"
#include "includes/properties.h"
#include "soft_hard_behavior/exponencial_softening.h"
#include <cmath>
#include <string>
#include <iostream>


namespace Kratos
{

Exponential_Softening::Exponential_Softening():SofteningHardeningCriteria() {}
Exponential_Softening::~Exponential_Softening() {}


double Exponential_Softening::Calculate(const Vector& Imput_Parameters)
{
    const double& Ft   = (*mprops)[FT];
    const double& Ec   = (*mprops)[YOUNG_MODULUS];
    const double& GE   = (*mprops)[FRACTURE_ENERGY];
    const double& l    = Imput_Parameters[0];
    const double& r_o  = Imput_Parameters[1];  // Ft/sqrt(Ec);
    const double& r    = Imput_Parameters[2];

    const double Hs_barra = Ft * Ft /(  2.00 * Ec * GE );
    const double ls_barra = 1.00 / Hs_barra;
    double Hs       = l / (ls_barra - l);
    if (Hs < 0.00)
    {
        Hs  = 0.00;
        std::cout<<"Warning: Softening Parameters is less than zero. Please refine more."<<std::endl;
    }
    double elev = -2.00 * Hs *( (r-r_o) /r_o );
    double q_r  =  r_o*exp(elev);
    double a    =  1.00 - (q_r/r);

    return a;
}

double  Exponential_Softening::FunctionBehavior(const Vector& Imput_Parameters)
{

    /// Using K_punto
    const double& Ft   = (*mprops)[FT];
    const double& Kp   = Imput_Parameters[3];
    double value = 0.00;
    if(Kp<1.00)
        value    = (1.00 - Kp) * Ft;
    return value;

    /*
    const double& he    = Imput_Parameters[0]; /// Longituf del elemento
    const double& Ep    = Imput_Parameters[1]; /// Deformacion plastica efectiva
    const double& Ft    = (*mprops)[FT];
    const double& Ec    = (*mprops)[YOUNG_MODULUS];
    const double& GE    = (*mprops)[FRACTURE_ENERGY];
            const double H      = Ft*Ft*he/(2.00 * GE);
    const double elev   = -2.00 * H * (Imput_Parameters[1]/Ft);
    const double ft     = Ft * exp(elev);
    //const double param  = 0.30 *(*mprops)[FT];
    const double result = ft;  //> param ?  ft:param;
    return result;
    */

}

double  Exponential_Softening::FirstDerivateFunctionBehavior(const Vector& Imput_Parameters)
{


    /// Using K_punto
    const double& Ft   = (*mprops)[FT];
    const double Kp    = Imput_Parameters[3] - Imput_Parameters[2];
    double value = 0.00;
    if(Kp<1.00)
        value    = -Ft * Kp;
    return value;

    /*
    const double& he    = Imput_Parameters[0]; /// Longituf del elemento
    const double& Ep    = Imput_Parameters[1]; /// Deformacion plastica efectiva
    const double& Ft    = (*mprops)[FT];
    const double& Ec    = (*mprops)[YOUNG_MODULUS];
    const double& GE    = (*mprops)[FRACTURE_ENERGY];
            const double H      = Ft*Ft*he/(2.00 * GE);
    const double elev   = -2.00 * H * (Imput_Parameters[1]/Ft);
    //const double ft     = Ft * exp(elev);
    //const double param = 0.30 *(*mprops)[FT];
    //const double result= -2.00 * Hs * exp(elev);
    const double result = -2.00 * H * exp(elev);  //(ft>param) ? -2.00 * Hs * exp(elev): 0.00;
    return result;
    */

}


/**
 * definition of CONSTITUTIVE_LAW variable
 */
}  /* namespace Kratos.*/


