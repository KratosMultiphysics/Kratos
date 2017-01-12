/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-06-06 10:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_CONFIG_H_INCLUDED)
#define RVE_CONFIG_H_INCLUDED

#include "includes/constitutive_law.h"

namespace Kratos
{

/*
macros for the homogenizer class. changes the way the homogenized stress is computed.
0 = classical approach: volume avarage
1 = integral of boundary reactions tensor product with the position vector
2 = same as 1, but optimized performing the build_rhs only on boudary elements
*/
#define RVE_HOMOGENIZER_OPTIMIZATION 0 // sigma = 1/V(int(sigma_i))
//#define RVE_HOMOGENIZER_OPTIMIZATION 1 // sigma = 1/V(sum(boundary_F X))
//#define  RVE_HOMOGENIZER_OPTIMIZATION 2 // sigma = 1/V(sum(boundary_F X))(with build rhs reduced)


/*
macro for the rve linear system of equations class. activates timing (and prints results on cmd line)
*/
//#define RVE_SOE_TIMING
/*
macro for the rve linear system of equations class. does not asseble zero entries
*/
#define RVE_SOE_NO_ASM_ZEROS
#define RVE_SOE_SMALL_NUMBER 1.0E-10

/*
defines the iteration number before which the initial guess of the rve newton iteration is reverted to the previous step
*/
#define RVE_ITER_THRESHOLD_FOR_REVERT 1

/*
defines if the rve should generate random imperfaction on the micro-scale
*/
#define RVE_USE_MICRO_RANDOM_IMPERFECTION

/*
activates the use of rve predictor calculator
*/
//#define RVE_USE_PREDICTOR_CALCULATOR
#ifdef RVE_USE_PREDICTOR_CALCULATOR
#define RVE_PREDICTOR_INFO
#endif // RVE_USE_PREDICTOR_CALCULATOR

}


#endif // RVE_CONFIG_H_INCLUDED
