/*
==============================================================================
KratosMeshlessApplication 
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
//   Last modified by:    $Author:  $
//   Date:                $Date: $
//   Revision:            $Revision: 1.3 $
//
//

// System includes 

#if defined(KRATOS_PYTHON)
// External includes 
#include <boost/python.hpp>


// Project includes 
#include "includes/define.h"
#include "meshless_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosMeshlessApplication)
{

    class_<KratosMeshlessApplication,
            KratosMeshlessApplication::Pointer,
            bases<KratosApplication>, boost::noncopyable >("KratosMeshlessApplication")
            ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(EFFECTIVE_RADIUS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DENSITY_NORM_PARAM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DIV_OF_VEL);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(OLD_DENSITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DENS_VARIATION);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DENS_DIFF);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(SOUND_VELOCITY);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(VER_WALL_LEFT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(VER_WALL_RIGHT);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(HOR_WALL_BOTTOM);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(IS_WET);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(INI_PRESSURE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(OUT_OF_SYSTEM);



    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_NORMALIZE_RHS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_APPLY_XSPH);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_BOUNDARY_PRESSURES);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_CATCH_FREESURFACE);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DUMMY_INTERMEDIATE_RHS);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(DELTA_TIME_ISPH);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_ACC );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( VISCOUS_ACC );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( BODYFORCE_ACC );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( BOUNDARY_ACC );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( XSPH_VELOCITY );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TEMP_POS );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TEMP_VEL );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TEMP_RHS );
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( TEMP_DISP );

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS( OLD_VEL );


}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
