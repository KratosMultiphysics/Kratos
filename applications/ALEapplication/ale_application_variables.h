// =============================================================================
/*
 KratosALEApllication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Pooyan Dadvand, Riccardo Rossi, Andreas Winterstein
                     pooyan@cimne.upc.edu
                     rrossi@cimne.upc.edu
                     a.winterstein@tum.de
- CIMNE (International Center for Numerical Methods in Engineering),
  Gran Capita' s/n, 08034 Barcelona, Spain
- Chair of Structural Analysis, Technical University of Munich
  Arcisstrasse 21 80333 Munich, Germany

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
*/
//==============================================================================

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: Dec 2016 $
 *  Revision:            $Revision: 0.0 $
 * ***************************************************************************/


 #if !defined(KRATOS_ALE_APPLICATION_VARIABLES_H_INCLUDED)
 #define  KRATOS_ALE_APPLICATION_VARIABLES_H_INCLUDED

 // System includes

 // External includes


 // Project includes
 #include "includes/define.h"
 //#include "includes/kratos_application.h"
 #include "includes/variables.h"

 namespace Kratos
 {
   //MESH_VELOCITY currently put to the core since used in other applications
   //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_VELOCITY)
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);


 } // namespace Kratos.

 #endif // KRATOS_ALE_APPLICATION_VARIABLES_H_INCLUDED
