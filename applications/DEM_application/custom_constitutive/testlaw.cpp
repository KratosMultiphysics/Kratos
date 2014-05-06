/* 
 * File:   DEMConstitutive.h
 *
 * Created on 21 de abril de 2014, 19:40
 */

//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (DEMCONSTITUTIVE_H_INCLUDED)
#define  DEMCONSTITUTIVE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/testlaw.h"
#include "DEM_application.h"
#include "testlaw.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TestLaw::TestLaw()

{
    cout << " TestLaw constructor..";

} // Class TestLaw

int main()

{
    TestLaw test;
    std::cout << " Testlaw print on screen..." << std::endl;
    return 0;
}

}  // namespace Kratos.
#endif //TestLaW_H  defined
