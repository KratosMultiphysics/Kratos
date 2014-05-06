/* 
 * File:   testlaw.h
 *
 * Created on 21 de abril de 2014, 19:40
 */

//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (TESTLAW_H_INCLUDED)
#define  TESTLAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"



namespace Kratos
{

/**
 ** Defines a test law that does nothing interesting
 **/

class TestLaw
{
protected:

public:

    TestLaw() { std::cout << "Law constructor..." ; }
    virtual ~TestLaw() { std::cout << "Law destructor..." ; }

    void Move() const { std::cout << "Move function..." ;}

}; // Class TestLAw
}  // namespace Kratos.
#endif // TESTLAW_H  defined
