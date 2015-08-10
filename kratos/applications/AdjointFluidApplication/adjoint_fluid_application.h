/*
==============================================================================
KratosAdjointFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
(Released on march 05, 2007).

Copyright 2015
Mate Pentek, Michael Andre
mate.pentek@tum.de
michael.andre@tum.de
- Lehrstuhl fuer Statik, Technische Universitaet Muenchen, Arcisstrasse
21 80333 Munich, Germany

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
//   Project Name:        KratosAdjointFluidApplication $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:          February 2015 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_ADJOINT_FLUID_APPLICATION_H_INCLUDED)
#define  KRATOS_ADJOINT_FLUID_APPLICATION_H_INCLUDED

///@defgroup AdjointFluidApplication Adjoint Fluid Application
///@brief Adjoint solvers for computing drag sensitivities.

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes
#include "adjoint_fluid_application_variables.h"
#include "custom_elements/vms_adjoint_element.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// Main class of the Adjoint Fluid Application
class KratosAdjointFluidApplication : public KratosApplication
{
public:

  ///@name Type Definitions
  ///@{
    
  /// Pointer definition of KratosAdjointFluidApplication
  KRATOS_CLASS_POINTER_DEFINITION(KratosAdjointFluidApplication);

  ///@}
  ///@name Life Cycle
  ///@{

  KratosAdjointFluidApplication();

  virtual ~KratosAdjointFluidApplication() {}

  ///@}
  ///@name Operations
  ///@{

  virtual void Register();

  ///@}
  ///@name Input and output
  ///@{

  virtual std::string Info() const
  {
    return "KratosAdjointFluidApplication";
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const
  {
    rOStream << Info();
    PrintData(rOStream);
  }

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const
  {
    std::cout << "in KratosAdjointFluidApplication" << std::endl;
    KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
        KratosApplication::PrintData(rOStream);
  }

  ///@}

private:

  ///@name Member Variables
  ///@{

  const VMSAdjointElement<2> mVMSAdjointElement2D;
  const VMSAdjointElement<3> mVMSAdjointElement3D;

  ///@}
  ///@name Unaccessible methods
  ///@{
    
  KratosAdjointFluidApplication& operator=(KratosAdjointFluidApplication const& rOther);

  KratosAdjointFluidApplication(KratosAdjointFluidApplication const& rOther);

  ///@}
    
}; // class KratosAdjointFluidApplication

///@} Kratos classes

///@} AdjointFluidApplication group
  
} // namespace Kratos

#endif // KRATOS_ADJOINT_FLUID_APPLICATION_H_INCLUDED defined
