//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
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
