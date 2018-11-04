//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

#if !defined(KRATOS_MESH_MOVING_APPLICATION_H_INCLUDED)
#define KRATOS_MESH_MOVING_APPLICATION_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/laplacian_meshmoving_element.h"
#include "custom_elements/structural_meshmoving_element.h"

#include "includes/mesh_moving_variables.h"
#include "includes/variables.h"

namespace Kratos {

///@name Kratos Globals
///@{

// Variables definition

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

/// Short class definition.
/** Detail class definition.
*/
class KratosMeshMovingApplication : public KratosApplication {
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of KratosMeshMovingApplication
  KRATOS_CLASS_POINTER_DEFINITION(KratosMeshMovingApplication);

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  KratosMeshMovingApplication();

  /// Destructor.
  ~KratosMeshMovingApplication() override {}

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  void Register() override;

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
  std::string Info() const override { return "KratosMeshMovingApplication"; }

  /// Print information about this object.
  void PrintInfo(std::ostream &rOStream) const override {
    rOStream << Info();
    PrintData(rOStream);
  }

  ///// Print object's data.
  void PrintData(std::ostream &rOStream) const override {
    KRATOS_WATCH("in my application");
    KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size());
    rOStream << "Variables:" << std::endl;
    KratosComponents<VariableData>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Elements:" << std::endl;
    KratosComponents<Element>().PrintData(rOStream);
    rOStream << std::endl;
    rOStream << "Conditions:" << std::endl;
    KratosComponents<Condition>().PrintData(rOStream);
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

  //       static const ApplicationCondition  msApplicationCondition;

  ///@}
  ///@name Member Variables
  ///@{
  const LaplacianMeshMovingElement mLaplacianMeshMovingElement2D3N;
  const LaplacianMeshMovingElement mLaplacianMeshMovingElement2D4N;
  const LaplacianMeshMovingElement mLaplacianMeshMovingElement3D4N;
  const LaplacianMeshMovingElement mLaplacianMeshMovingElement3D8N;
  const StructuralMeshMovingElement mStructuralMeshMovingElement2D3N;
  const StructuralMeshMovingElement mStructuralMeshMovingElement2D4N;
  const StructuralMeshMovingElement mStructuralMeshMovingElement3D4N;
  const StructuralMeshMovingElement mStructuralMeshMovingElement3D8N;
  const StructuralMeshMovingElement mStructuralMeshMovingElement3D6N;
  const StructuralMeshMovingElement mStructuralMeshMovingElement3D15N;
  const LaplacianMeshMovingElement mLaplacianMeshMovingElement;
  const StructuralMeshMovingElement mStructuralMeshMovingElement;
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
  KratosMeshMovingApplication &operator=(KratosMeshMovingApplication const &rOther);

  /// Copy constructor.
  KratosMeshMovingApplication(KratosMeshMovingApplication const &rOther);

  ///@}

}; // Class KratosMeshMovingApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_MESH_MOVING_APPLICATION_H_INCLUDED  defined
