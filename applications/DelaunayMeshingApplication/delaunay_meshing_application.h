//----------------------------------------------------------
//         ___      _                                      .
//  KRATOS|   \ ___| |__ _ _  _ _ _  __ _ _  _             .
//        | |) / -_| / _` | || | ' \/ _` | || |            .
//        |___/\___|_\__,_|\_,_|_||_\__,_|\_, |MESHING     .
//                                        |__/             .
//                                                         .
//  License:(BSD)   DelaunayMeshingApplication/license.txt .
//  Main authors:   Josep Maria Carbonell                  .
//                        ..                               .
//----------------------------------------------------------
//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_DELAUNAY_MESHING_APPLICATION_H_INCLUDED )
#define  KRATOS_DELAUNAY_MESHING_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes

// Core applications

//conditions
#include "custom_conditions/composite_condition.hpp"

#include "includes/kratos_application.h"

#include "delaunay_meshing_application_variables.h"

namespace Kratos
{
  ///@name Type	Definitions
  ///@{

  ///@name Kratos Globals
  ///@{

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
  class KRATOS_API(DELAUNAY_MESHING_APPLICATION) KratosDelaunayMeshingApplication : public KratosApplication
  {
  public:


    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosDelaunayMeshingApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosDelaunayMeshingApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosDelaunayMeshingApplication();

    /// Destructor.
    virtual ~KratosDelaunayMeshingApplication(){}


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
    std::string Info() const override
      {
	return "KratosDelaunayMeshingApplication";
      }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << Info();
      PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      KRATOS_WATCH( "in KratosDelaunayMeshingApplication" )
      KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
      rOStream << "Variables:" << std::endl;
      KratosComponents<VariableData>().PrintData(rOStream);
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
    ///@}
    ///@name Member Variables
    ///@{

    const CompositeCondition mCompositeCondition2D2N;
    const CompositeCondition mCompositeCondition3D3N;

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
    KratosDelaunayMeshingApplication& operator=(KratosDelaunayMeshingApplication const& rOther);

    /// Copy constructor.
    KratosDelaunayMeshingApplication(KratosDelaunayMeshingApplication const& rOther);

    ///@}

  }; // Class KratosDelaunayMeshingApplication

  ///@}


  ///@name Type Definitions
  ///@{
  ///@}
  ///@name Input and output
  ///@{
  ///@}


}  // namespace Kratos.

#endif // KRATOS_DELAUNAY_MESHING_APPLICATION_H_INCLUDED  defined


