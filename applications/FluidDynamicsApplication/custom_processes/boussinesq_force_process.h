//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#if !defined(KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED )
#define  KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
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

  /// Auxiliary process to set Boussinesq buoyancy forces in variable temperature flows.
  /** This process sets the BODY_FORCE variable to (1 - temp_fluctuation/reference_temp)*g,
      so that the fluid element takes buoyancy into account when solving the flow.
   */
  class BoussinesqForceProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of BoussinesqForceProcess
      KRATOS_CLASS_POINTER_DEFINITION(BoussinesqForceProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      BoussinesqForceProcess(ModelPart::Pointer pModelPart, Parameters& rParameters);

      /// Destructor.
      virtual ~BoussinesqForceProcess() override;


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      virtual void Execute() override;

      virtual void ExecuteInitialize() override;

      virtual void ExecuteInitializeSolutionStep() override;

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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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

      void ValidateModelPart();

      void AssignBoussinesqForce();

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

      ModelPart::Pointer mpModelPart;

      array_1d<double,3> mrGravity;

      bool mUseAmbientTemperature;

      double mThermalExpansionCoefficient;

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
      BoussinesqForceProcess& operator=(BoussinesqForceProcess const& rOther);

      /// Copy constructor.
      BoussinesqForceProcess(BoussinesqForceProcess const& rOther);


      ///@}

    }; // Class BoussinesqForceProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const BoussinesqForceProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUSSINESQ_FORCE_PROCESS_H_INCLUDED  defined
