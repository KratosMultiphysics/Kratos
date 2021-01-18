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
  /** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis
      so that the fluid element can take natural convection into account.

      This process makes use of the following data:
      - TEMPERATURE from the nodal solution step data: current temperature for the node (mandatory).
      - AMBIENT_TEMPERATURE from ProcessInfo: The reference temperature for the simulation (mandatory).
      - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
      - thermal_expansion_coefficient from the Parameters: a double defining the thermal expansion coefficient for the fluid (optional).

      With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
      The force is set to (1 + thermal_expansion_coefficient*(temperature - ambient_temperature) ) * g

      If the thermal expansion coefficient is not provided, it is assumed to be (1/ambient_temperature).
      This is the usual value for perfect gases (if the temperature is given in Kelvin).
   */
  class KRATOS_API(FLUID_DYNAMICS_APPLICATION) BoussinesqForceProcess: public Process
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
      BoussinesqForceProcess(ModelPart& rModelPart, Parameters& rParameters);

      /// Destructor.
      ~BoussinesqForceProcess() override;


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute() override;

      void ExecuteInitialize() override;

      void ExecuteInitializeSolutionStep() override;

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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


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

      ModelPart& mrModelPart;

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
