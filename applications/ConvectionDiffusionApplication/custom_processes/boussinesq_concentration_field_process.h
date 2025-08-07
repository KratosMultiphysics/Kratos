//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Aniol Sala Pascual
//

#if !defined(KRATOS_BOUSSINESQ_CONCENTRATION_FIELD_PROCESS_H_INCLUDED )
#define  KRATOS_BOUSSINESQ_CONCENTRATION_FIELD_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{
  ///@addtogroup ConvectionDiffusionApplication
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

  /// Auxiliary process to set Boussinesq buoyancy forces in variable density flows due to a non-constant concentration field.
  /** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis.
      so that the fluid element can take natural convection into account.

      This process makes use of the following data:
      - base_fluid_density from the Parameters in the constructor: double that defines the density of the base fluid (mandatory)
      - particles_density from the Parameters in the constructor: double that defines the density of the particles' suspension (mandatory)
      - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
      - modify_pressure from the Parameters passed in the constructor: bool that asserts if the constant term of the body force is absorbed within the pressure

      With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
      The force is set to (1 + (rho_p - rho_0)/rho_0 * phi) * g ,
      where rho_0 is the base fluid's density, rho_p is the particle's density, and phi is the concentration field.

   */
  class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) BoussinesqConcentrationFieldProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of BoussinesqConcentrationFieldProcess
      KRATOS_CLASS_POINTER_DEFINITION(BoussinesqConcentrationFieldProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      BoussinesqConcentrationFieldProcess(ModelPart& rModelPart);
      BoussinesqConcentrationFieldProcess(ModelPart& rModelPart, Parameters& rParameters);

      /// Destructor.
      ~BoussinesqConcentrationFieldProcess() override;


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
      array_1d<double,3> mrR0;

      double mRho0;
      double mRhoP;
      bool mModifyPressure;
      bool mModifyDensity;

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
      BoussinesqConcentrationFieldProcess& operator=(BoussinesqConcentrationFieldProcess const& rOther);

      /// Copy constructor.
      BoussinesqConcentrationFieldProcess(BoussinesqConcentrationFieldProcess const& rOther);


      ///@}

    }; // Class BoussinesqConcentrationFieldProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const BoussinesqConcentrationFieldProcess& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOUSSINESQ_CONCENTRATION_FIELD_PROCESS_H_INCLUDED  defined
