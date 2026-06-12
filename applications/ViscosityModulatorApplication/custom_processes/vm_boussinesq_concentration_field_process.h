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

#if !defined(KRATOS_VM_BOUSSINESQ_CONCENTRATION_FIELD_PROCESS_H_INCLUDED )
#define  KRATOS_VM_BOUSSINESQ_CONCENTRATION_FIELD_PROCESS_H_INCLUDED



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

  /// Auxiliary process to set Boussinesq buoyancy forces in variable density flows due to a non-constant concentration field.
  class KRATOS_API(VISCOSITY_MODULATOR_APPLICATION) VmBoussinesqConcentrationFieldProcess: public Process
  {
  public:
      KRATOS_CLASS_POINTER_DEFINITION(VmBoussinesqConcentrationFieldProcess);

      /// Constructor
      VmBoussinesqConcentrationFieldProcess(ModelPart& rModelPart);
      VmBoussinesqConcentrationFieldProcess(ModelPart& rModelPart, Parameters& rParameters);

      /// Destructor.
      ~VmBoussinesqConcentrationFieldProcess() override;

      void Execute() override;

      void ExecuteInitialize() override;

      void ExecuteInitializeSolutionStep() override;

      void ExecuteFinalizeSolutionStep() override;

      std::string Info() const override;

      void PrintInfo(std::ostream& rOStream) const override;

      void PrintData(std::ostream& rOStream) const override;

    protected:

      void ValidateModelPart();

      void AssignBoussinesqForce();

    private:

      ModelPart& mrModelPart;

      array_1d<double,3> mrGravity;
      array_1d<double,3> mrR0;

      double mRho0;
      double mRhoMax;
      bool mModifyPressure;
      bool mModifyDensity;

      VmBoussinesqConcentrationFieldProcess& operator=(VmBoussinesqConcentrationFieldProcess const& rOther);

      VmBoussinesqConcentrationFieldProcess(VmBoussinesqConcentrationFieldProcess const& rOther);

    }; // Class VmBoussinesqConcentrationFieldProcess

  inline std::ostream& operator << (
      std::ostream& rOStream,
      const VmBoussinesqConcentrationFieldProcess& rThis);

}  // namespace Kratos.

#endif // KRATOS_VM_BOUSSINESQ_CONCENTRATION_FIELD_PROCESS_H_INCLUDED  defined
