//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//


// Project includes
#include "cylinder_particle.h"


namespace Kratos
{
     // using namespace GeometryFunctions;

      CylinderParticle::CylinderParticle()
      : SphericParticle(){/*mInitializedVariablesFlag = 0;*/}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : SphericParticle(NewId, pGeometry){/*mInitializedVariablesFlag = 0;*/}

      CylinderParticle::CylinderParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericParticle(NewId, pGeometry, pProperties){/*mInitializedVariablesFlag = 0;*/}

      CylinderParticle::CylinderParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericParticle(NewId, ThisNodes){/*mInitializedVariablesFlag = 0;*/}

      Element::Pointer CylinderParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return SphericParticle::Pointer(new CylinderParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      }

      /// Destructor.
      CylinderParticle::~CylinderParticle() {}

      double CylinderParticle::CalculateVolume(){
          return Globals::Pi * GetRadius() * GetRadius();
      }

      double CylinderParticle::CalculateMomentOfInertia() {
          return 0.5 * GetMass() * GetRadius() * GetRadius();
      }

      void CylinderParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info){}
      void CylinderParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info){}
      void CylinderParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
      void CylinderParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}

      void CylinderParticle::FinalizeStressTensor(ProcessInfo& r_process_info, double& rRepresentative_Volume){

        KRATOS_TRY
        SphericParticle::FinalizeStressTensor(r_process_info, rRepresentative_Volume);

        if (!r_process_info[IMPOSED_Z_STRAIN_OPTION]) return;

        double z_strain_value = r_process_info[IMPOSED_Z_STRAIN_VALUE];
        double myYoung = GetYoung();
        double myPoisson = GetPoisson();

        // (*mStressTensor)(2,2) += E*z_displacement - poisson*(sigma_xx + sigma_yy);
        (*mStressTensor)(2, 2) = myYoung*z_strain_value - myPoisson*((*mStressTensor)(0, 0) + (*mStressTensor)(1, 1));

        KRATOS_CATCH("")
    }

}  // namespace Kratos.

