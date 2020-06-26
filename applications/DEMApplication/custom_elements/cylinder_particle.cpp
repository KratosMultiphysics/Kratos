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

}  // namespace Kratos.

