//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2014-09-25 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_CLUSTER3D_H_INCLUDED )
#define  KRATOS_CLUSTER3D_H_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cmath>

// External includes 
//#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "containers/weak_pointer_vector.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"
#include "custom_elements/spheric_particle.h"
#include "custom_utilities/create_and_destroy.h"
#include "utilities/quaternion.h"

namespace Kratos
{
    
    class Cluster3D : public Element {
        
    public:
        /// Pointer definition of Cluster3D
        KRATOS_CLASS_POINTER_DEFINITION(Cluster3D);
       
        Cluster3D( );
        Cluster3D( IndexType NewId, GeometryType::Pointer pGeometry );
        Cluster3D( IndexType NewId, NodesArrayType const& ThisNodes);
        Cluster3D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;      

        /// Destructor.
        virtual ~Cluster3D();
      
        using Element::Initialize;
        virtual void Initialize(ProcessInfo& r_process_info);
        virtual void SetIntegrationScheme(DEMIntegrationScheme::Pointer& integration_scheme);
        virtual void InitializeSolutionStep(ProcessInfo& r_process_info) override {};
        virtual void FinalizeSolutionStep(ProcessInfo& r_process_info) override {};
        virtual void CustomInitialize(ProcessInfo& r_process_info);
        virtual void SetOrientation(const Quaternion<double> Orientation);
        virtual void CreateParticles(ParticleCreatorDestructor* p_creator_destructor, ModelPart& dem_model_part, PropertiesProxy* p_fast_properties, const bool continuum_strategy);
        virtual void UpdatePositionOfSpheres();
        virtual void UpdateLinearDisplacementAndVelocityOfSpheres();
        virtual void GetClustersForce(const array_1d<double,3>& gravity);
        virtual void CollectForcesAndTorquesFromSpheres();
        virtual void ComputeAdditionalForces(const array_1d<double,3>& gravity);
        unsigned int GetNumberOfSpheres() { return mListOfSphericParticles.size(); };
        std::vector<SphericParticle*>  GetSpheres() { return mListOfSphericParticles; }; 
        virtual void SetContinuumGroupToBreakableClusterSpheres(const int Id);
        virtual void SetInitialConditionsToSpheres(const array_1d<double,3>& velocity);
        virtual void SetInitialNeighbours(const double search_increment);
        virtual void CreateContinuumConstitutiveLaws();
        virtual void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) override;
        
        virtual void Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag);
        virtual DEMIntegrationScheme& GetIntegrationScheme() { return *mpIntegrationScheme; }
           
        virtual double GetMass();
        virtual double SlowGetDensity();
        virtual int SlowGetParticleMaterial();

        virtual std::string Info() const override
        {
	    std::stringstream buffer;
	    buffer << "Discrete Element #" << Id();
	    return buffer.str();
        }
      
        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override
        {
	    rOStream << "Discrete Element #" << Id();
        }
      
        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override
        {
	    //mpGeometry->PrintData(rOStream);
        }
 
    protected:
       
        std::vector<double>                mListOfRadii;
        std::vector<array_1d<double, 3> >  mListOfCoordinates;        
        std::vector<SphericParticle*>      mListOfSphericParticles; 
        DEMIntegrationScheme* mpIntegrationScheme;        
      
    private:
       
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        }

    }; // Class Cluster3D
   
    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, Cluster3D& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const Cluster3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
}  // namespace Kratos.

#endif // KRATOS_CLUSTER3D_INCLUDED  defined
