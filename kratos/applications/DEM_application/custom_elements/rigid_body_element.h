//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2017-01-09 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//

#if !defined(KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED)
#define KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED

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

namespace Kratos {
    
    class RigidBodyElement : public Element {
        
    public:
        /// Pointer definition of RigidBodyElement
        KRATOS_CLASS_POINTER_DEFINITION(RigidBodyElement);
       
        RigidBodyElement();
        RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry);
        RigidBodyElement(IndexType NewId, NodesArrayType const& ThisNodes);
        RigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;      

        /// Destructor
        virtual ~RigidBodyElement();
      
        using Element::Initialize;
        virtual void Initialize(ProcessInfo& r_process_info);
        virtual void SetIntegrationScheme(DEMIntegrationScheme::Pointer& integration_scheme);
        virtual void InitializeSolutionStep(ProcessInfo& r_process_info){};
        virtual void FinalizeSolutionStep(ProcessInfo& r_process_info){};
        virtual void CustomInitialize(ProcessInfo& r_process_info);
        virtual void SetOrientation(const Quaternion<double> Orientation);
        virtual void UpdatePositionOfNodes();
        virtual void UpdateLinearDisplacementAndVelocityOfNodes();
        virtual void GetRigidBodyElementForce(const array_1d<double,3>& gravity);
        virtual void CollectForcesAndTorquesFromNodes();
        virtual void ComputeAdditionalForces(const array_1d<double,3>& gravity);
        virtual void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info);
        
        virtual void Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag);
        virtual DEMIntegrationScheme& GetIntegrationScheme() { return *mpIntegrationScheme; }
   
        double GetSqrtOfRealMass();
        virtual double SlowGetDensity();
        virtual int SlowGetRigidBodyElementMaterial();

        virtual std::string Info() const
        {
	    std::stringstream buffer;
	    buffer << "Discrete Element #" << Id();
	    return buffer.str();
        }
      
        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
	    rOStream << "Discrete Element #" << Id();
        }
      
        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const
        {
	    //mpGeometry->PrintData(rOStream);
        }
 
    protected:

        DEMIntegrationScheme* mpIntegrationScheme;        
      
    private:
       
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        }

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        }

    }; // Class RigidBodyElement
   
    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, RigidBodyElement& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const RigidBodyElement& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
} // namespace Kratos

#endif // KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED  defined
