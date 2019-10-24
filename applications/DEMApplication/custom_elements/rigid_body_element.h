// Created by: Salva Latorre, latorre@cimne.upc.edu

#if !defined KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED
#define KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"
#include "containers/global_pointers_vector.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"
#include "custom_utilities/create_and_destroy.h"
#include "utilities/quaternion.h"
#include "custom_conditions/RigidFace.h"
#include "../custom_strategies/schemes/dem_integration_scheme.h"
#include "../custom_strategies/schemes/symplectic_euler_scheme.h"


namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) RigidBodyElement3D : public Element {

    public:
        /// Pointer definition of RigidBodyElement3D
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RigidBodyElement3D);

        RigidBodyElement3D();
        RigidBodyElement3D(IndexType NewId, GeometryType::Pointer pGeometry);
        RigidBodyElement3D(IndexType NewId, NodesArrayType const& ThisNodes);
        RigidBodyElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~RigidBodyElement3D();

        using Element::Initialize;
        virtual void Initialize(ProcessInfo& r_process_info);
        virtual void SetIntegrationScheme(DEMIntegrationScheme::Pointer& translational_integration_scheme, DEMIntegrationScheme::Pointer& rotational_integration_scheme);
        virtual void CustomInitialize(ModelPart& rigid_body_element_sub_model_part);
        virtual void SetOrientation(const Quaternion<double> Orientation);
        virtual void UpdateLinearDisplacementAndVelocityOfNodes();
        virtual void UpdateAngularDisplacementAndVelocityOfNodes();
        virtual void GetRigidBodyElementsForce(const array_1d<double,3>& gravity);
        virtual void CollectForcesAndTorquesFromTheNodes();
        virtual void ComputeExternalForces(const array_1d<double,3>& gravity);
        virtual void SetInitialConditionsToNodes(const array_1d<double,3>& velocity);

        virtual void Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag);
        virtual DEMIntegrationScheme& GetTranslationalIntegrationScheme() { return *mpTranslationalIntegrationScheme; }
        virtual DEMIntegrationScheme& GetRotationalIntegrationScheme() { return *mpRotationalIntegrationScheme; }

        virtual double GetMass();

        double GetSqrtOfRealMass();

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

        std::vector<array_1d<double, 3> > mListOfCoordinates;
        std::vector<Node<3>::Pointer > mListOfNodes;
        DEMIntegrationScheme* mpTranslationalIntegrationScheme;
        DEMIntegrationScheme* mpRotationalIntegrationScheme;

    //protected:

        std::vector<RigidFace3D*> mListOfRigidFaces;
        array_1d<double,3> mInertias;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override;

        virtual void load(Serializer& rSerializer) override;

    }; // Class RigidBodyElement3D

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, RigidBodyElement3D& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const RigidBodyElement3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos

#endif // KRATOS_RIGID_BODY_ELEMENT_H_INCLUDED  defined
