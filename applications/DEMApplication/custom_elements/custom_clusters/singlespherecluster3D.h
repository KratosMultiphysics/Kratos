//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2014-09-25 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_SINGLE_SPHERE_CLUSTER3D_H_INCLUDED )
#define KRATOS_SINGLE_SPHERE_CLUSTER3D_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/indexed_object.h"
#include "containers/global_pointers_vector.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/create_and_destroy.h"

#include "includes/condition.h"
#include "custom_elements/cluster3D.h"

namespace Kratos
{
    class Element;
    class ProcessInfo;

    class SingleSphereCluster3D : public Cluster3D
    {
    public:

        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SingleSphereCluster3D);

        SingleSphereCluster3D( );
        SingleSphereCluster3D( IndexType NewId, GeometryType::Pointer pGeometry );
        SingleSphereCluster3D( IndexType NewId, NodesArrayType const& ThisNodes);
        SingleSphereCluster3D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

      /// Destructor.
        virtual ~SingleSphereCluster3D();

        void Initialize(const ProcessInfo& r_process_info) override;
        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& r_process_info) override;

        double SlowGetDensity() override;
        int SlowGetParticleMaterial() override;

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
        virtual void PrintData(std::ostream& rOStream) const override {}

    protected:

    private:

    }; // Class SingleSphereCluster3D

  /// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
                    SingleSphereCluster3D& rThis){
        return rIStream;
    }

  /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                    const SingleSphereCluster3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
}  // namespace Kratos.

#endif // KRATOS_SINGLE_SPHERE_CLUSTER3D_INCLUDED defined
