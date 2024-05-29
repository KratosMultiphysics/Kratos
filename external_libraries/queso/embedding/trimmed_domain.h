// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef TRIMMED_DOMAIN_INCLUDE_H
#define TRIMMED_DOMAIN_INCLUDE_H

/// STL includes
#include <memory>
/// Project includes
#include "embedding/geometry_query.h"
#include "utilities/mesh_utilities.h"
#include "embedding/trimmed_domain_on_plane.h"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  TrimmedDomain
 * @author Manuel Messmer
 * @brief  Provides geometrical operations e.g. to compute closed triangle mesh of trimmed domain (based on a clipped triangle mesh).
 * @details Uses AABB Tree for fast search.
*/
class TrimmedDomain {

public:
    ///@name Type Definitions
    ///@{

    typedef Unique<TriangleMeshInterface> TriangleMeshPtrType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor for trimmed domain.
    /// @param pTriangleMesh Clipped triangle mesh. Must contain edges on boundary planes of aabb (see: TrimmedomainOnPlane).
    /// @param rLowerBound Lower bound of trimmed domain.
    /// @param rUpperBound Upper bound of trimmed domain.
    /// @param pOperator Pointer to BrepOperator to perform IsInside()-check. Only used as safety procedure, if IsInsideTrimmedDomain() failes.
    /// @param rParameters QuESo Parameters.
    /// @param SwitchPlaneOrientation If true, orientation of edges on TrimmedDomainOnPlane are switched.
    TrimmedDomain(TriangleMeshPtrType pTriangleMesh, const PointType& rLowerBound, const PointType& rUpperBound,
            const BRepOperator* pOperator, IndexType MinNumberOfTriangles = 100, bool SwitchPlaneOrientation = false )
        : mpTriangleMesh(std::move(pTriangleMesh)), mLowerBound(rLowerBound), mUpperBound(rUpperBound),
          mpClippedMesh(mpTriangleMesh->Clone()), mGeometryQuery(*mpClippedMesh, false)
    {
        // Set relative snap tolerance.
        mSnapTolerance = RelativeSnapTolerance(mLowerBound, mUpperBound);

        // Construct trimmed domain on plane upper bound of AABB.
        bool upper_bound = true;
        auto p_trimmed_domain_upper_x = MakeUnique<TrimmedDomainOnPlane>(0, upper_bound, mLowerBound, mUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_upper_y = MakeUnique<TrimmedDomainOnPlane>(1, upper_bound, mLowerBound, mUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_upper_z = MakeUnique<TrimmedDomainOnPlane>(2, upper_bound, mLowerBound, mUpperBound, this, SwitchPlaneOrientation);
        // Construct trimmed domain on plane lower bound of AABB.
        upper_bound = false;
        auto p_trimmed_domain_lower_x = MakeUnique<TrimmedDomainOnPlane>(0, upper_bound, mLowerBound, mUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_lower_y = MakeUnique<TrimmedDomainOnPlane>(1, upper_bound, mLowerBound, mUpperBound, this, SwitchPlaneOrientation);
        auto p_trimmed_domain_lower_z = MakeUnique<TrimmedDomainOnPlane>(2, upper_bound, mLowerBound, mUpperBound, this, SwitchPlaneOrientation);

        if( mpTriangleMesh->NumOfTriangles() > 0 ){
            auto p_t1 = p_trimmed_domain_lower_x->pGetTriangulation( *(mpTriangleMesh.get()), pOperator );
            auto p_t2 = p_trimmed_domain_upper_x->pGetTriangulation( *(mpTriangleMesh.get()), pOperator );
            auto p_t3 = p_trimmed_domain_lower_y->pGetTriangulation( *(mpTriangleMesh.get()), pOperator );
            auto p_t4 = p_trimmed_domain_upper_y->pGetTriangulation( *(mpTriangleMesh.get()), pOperator );
            auto p_t5 = p_trimmed_domain_lower_z->pGetTriangulation( *(mpTriangleMesh.get()), pOperator );
            auto p_t6 = p_trimmed_domain_upper_z->pGetTriangulation( *(mpTriangleMesh.get()), pOperator );

            const IndexType num_triangles = p_t1->NumOfTriangles() + p_t2->NumOfTriangles() + p_t3->NumOfTriangles()
                + p_t4->NumOfTriangles() + p_t5->NumOfTriangles() + p_t6->NumOfTriangles();

            mpTriangleMesh->Reserve(2UL*num_triangles);

            MeshUtilities::Append(*(mpTriangleMesh), *(p_t1));
            MeshUtilities::Append(*(mpTriangleMesh), *(p_t2));
            MeshUtilities::Append(*(mpTriangleMesh), *(p_t3));
            MeshUtilities::Append(*(mpTriangleMesh), *(p_t4));
            MeshUtilities::Append(*(mpTriangleMesh), *(p_t5));
            MeshUtilities::Append(*(mpTriangleMesh), *(p_t6));

            MeshUtilities::Refine(*(mpTriangleMesh.get()), MinNumberOfTriangles);
        }
    }

    ///@}
    ///@name Operations
    ///@{

    ///@brief Returns true if point is inside TrimmedDomain. Returns false, if test was not successful.
    ///@note Calls: IsInsideTrimmedDomain(rPoint, rSuccess). If you want to know if test was successful, you should directly call
    ///      IsInsideTrimmedDomain(rPoint, rSuccess).
    ///@param rPoint
    ///@return bool
    bool IsInsideTrimmedDomain(const PointType& rPoint) const {
        bool success = true;
        return IsInsideTrimmedDomain(rPoint, success);
    }

    ///@brief Returns true if point is inside TrimmedDomain. Expects point to be inside AABB. Check is omitted.
    ///@brief Performs ray tracing in direction of the first triangle. Search for all intersection of ray. Inside/Outside is detected
    ///       based on the orientation of the closest intersected triangle (forward or backward facing).
    ///@param rPoint
    ///@param[out] rSuccess Is set to false, if e.g. all triangles are detected as parallel and therefore give ambiguous results.
    ///@return bool
    bool IsInsideTrimmedDomain(const PointType& rPoint, bool& rSuccess) const;

    /// @brief Return ptr to Triangle mesh (Raw Ptr)
    /// @return const TriangleMeshInterface*
    const TriangleMeshInterface* const pGetTriangleMesh() const{
        return mpTriangleMesh.get();
    }

    /// @brief Return reference to triangle mesh
    /// @return const TriangleMeshInterface&
    const TriangleMeshInterface& GetTriangleMesh() const{
        return *(mpTriangleMesh.get());
    }

    ///@brief Triangulates trimmed domain (Surface mesh of outer hull) and return boundary integration points.
    /// @tparam BoundaryIntegrationPointType
    ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
    template<typename BoundaryIntegrationPointType>
    Unique<std::vector<BoundaryIntegrationPointType>> pGetBoundaryIps() const {
        // Pointer to boundary integration points
        auto p_boundary_ips = MakeUnique<std::vector<BoundaryIntegrationPointType>>();

        p_boundary_ips->reserve(mpTriangleMesh->NumOfTriangles()*12UL);
        for( IndexType triangle_id = 0; triangle_id < mpTriangleMesh->NumOfTriangles(); ++triangle_id ){
            IndexType method = 3; // Creates 12 points per triangle.
            auto p_new_points = mpTriangleMesh->pGetIPsGlobal<BoundaryIntegrationPointType>(triangle_id, method);
            p_boundary_ips->insert(p_boundary_ips->end(), p_new_points->begin(), p_new_points->end());
        }

        return p_boundary_ips;
    }

    ///@brief Returns intersections state of AABB. This is an interface for the octree.
    ///@note This test is only performed on the mClippedMesh to be more efficient.
    ///@param rLowerBound Lower bound of AABB to be tested. Expected to be inside trimmed domain.
    ///@param rUpperBound Upper bound of AABB to be tested. Expected to be inside trimmed domain.
    ///@param Tolerance Default: SNAPTOL. Tolerance reduces AABB slightly. If Tolerance=0 touch is detected as intersection.
    ///                 If Tolerance>0, touch is not detected as intersection.
    ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
    IntersectionStatusType GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance=SNAPTOL) const;

    /// @brief Returns bounding box of trimmed domain. (Might be smaller than the actual domain of element.)
    /// @return BoundingBox (std::pair: first - lower_bound, second - upper_bound)
    const BoundingBoxType GetBoundingBoxOfTrimmedDomain() const;

    ///@}
private:

    ///@}
    ///@name Private Members
    ///@{
    TriangleMeshPtrType mpTriangleMesh;
    PointType mLowerBound;
    PointType mUpperBound;

    Unique<TriangleMeshInterface> mpClippedMesh;
    GeometryQuery mGeometryQuery;
    double mSnapTolerance;

    ///@}
};

///@}

} // End namespace queso

#endif // TRIMMED_DOMAIN_INCLUDE_H