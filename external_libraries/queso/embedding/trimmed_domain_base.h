// // Author: Manuel Meßmer
// // Email: manuel.messmer@tum.de

// #ifndef TRIMMED_DOMAIN_BASE_INCLUDE_H
// #define TRIMMED_DOMAIN_BASE_INCLUDE_H

// //// STL includes
// #include <memory>
// #include <functional>
// //// Project includes
// #include "queso/containers/boundary_integration_point.hpp"
// #include "queso/containers/triangle_mesh.hpp"
// #include "queso/utilities/mesh_utilities.h"
// #include "queso/includes/parameters.h"

// namespace queso {

// ///@name QuESo Classes
// ///@{

// /// Forward declarations
// class Element;

// /**
//  * @class  TrimmedDomainBase
//  * @author Manuel Messmer
//  * @brief  Base class for TrimmedDomain. Stores boundary mesh of trimmed domain.
//  * @todo Remove trimmed domain base.
// */
// class TrimmedDomainBase {

// public:
//     ///@name Type Definitions
//     ///@{

//     typedef std::vector<BoundaryIntegrationPoint> BoundaryIPVectorType;
//     typedef Unique<BoundaryIPVectorType> BoundaryIPVectorPtrType;
//     typedef Unique<TriangleMesh> TriangleMeshPtrType;
//     typedef std::pair<PointType, PointType> BoundingBox;

//     ///@}
//     ///@name Life Cycle
//     ///@{

//     /// Constructor
//     TrimmedDomainBase(const PointType& rLowerBound, const PointType& rUpperBound)
//         : mLowerBound(rLowerBound), mUpperBound(rUpperBound)
//     {
//     }

//     /// Constructor
//     TrimmedDomainBase(TriangleMeshPtrType pTriangleMesh, const PointType& rLowerBound, const PointType& rUpperBound)
//         : mpTriangleMesh(std::move(pTriangleMesh)), mLowerBound(rLowerBound), mUpperBound(rUpperBound)
//     {
//     }

//     /// Destructor
//     virtual ~TrimmedDomainBase() = default;

//     ///@}
//     ///@name Operations
//     ///@{

//     ///@brief Returns true if point is inside TrimmedDomain. Returns false, if test was not successful.
//     ///@note Calls: IsInsideTrimmedDomain(rPoint, rSuccess). If you want to know if test was successful, you should directly call
//     ///      IsInsideTrimmedDomain(rPoint, rSuccess).
//     ///@param rPoint
//     ///@return bool
//     bool IsInsideTrimmedDomain(const PointType& rPoint) const {
//         bool success = true;
//         return IsInsideTrimmedDomain(rPoint, success);
//     }

//     ///@brief Returns true if point is inside TrimmedDomain.
//     ///@param rPoint
//     ///@param[out] rSuccess is true if operation was successful.
//     ///@return bool
//     virtual bool IsInsideTrimmedDomain(const PointType& rPoint, bool& rSuccess ) const = 0;

//     ///@brief Returns boundary integration points of TrimmedDomain.
//     ///@return BoundaryIPVectorPtrType. Boundary integration points to be used for ConstantTerms::Compute.
//     virtual BoundaryIPVectorPtrType pGetBoundaryIps() const = 0;

//     /// @brief Returns bounding box of trimmed domain. (Might be smaller than the actual domain of element.)
//     /// @return BoundingBox (std::pair: first - lower_bound, second - upper_bound)
//     virtual const BoundingBox GetBoundingBoxOfTrimmedDomain() const = 0;

//     /// @brief Return ptr to Triangle mesh (Raw Ptr)
//     /// @return const TriangleMesh*
//     const TriangleMesh* const pGetTriangleMesh() const{
//         return mpTriangleMesh.get();
//     }

//     /// @brief Return reference to triangle mesh
//     /// @return const TriangleMesh&
//     const TriangleMesh& GetTriangleMesh() const{
//         return *(mpTriangleMesh.get());
//     }

//     ///@brief Returns intersections state of AABB (This is an Interface for the Octree).
//     ///@param rLowerBound
//     ///@param rUpperBound
//     ///@param Tolerance Tolerance reduces AABB slightly. If Tolerance=0 touch is detected as intersection.
//     ///                 If Tolerance>0, touch is not detected as intersection.
//     ///@return IntersectionStatus, enum: (0-Inside, 1-Outside, 2-Trimmed).
//     virtual IntersectionStatusType GetIntersectionState(const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance=EPS0) const = 0;

//     /// @brief Returns part of triangle mesh that IsInDomain.
//     /// @param IsInDomain std::function
//     /// @return TriangleMeshPtrType (Unique)
//     TriangleMeshPtrType pGetTriangleMesh(std::function<bool(double, double,double)> &IsInDomain) const {
//         // Get Ids of all triangles that are inside given domain.
//         std::vector<IndexType> triangle_ids;
//         const IndexType num_triangles = mpTriangleMesh->NumOfTriangles();
//         triangle_ids.reserve(num_triangles);
//         for( IndexType triangle_id = 0; triangle_id < num_triangles; ++triangle_id ){
//             const auto& p1 = mpTriangleMesh->P1(triangle_id);
//             const auto& p2 = mpTriangleMesh->P2(triangle_id);
//             const auto& p3 = mpTriangleMesh->P3(triangle_id);
//             if( IsInDomain(p1[0], p1[1], p1[2]) ){
//                 if( IsInDomain(p2[0], p2[1], p2[2]) ){
//                     if( IsInDomain(p3[0], p3[1], p3[2]) ){
//                         triangle_ids.push_back(triangle_id);
//                     }
//                 }
//             }
//         }
//         // Copy all triangles in (triangle_ids) to new mesh.
//         auto p_new_mesh = MakeUnique<TriangleMesh>();
//         p_new_mesh->Reserve(mpTriangleMesh->NumOfTriangles());
//         MeshUtilities::Append(*p_new_mesh, *mpTriangleMesh, triangle_ids);
//         return p_new_mesh;
//     }

//     ///@}

// protected:
//     ///@name Protected member variables
//     ///@{
//     TriangleMeshPtrType mpTriangleMesh;
//     PointType mLowerBound;
//     PointType mUpperBound;

//     ///@}
// }; // End TrimmedDomainBase
// ///@}

// } // End namespace queso

// #endif // TRIMMED_DOMAIN_BASE_INCLUDE_H