//   ____        ______  _____
//  / __ \      |  ____|/ ____|
// | |  | |_   _| |__  | (___   ___
// | |  | | | | |  __|  \___ \ / _ \'
// | |__| | |_| | |____ ____) | (_) |
//  \___\_\\__,_|______|_____/ \___/
//         Quadrature for Embedded Solids
//
//  License:    BSD 4-Clause License
//              See: https://github.com/manuelmessmer/QuESo/blob/main/LICENSE
//
//  Authors:    Manuel Messmer

//// STL includes
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <random>
//// Project includes
#include "queso/embedding/brep_operator.h"
#include "queso/embedding/ray_aabb_primitive.h"


namespace queso {

typedef BRepOperator::TrimmedDomainPtrType TrimmedDomainPtrType;
typedef BRepOperator::StatusVectorType StatusVectorType;

bool BRepOperator::IsInside(const PointType& rPoint) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> drandon(0.5, 1.5);

    // Rough test, if point is actually within the bounding box of the mesh.
    if( mGeometryQuery.IsWithinBoundingBox(rPoint)) {
        bool success = false;
        IndexType iteration = 0UL;
        const IndexType max_iteration = 100UL;
        bool is_inside = false;
        while( !success ){
            if( iteration >= max_iteration){ return false; }
            iteration++;
            // Get random direction. Must be postive! -> x>0, y>0, z>0
            Vector3d direction{drandon(gen), drandon(gen), drandon(gen)};

            // Normalize
            const double norm_direction = Math::Norm( direction );
            Math::DivideSelf( direction, norm_direction);

            // Construct ray
            Ray_AABB_primitive ray(rPoint, direction);

            // Query ray
            std::tie(is_inside, success) = mGeometryQuery.IsInside(ray);
        }
        return is_inside;
    }
    return false;
}

bool BRepOperator::IsTrimmed(const PointType& rLowerBound,  const PointType& rUpperBound, double Tolerance) const {
    return mGeometryQuery.DoIntersect(rLowerBound, rUpperBound, Tolerance);
}


bool BRepOperator::OnBoundedSideOfClippedSection( const PointType& rPoint, const PointType& rLowerBound, const PointType& rUpperBound ) const {
    double tolerance = RelativeSnapTolerance(rLowerBound, rUpperBound, SNAPTOL);

    auto p_clipped_mesh = pClipTriangleMesh(rLowerBound, rUpperBound);
    auto& clipped_mesh = *p_clipped_mesh;

    GeometryQuery geometry_query_local(clipped_mesh, false);

    IndexType current_id = 0;
    const IndexType num_triangles = clipped_mesh.NumOfTriangles();
    if( num_triangles == 0 ){ return false; }

    IndexType success_count = 0;
    int inside_count = 0;
    while( success_count < 10 && current_id < num_triangles ){
        // Get direction
        const auto center_triangle = clipped_mesh.Center(current_id);
        Vector3d direction = Math::Subtract( center_triangle, rPoint );

        // Normalize
        double norm_direction = Math::Norm( direction );
        Math::DivideSelf( direction, norm_direction );

        // Construct ray
        Ray_AABB_primitive ray(rPoint, direction);

        // Get vertices of current triangle
        const auto& p1 = clipped_mesh.P1(current_id);
        const auto& p2 = clipped_mesh.P2(current_id);
        const auto& p3 = clipped_mesh.P3(current_id);

        // Make sure target triangle is not parallel and has a significant area.
        const double area = clipped_mesh.Area(current_id);
        if( !ray.is_parallel(p1, p2, p3, 100.0*tolerance) && area >  100*ZEROTOL) {
            auto [is_inside, success] = geometry_query_local.IsInside(ray);
            if( success ){
                ++success_count;
                if( is_inside ){
                    ++inside_count;
                } else {
                    --inside_count;
                }
            }
        }
        ++current_id;
    }
    if( inside_count > 0){
        return true;
    } else {
        return false;
    }

}

IntersectionStatus BRepOperator::GetIntersectionState(
        const PointType& rLowerBound, const PointType& rUpperBound, double Tolerance) const
{

    if( mGeometryQuery.DoIntersect(rLowerBound, rUpperBound, Tolerance) ){
        return IntersectionStatus::Trimmed;
    }

    // Multiple test do not seem to be necessary.
    // Note that for robust analysis, use flood flow.
    // IntersectionStatus status_confirm = (IsInside(center)) ? Inside : Outside;
    // while( status != status_confirm){
    //     status = (IsInside(center)) ? Inside : Outside;
    //     status_confirm = (IsInside(center)) ? Inside : Outside;
    // }

    // Test if center is inside or outside.
    const PointType center = Math::AddAndMult(0.5, rUpperBound, rLowerBound);
    IntersectionStatus status = (IsInside(center)) ? Inside : Outside;

    // If triangle is not intersected, center location will determine if inside or outside.
    return status;
}


Unique<StatusVectorType> BRepOperator::pGetElementClassifications(const Parameters& rParameters) const {
    FloodFill flood_fill(this, rParameters);
    return flood_fill.ClassifyElements();
}

TrimmedDomainPtrType BRepOperator::pGetTrimmedDomain(const PointType& rLowerBound, const PointType& rUpperBound,
        double MinElementVolumeRatio, IndexType MinNumberOfBoundaryTriangles ) const {
    // Instantiate random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> drandon(1, 100);

    // Make copy of lower and upper bound.
    auto lower_bound = rLowerBound;
    auto upper_bound = rUpperBound;

    const double snap_tolerance = RelativeSnapTolerance(lower_bound, upper_bound);
    const auto delta = Math::Subtract(upper_bound, lower_bound);
    const double volume_non_trimmed_domain = delta[0]*delta[1]*delta[2];

    Unique<TrimmedDomain> best_prev_solution = nullptr;
    double best_error = 5e-1;
    bool switch_plane_orientation = false;
    IndexType iteration = 1UL;
    while( iteration < 15){
        auto p_new_mesh = pClipTriangleMesh(lower_bound, upper_bound);
        if( p_new_mesh->NumOfTriangles() > 0) {
            auto p_trimmed_domain = MakeUnique<TrimmedDomain>(std::move(p_new_mesh), lower_bound, upper_bound, this, MinNumberOfBoundaryTriangles, switch_plane_orientation);
            const auto& r_trimmed_domain_mesh = p_trimmed_domain->GetTriangleMesh();

            double volume = MeshUtilities::Volume( r_trimmed_domain_mesh);
            bool volume_ratio = volume / volume_non_trimmed_domain > MinElementVolumeRatio;

            const double epsilon = MeshUtilities::EstimateQuality(r_trimmed_domain_mesh);
            if( epsilon < 1e-5  ){
                if( volume_ratio )
                    return p_trimmed_domain;
                else
                    return nullptr;
            }
            else {
                if( epsilon < best_error && volume_ratio ){
                    best_error = epsilon;
                    best_prev_solution = std::move(p_trimmed_domain);
                }
                if( switch_plane_orientation ){
                    // Perturb AABB slightly.
                    PointType lower_perturbation{drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance};
                    PointType upper_perturbation{drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance, drandon(gen)*snap_tolerance};

                    Math::SubstractSelf( lower_bound, lower_perturbation );
                    Math::AddSelf(upper_bound, upper_perturbation);
                }
            }
        }
        if(switch_plane_orientation ){
            ++iteration;
            switch_plane_orientation = false;
        }
        else {
            switch_plane_orientation = true;
        }
    }

    return best_prev_solution;
}


Unique<TriangleMeshInterface> BRepOperator::pClipTriangleMesh(
        const PointType& rLowerBound, const PointType& rUpperBound ) const {

    const double snap_tolerance = 0.1*RelativeSnapTolerance(rUpperBound, rLowerBound);
    auto p_intersected_triangle_ids = mGeometryQuery.GetIntersectedTriangleIds(rLowerBound, rUpperBound, snap_tolerance);
    auto p_triangle_mesh = MakeUnique<TriangleMesh>();
    p_triangle_mesh->Reserve( 2*p_intersected_triangle_ids->size() );
    p_triangle_mesh->ReserveEdgesOnPlane( p_intersected_triangle_ids->size() );
    for( auto triangle_id : (*p_intersected_triangle_ids) ){
        const auto& P1 = mTriangleMesh.P1(triangle_id);
        const auto& P2 = mTriangleMesh.P2(triangle_id);
        const auto& P3 = mTriangleMesh.P3(triangle_id);
        const auto& normal = mTriangleMesh.Normal(triangle_id);
        auto p_polygon = Clipper::ClipTriangle(P1, P2, P3, normal, rLowerBound, rUpperBound );
        if( p_polygon ){
            p_polygon->AddToTriangleMesh(*p_triangle_mesh.get());
        }
    }
    if( p_triangle_mesh->NumOfTriangles() > 0){
        p_triangle_mesh->Check();
    }
    return p_triangle_mesh;
}


Unique<TriangleMeshInterface> BRepOperator::pClipTriangleMeshUnique(const PointType& rLowerBound, const PointType& rUpperBound ) const {
    const PointType offset{30*ZEROTOL, 30*ZEROTOL, 30*ZEROTOL};
    const auto lower_bound = Math::Add(rLowerBound, offset);
    const auto upper_bound = Math::Add(rUpperBound, offset);
    double snap_tolerance = 1.0*ZEROTOL;

    auto p_intersected_triangle_ids = mGeometryQuery.GetIntersectedTriangleIds(lower_bound, upper_bound, snap_tolerance);
    auto p_triangle_mesh = MakeUnique<TriangleMesh>();
    p_triangle_mesh->Reserve( 2*p_intersected_triangle_ids->size() );
    p_triangle_mesh->ReserveEdgesOnPlane( p_intersected_triangle_ids->size() );

    for( auto triangle_id : (*p_intersected_triangle_ids) ){
        const auto& P1 = mTriangleMesh.P1(triangle_id);
        const auto& P2 = mTriangleMesh.P2(triangle_id);
        const auto& P3 = mTriangleMesh.P3(triangle_id);
        const auto& normal = mTriangleMesh.Normal(triangle_id);
        auto p_polygon = Clipper::ClipTriangle(P1, P2, P3, normal, lower_bound, upper_bound);
        if( p_polygon ){
            p_polygon->AddToTriangleMesh(*p_triangle_mesh.get());
        }
    }
    if( p_triangle_mesh->NumOfTriangles() > 0){
        p_triangle_mesh->Check();
    }
    return p_triangle_mesh;
}

} // End namespace queso

// Winding numbers algorithm: It actually works!!!

// QuESo_INFO << "Done" << std::endl;
// QuESo_INFO << "Tree timer: " << tree_timer << std::endl;
// QuESo_INFO << "ray timer: " << ray_timer << std::endl;
// return MakeUnique<std::vector<bool>>(rr);
// for( int triangle_id = 0; triangle_id < rTriangleMesh.NumOfTriangles(); ++triangle_id ){
//     const auto& p1 = rTriangleMesh.P1(triangle_id);
//     const auto& p2 = rTriangleMesh.P2(triangle_id);
//     const auto& p3 = rTriangleMesh.P3(triangle_id);

//     IndexType count = 0;
//     for( const auto& point : rPoints){
//     //const auto point = rPoints[0];
//         PointType a{p1[0] - point[0], p1[1] - point[1], p1[2] - point[2]};
//         PointType b{p2[0] - point[0], p2[1] - point[1], p2[2] - point[2]};
//         PointType c{p3[0] - point[0], p3[1] - point[1], p3[2] - point[2]};

//         // Compute determinant of:
//         // [a0 b0 c0]
//         // [a1 b1 c1]
//         // [a2 b2 c2]
//         const long double omega =  a[0] * (b[1]*c[2] - b[2]*c[1])
//                                   -b[0] * (a[1]*c[2] - a[2]*c[1])
//                                   +c[0] * (a[1]*b[2] - a[2]*b[1]);


//         const long double anorm = std::sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
//         const long double bnorm = std::sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
//         const long double cnorm = std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

//         long double k = anorm*bnorm*cnorm;

//         k += cnorm * ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
//         k += anorm * ( b[0]*c[0] + b[1]*c[1] + b[2]*c[2]);
//         k += bnorm * ( a[0]*c[0] + a[1]*c[1] + a[2]*c[2]);

//         ret[count] += std::atan2(omega,k);
//         count++;
//     }

// }

// Unique<std::vector<bool>> p_result = MakeUnique<std::vector<bool>>(rPoints.size(), false);
// auto& r_result = *p_result;
// for( int i = 0; i < ret.size(); ++i ){
//     if( ret[i] >= (2.0*My_PI-EPS2) ) // Due to limited precision of double small treshhold (EPS2) is required.
//         r_result[i] = true;
// }

// return std::move(p_result);

// def is_inside_turbo(triangles, X):
// 	# Compute euclidean norm along axis 1
// 	def anorm2(X):
// 		return numpy.sqrt(numpy.sum(X ** 2, axis = 1))



// 	# Compute 3x3 determinant along axis 1
// 	def adet(X, Y, Z):
// 		ret  = numpy.multiply(numpy.multiply(X[:,0], Y[:,1]), Z[:,2])
// 		ret += numpy.multiply(numpy.multiply(Y[:,0], Z[:,1]), X[:,2])
// 		ret += numpy.multiply(numpy.multiply(Z[:,0], X[:,1]), Y[:,2])
// 		ret -= numpy.multiply(numpy.multiply(Z[:,0], Y[:,1]), X[:,2])
// 		ret -= numpy.multiply(numpy.multiply(Y[:,0], X[:,1]), Z[:,2])
// 		ret -= numpy.multiply(numpy.multiply(X[:,0], Z[:,1]), Y[:,2])
// 		return ret



// 	# One generalized winding number per input vertex
// 	ret = numpy.zeros(X.shape[0], dtype = X.dtype)

// 	# Accumulate generalized winding number for each triangle
// 	for U, V, W in triangles:
// 		A, B, C = U - X, V - X, W - X
// 		omega = adet(A, B, C)

// 		a, b, c = anorm2(A), anorm2(B), anorm2(C)
// 		k  = a * b * c
// 		k += c * numpy.sum(numpy.multiply(A, B), axis = 1)
// 		k += a * numpy.sum(numpy.multiply(B, C), axis = 1)
// 		k += b * numpy.sum(numpy.multiply(C, A), axis = 1)

// 		ret += numpy.arctan2(omega, k)

// 	# Job done
// 	return ret >= 2 * numpy.pi