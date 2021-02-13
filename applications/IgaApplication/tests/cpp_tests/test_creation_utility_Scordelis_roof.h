//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_IGA_TEST_CREATION_UTILITY_SCORDELIS_ROOF )
#define  KRATOS_IGA_TEST_CREATION_UTILITY_SCORDELIS_ROOF

#include "geometries/geometry.h"
#include "geometries/nurbs_surface_geometry.h"

namespace Kratos
{
namespace TestCreationUtilityScordelisRoof
{
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node<3> NodeType;
    typedef PointerVector<NodeType> NodeVector;
    typedef Geometry<NodeType> GeometryType;
    typedef NurbsSurfaceGeometry<3, NodeVector> NurbsSurfaceType;

    ///@}
    ///@name Operations
    ///@{

    void AddDisplacementDofs(ModelPart& rModelPart) {
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }
    }

    NurbsSurfaceType GenerateNurbsSurface(ModelPart& rModelPart, SizeType PolynomialDegree) {
        
        SizeType p = PolynomialDegree;
        SizeType q = PolynomialDegree;

        Vector knot_u = ZeroVector(2*(p+1));
        for(IndexType index = knot_u.size()/2; index < knot_u.size(); ++index) {
            knot_u[index] = 1.0;
        }

        Vector knot_v = ZeroVector(2*(q+1));
        for(IndexType index = knot_u.size()/2; index < knot_u.size(); ++index) {
            knot_v[index] = 1.0;
        }
        
        NodeVector points((p+1)*(q+1));

        //define a vector

        if (p == 4)
        {
            points(0) = rModelPart.CreateNewNode(1, -25, -16.069690242163482, 19.151111077974452);
            points(1) = rModelPart.CreateNewNode(2, -12.5, -16.069690242163482, 19.151111077974452);
            points(2) = rModelPart.CreateNewNode(3, 0.0, -16.069690242163482, 19.151111077974452);
            points(3) = rModelPart.CreateNewNode(4, 12.5, -16.069690242163482, 19.151111077974452);
            points(4) = rModelPart.CreateNewNode(5, 25, -16.069690242163482, 19.151111077974452);

            points(5) = rModelPart.CreateNewNode(6, -25, -9.099255856655059, 25);
            points(6) = rModelPart.CreateNewNode(7, -12.5, -9.099255856655059, 25);
            points(7) = rModelPart.CreateNewNode(8, 0.0, -9.099255856655059, 25);
            points(8) = rModelPart.CreateNewNode(9, 12.5, -9.099255856655059, 25);
            points(9) = rModelPart.CreateNewNode(10, 25, -9.099255856655059, 25);

            points(10) = rModelPart.CreateNewNode(11, -25, 0.0, 27.309906636301196);
            points(11) = rModelPart.CreateNewNode(12, -12.5, 0.0, 27.309906636301196);
            points(12) = rModelPart.CreateNewNode(13, 0.0, 0.0, 27.309906636301196);
            points(13) = rModelPart.CreateNewNode(14, 12.5, 0.0, 27.309906636301196);
            points(14) = rModelPart.CreateNewNode(15, 25, 0.0, 27.309906636301196);

            points(15) = rModelPart.CreateNewNode(16, -25, 9.099255856655059, 25);
            points(16) = rModelPart.CreateNewNode(17, -12.5, 9.099255856655059, 25);
            points(17) = rModelPart.CreateNewNode(18, 0.0, 9.099255856655059, 25);
            points(18) = rModelPart.CreateNewNode(19, 12.5, 9.099255856655059, 25);
            points(19) = rModelPart.CreateNewNode(20, 25, 9.099255856655059, 25);

            points(20) = rModelPart.CreateNewNode(21, -25, 16.069690242163482, 19.151111077974452);
            points(21) = rModelPart.CreateNewNode(22, -12.5, 16.069690242163482, 19.151111077974452);
            points(22) = rModelPart.CreateNewNode(23, 0.0, 16.069690242163482, 19.151111077974452);
            points(23) = rModelPart.CreateNewNode(24, 12.5, 16.069690242163482, 19.151111077974452);
            points(24) = rModelPart.CreateNewNode(25, 25, 16.069690242163482, 19.151111077974452);
        }

        if (p == 5)
        {
            points(0) = rModelPart.CreateNewNode(1, -25, -16.069690242163482, 19.151111077974452);
            points(1) = rModelPart.CreateNewNode(2, -15, -16.069690242163482, 19.151111077974452);
            points(2) = rModelPart.CreateNewNode(3, -5, -16.069690242163482, 19.151111077974452);
            points(3) = rModelPart.CreateNewNode(4, 5, -16.069690242163482, 19.151111077974452);
            points(4) = rModelPart.CreateNewNode(5, 15, -16.069690242163482, 19.151111077974452);
            points(5) = rModelPart.CreateNewNode(6, 25, -16.069690242163482, 19.151111077974452);

            points(6) = rModelPart.CreateNewNode(7, -25, -10.637273878912895, 23.709449644779440);
            points(7) = rModelPart.CreateNewNode(8, -15, -10.637273878912895, 23.709449644779440);
            points(8) = rModelPart.CreateNewNode(9, -5, -10.637273878912895, 23.709449644779440);
            points(9) = rModelPart.CreateNewNode(10, 5, -10.637273878912895, 23.709449644779440);
            points(10) = rModelPart.CreateNewNode(11, 15, -10.637273878912895, 23.709449644779440);
            points(11) = rModelPart.CreateNewNode(12, 25, -10.637273878912895, 23.709449644779440);
            
            points(12) = rModelPart.CreateNewNode(13, -25, -3.738760296802554, 26.360797461092400);
            points(13) = rModelPart.CreateNewNode(14, -15, -3.738760296802554, 26.360797461092400);
            points(14) = rModelPart.CreateNewNode(15, -5, -3.738760296802554, 26.360797461092400);
            points(15) = rModelPart.CreateNewNode(16, 5, -3.738760296802554, 26.360797461092400);
            points(16) = rModelPart.CreateNewNode(17, 15, -3.738760296802554, 26.360797461092400);
            points(17) = rModelPart.CreateNewNode(18, 25, -3.738760296802554, 26.360797461092400);

            points(18) = rModelPart.CreateNewNode(19, -25, 3.7387602968025549, 26.360797461092400);
            points(19) = rModelPart.CreateNewNode(20, -15, 3.738760296802554, 26.360797461092400);
            points(20) = rModelPart.CreateNewNode(21, -5, 3.738760296802554, 26.360797461092400);
            points(21) = rModelPart.CreateNewNode(22, 5, 3.738760296802554, 26.360797461092400);
            points(22) = rModelPart.CreateNewNode(23, 15, 3.738760296802554, 26.360797461092400);
            points(23) = rModelPart.CreateNewNode(24, 25, 3.738760296802554, 26.360797461092400);

            points(24) = rModelPart.CreateNewNode(25, -25, 10.6372738789128952, 23.709449644779440);
            points(25) = rModelPart.CreateNewNode(26, -15, 10.637273878912895, 23.709449644779440);
            points(26) = rModelPart.CreateNewNode(27, -5, 10.637273878912895, 23.709449644779440);
            points(27) = rModelPart.CreateNewNode(28, 5, 10.637273878912895, 23.709449644779440);
            points(28) = rModelPart.CreateNewNode(29, 15, 10.637273878912895, 23.709449644779440);
            points(29) = rModelPart.CreateNewNode(30, 25, 10.637273878912895, 23.709449644779440);

            points(30) = rModelPart.CreateNewNode(31, -25, 16.069690242163482, 19.151111077974452);
            points(31) = rModelPart.CreateNewNode(32, -15, 16.069690242163482, 19.151111077974452);
            points(32) = rModelPart.CreateNewNode(33, -5, 16.069690242163482, 19.151111077974452);
            points(33) = rModelPart.CreateNewNode(34, 5, 16.069690242163482, 19.151111077974452);
            points(34) = rModelPart.CreateNewNode(35, 15, 16.069690242163482, 19.151111077974452);
            points(35) = rModelPart.CreateNewNode(36, 25, 16.069690242163482, 19.151111077974452);
        }

        Vector weight = ZeroVector((p+1)*(q+1));

        if (p == 4)
        {
            weight[0] = 1.000000000000000;
            weight[1] = 1.000000000000000;
            weight[2] = 1.000000000000000;
            weight[3] = 1.000000000000000;
            weight[4] = 1.000000000000000;

            weight[5] = 0.883022221559489;
            weight[6] = 0.883022221559489;
            weight[7] = 0.883022221559489;
            weight[8] = 0.883022221559489;
            weight[9] = 0.883022221559489;

            weight[10] = 0.844029628745985;
            weight[11] = 0.844029628745985;
            weight[12] = 0.844029628745985;
            weight[13] = 0.844029628745985;
            weight[14] = 0.844029628745985;

            weight[15] = 0.883022221559489;
            weight[16] = 0.883022221559489;
            weight[17] = 0.883022221559489;
            weight[18] = 0.883022221559489;
            weight[19] = 0.883022221559489;

            weight[20] = 1.000000000000000;
            weight[21] = 1.000000000000000;
            weight[22] = 1.000000000000000;
            weight[23] = 1.000000000000000;
            weight[24] = 1.000000000000000;
        }
        if (p == 5)
        {
            weight[0] = 1.000000000000000;
            weight[1] = 1.000000000000000;
            weight[2] = 1.000000000000000;
            weight[3] = 1.000000000000000;
            weight[4] = 1.000000000000000;
            weight[5] = 1.000000000000000;

            weight[6] = 0.906417777247591;
            weight[7] = 0.906417777247591;
            weight[8] = 0.906417777247591;
            weight[9] = 0.906417777247591;
            weight[10] = 0.906417777247591;
            weight[11] = 0.906417777247591;

            weight[12] = 0.859626665871387;
            weight[13] = 0.859626665871387;
            weight[14] = 0.859626665871387;
            weight[15] = 0.859626665871387;
            weight[16] = 0.859626665871387;
            weight[17] = 0.859626665871387;

            weight[18] = 0.859626665871387;
            weight[19] = 0.859626665871387;
            weight[20] = 0.859626665871387;
            weight[21] = 0.859626665871387;
            weight[22] = 0.859626665871387;
            weight[23] = 0.859626665871387;

            weight[24] = 0.906417777247591;
            weight[25] = 0.906417777247591;
            weight[26] = 0.906417777247591;
            weight[27] = 0.906417777247591;
            weight[28] = 0.906417777247591;
            weight[29] = 0.906417777247591;

            weight[30] = 1.000000000000000;
            weight[31] = 1.000000000000000;
            weight[32] = 1.000000000000000;
            weight[33] = 1.000000000000000;
            weight[34] = 1.000000000000000;
            weight[35] = 1.000000000000000;
        }

        return NurbsSurfaceType(
            points, p, q, knot_u, knot_v, weight);
    }

    typename Geometry<NodeType>::Pointer GetQuadraturePointGeometry(
        ModelPart& rModelPart, SizeType PolynomialDegree, IntegrationPoint<3> IntegrationPoint)
    {
        typename GeometryType::IntegrationPointsArrayType integration_points(1);
        integration_points[0] = IntegrationPoint;
        typename GeometryType::GeometriesArrayType result_geometries;
        
        GenerateNurbsSurface(rModelPart,PolynomialDegree).CreateQuadraturePointGeometries(
                result_geometries, 3, integration_points);
        
        return result_geometries(0);
    }

    ///@}
}; /// namespace TestCreationUtilityScordelisRoof
} /// namespace Kratos.

#endif /// KRATOS_IGA_TEST_CREATION_UTILITY_SCORDELIS_ROOF
