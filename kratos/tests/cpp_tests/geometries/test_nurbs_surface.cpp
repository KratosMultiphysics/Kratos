//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//                   Andreas Apostolatos
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_surface_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3> NodeType;

    /// Factory functions
    NurbsSurfaceGeometry<3, Point> GenerateReferencePointSurface()
    {
        NurbsSurfaceGeometry<3, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(-10.0, -5.0, -1.0)));
        points.push_back(Point::Pointer(new Point(-12.0, 3.0, 3.0)));
        points.push_back(Point::Pointer(new Point(-9.0, 11.0, -0.0701928417)));
        points.push_back(Point::Pointer(new Point(-5.0, -3.0, 1.0)));
        points.push_back(Point::Pointer(new Point(-6.0, 4.0, -2.0)));
        points.push_back(Point::Pointer(new Point(-5.0, 7.0, 0.9298071583)));
        points.push_back(Point::Pointer(new Point(0.0, -4.0, -1.0)));
        points.push_back(Point::Pointer(new Point(1.0, 6.0, 5.0)));
        points.push_back(Point::Pointer(new Point(0.0, 13.0, -0.2350184214)));
        points.push_back(Point::Pointer(new Point(4.0, -2.0, 0.0)));
        points.push_back(Point::Pointer(new Point(5.0, 4.0, -1.0)));
        points.push_back(Point::Pointer(new Point(5.0, 11.0, 0.7649815786)));

        Vector knot_vector_u = ZeroVector(5);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 7.5;
        knot_vector_u[3] = 15.0;
        knot_vector_u[4] = 15.0;

        Vector knot_vector_v = ZeroVector(3);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 10.0;
        knot_vector_v[2] = 20.0;

        int p = 2;
        int q = 1;

        Vector weights = ZeroVector(12);
        weights[0] = 1.0;
        weights[1] = 1.0;
        weights[2] = 1.0;
        weights[3] = 1.0;
        weights[4] = 2.5;
        weights[5] = 1.0;
        weights[6] = 1.0;
        weights[7] = 1.0;
        weights[8] = 1.0;
        weights[9] = 1.0;
        weights[10] = 1.0;
        weights[11] = 1.0;

        return NurbsSurfaceGeometry<3, Point>(
                points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    NurbsSurfaceGeometry<3, Point> GenerateReferencePieceOfCylinderNurbsSurface()
    {
        NurbsSurfaceGeometry<3, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(0, 10, 0)));
        points.push_back(Point::Pointer(new Point(6.6817863791929888, 10, 0)));
        points.push_back(Point::Pointer(new Point(9.2387953251128678, 3.8268343236508979, 0)));
        points.push_back(Point::Pointer(new Point(11.795804271032745, -2.3463313526982033, 0)));
        points.push_back(Point::Pointer(new Point(7.0710678118654755, -7.0710678118654755, 0)));
        points.push_back(Point::Pointer(new Point(0, 10, 10)));
        points.push_back(Point::Pointer(new Point(6.6817863791929888, 10, 10)));
        points.push_back(Point::Pointer(new Point(9.2387953251128678, 3.8268343236508979, 10)));
        points.push_back(Point::Pointer(new Point(11.795804271032745, -2.3463313526982033, 10)));
        points.push_back(Point::Pointer(new Point(7.0710678118654755, -7.0710678118654755)));

        Vector knot_vector_u = ZeroVector(6);
        knot_vector_u[0] = 0.0;
        knot_vector_u[1] = 0.0;
        knot_vector_u[2] = 11.780972450961723;
        knot_vector_u[3] = 11.780972450961723;
        knot_vector_u[4] = 23.561944901923447;
        knot_vector_u[5] = 23.561944901923447;

        Vector knot_vector_v = ZeroVector(2);
        knot_vector_v[0] = 0.0;
        knot_vector_v[1] = 10.0;

        int p = 2;
        int q = 1;

        Vector weights = ZeroVector(10);
        weights[0] = 1.0;
        weights[1] = 0.83146961230254524;
        weights[2] = 1.0;
        weights[3] = 0.83146961230254524;
        weights[4] = 1.0;
        weights[5] = 1.0;
        weights[6] = 0.83146961230254524;
        weights[7] = 1.0;
        weights[8] = 0.83146961230254524;
        weights[9] = 1.0;

        return NurbsSurfaceGeometry<3, Point>(
            points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    NurbsSurfaceGeometry<3, NodeType> GenerateReferenceNodeSurface() {
        Geometry<NodeType>::PointsArrayType points;

        points.push_back(NodeType::Pointer(new NodeType(1, 0, 5, 0)));
        points.push_back(NodeType::Pointer(new NodeType(2, 5, 5, 0)));
        points.push_back(NodeType::Pointer(new NodeType(3, 10, 5, -4)));
        points.push_back(NodeType::Pointer(new NodeType(4, 0, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(5, 5, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(6, 10, 0, -4)));

        Vector knot_u = ZeroVector(4);
        knot_u[0] = 0.0;
        knot_u[1] = 0.0;
        knot_u[2] = 10.0;
        knot_u[3] = 10.0;
        Vector knot_v = ZeroVector(2); 
        knot_v[0] = 0.0;
        knot_v[1] = 5.0;

        int p = 2;
        int q = 1;

        return NurbsSurfaceGeometry<3, NodeType>(
            points, p, q, knot_u, knot_v);
    }

    NurbsSurfaceGeometry<3, Point> GenerateReferenceQuarterSphereGeometry()
    {
        NurbsSurfaceGeometry<3, Point>::PointsArrayType points;

        points.push_back(Point::Pointer(new Point(000000000000000, -7.500000000000000e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -7.500000000000000e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(1.335563126175092e-02, -7.435787115897656e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -7.127196936176004e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(4.491844452241197e-02, -6.261347715628729e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(6.261347715628729e-02, -4.491844452241197e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(7.127196936176004e-02, -2.721342545528786e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -1.335563126175091e-02, 000000000000000)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, -4.344437182340776e-03, 000000000000000)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -7.500000000000000e-02, 4.344437182340777e-03)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -7.500000000000000e-02, 4.344437182340778e-03)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -7.435787115897656e-02, 4.307241336836834e-03)));
        points.push_back(Point::Pointer(new Point(2.721342545528787e-02, -7.127196936176004e-02, 4.128487916718439e-03)));
        points.push_back(Point::Pointer(new Point(4.491844452241198e-02, -6.261347715628728e-02, 3.626937576978924e-03)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -4.491844452241198e-02, 2.601938140747707e-03)));
        points.push_back(Point::Pointer(new Point(7.127196936176004e-02, -2.721342545528786e-02, 1.576360232090821e-03)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -1.335563126175091e-02, 7.736360139624470e-04)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, -4.344437182340775e-03, 2.516551257507342e-04)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -7.435787115897656e-02, 1.335563126175092e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -7.435787115897657e-02, 1.335563126175092e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -7.372124004393278e-02, 1.324128411477432e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -7.066175886731038e-02, 1.269176189459301e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241197e-02, -6.207739822936977e-02, 1.114990017220583e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -4.453386547278862e-02, 7.998855758569994e-03)));
        points.push_back(Point::Pointer(new Point(7.127196936176004e-02, -2.698043178398277e-02, 4.846033009999609e-03)));
        points.push_back(Point::Pointer(new Point(7.435787115897655e-02, -1.324128411477432e-02, 2.378305151998111e-03)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, -4.307241336836834e-03, 7.736360139624470e-04)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -7.127196936176004e-02, 2.721342545528786e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -7.127196936176006e-02, 2.721342545528787e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -7.066175886731038e-02, 2.698043178398277e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -6.772924822271549e-02, 2.586072567037090e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241198e-02, -5.950114434021560e-02, 2.271902924118659e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628729e-02, -4.268568002372353e-02, 1.629846322104189e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176003e-02, -2.586072567037090e-02, 9.874273666806793e-03)));
        points.push_back(Point::Pointer(new Point(7.435787115897655e-02, -1.269176189459300e-02, 4.846033009999608e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -4.128487916718439e-03, 1.576360232090820e-03)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -6.261347715628729e-02, 4.491844452241197e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -6.261347715628729e-02, 4.491844452241196e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -6.207739822936976e-02, 4.453386547278861e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -5.950114434021559e-02, 4.268568002372353e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241198e-02, -5.227263362134548e-02, 3.750000000000001e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -3.750000000000001e-02, 2.690222211084004e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176002e-02, -2.271902924118658e-02, 1.629846322104189e-02)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -1.114990017220582e-02, 7.998855758569987e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -3.626937576978923e-03, 2.601938140747705e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -4.491844452241197e-02, 6.261347715628729e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -4.491844452241198e-02, 6.261347715628728e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -4.453386547278862e-02, 6.207739822936977e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -4.268568002372353e-02, 5.950114434021559e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241198e-02, -3.750000000000001e-02, 5.227263362134548e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -2.690222211084004e-02, 3.750000000000001e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176003e-02, -1.629846322104189e-02, 2.271902924118658e-02)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -7.998855758569991e-03, 1.114990017220582e-02)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, -2.601938140747705e-03, 3.626937576978923e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -2.721342545528786e-02, 7.127196936176004e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -2.721342545528787e-02, 7.127196936176004e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -2.698043178398278e-02, 7.066175886731040e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -2.586072567037090e-02, 6.772924822271549e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241198e-02, -2.271902924118659e-02, 5.950114434021560e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -1.629846322104189e-02, 4.268568002372353e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176003e-02, -9.874273666806793e-03, 2.586072567037090e-02)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -4.846033009999609e-03, 1.269176189459301e-02)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, -1.576360232090820e-03, 4.128487916718438e-03)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -1.335563126175091e-02, 7.435787115897656e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -1.335563126175091e-02, 7.435787115897657e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175091e-02, -1.324128411477432e-02, 7.372124004393280e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -1.269176189459301e-02, 7.066175886731038e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241197e-02, -1.114990017220583e-02, 6.207739822936977e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -7.998855758569991e-03, 4.453386547278862e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176004e-02, -4.846033009999609e-03, 2.698043178398277e-02)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -2.378305151998111e-03, 1.324128411477432e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -7.736360139624468e-04, 4.307241336836834e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, -4.344437182340776e-03, 7.499999999999998e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, -4.344437182340776e-03, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175092e-02, -4.307241336836833e-03, 7.435787115897656e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, -4.128487916718437e-03, 7.127196936176003e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241196e-02, -3.626937576978923e-03, 6.261347715628728e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628728e-02, -2.601938140747705e-03, 4.491844452241196e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176003e-02, -1.576360232090820e-03, 2.721342545528786e-02)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, -7.736360139624468e-04, 1.335563126175091e-02)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, -2.516551257507341e-04, 4.344437182340776e-03)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, 000000000000000, 000000000000000)));

        points.push_back(Point::Pointer(new Point(000000000000000, 000000000000000, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(4.344437182340777e-03, 000000000000000, 7.500000000000000e-02)));
        points.push_back(Point::Pointer(new Point(1.335563126175092e-02, 000000000000000, 7.435787115897656e-02)));
        points.push_back(Point::Pointer(new Point(2.721342545528786e-02, 000000000000000, 7.127196936176004e-02)));
        points.push_back(Point::Pointer(new Point(4.491844452241197e-02, 000000000000000, 6.261347715628729e-02)));
        points.push_back(Point::Pointer(new Point(6.261347715628729e-02, 000000000000000, 4.491844452241197e-02)));
        points.push_back(Point::Pointer(new Point(7.127196936176004e-02, 000000000000000, 2.721342545528786e-02)));
        points.push_back(Point::Pointer(new Point(7.435787115897656e-02, 000000000000000, 1.335563126175091e-02)));
        points.push_back(Point::Pointer(new Point(7.499999999999998e-02, 000000000000000, 4.344437182340776e-03)));
        points.push_back(Point::Pointer(new Point(7.500000000000000e-02, 000000000000000, 000000000000000)));

        Vector knot_vector_u = ZeroVector(14);
        knot_vector_u[0] = .0;
        knot_vector_u[1] = .0;
        knot_vector_u[2] = .0;
        knot_vector_u[3] = .0;
        knot_vector_u[4] = .0;
        knot_vector_u[5] = .2;
        knot_vector_u[6] = .4;
        knot_vector_u[7] = .6;
        knot_vector_u[8] = .8;
        knot_vector_u[9] = 1.0;
        knot_vector_u[10] = 1.0;
        knot_vector_u[11] = 1.0;
        knot_vector_u[12] = 1.0;
        knot_vector_u[13] = 1.0;

        Vector knot_vector_v = ZeroVector(14);
        knot_vector_v[0] = .0;
        knot_vector_v[1] = .0;
        knot_vector_v[2] = .0;
        knot_vector_v[3] = .0;
        knot_vector_v[4] = .0;
        knot_vector_v[5] = .2;
        knot_vector_v[6] = .4;
        knot_vector_v[7] = .6;
        knot_vector_v[8] = .8;
        knot_vector_v[9] = 1.0;
        knot_vector_v[10] = 1.0;
        knot_vector_v[11] = 1.0;
        knot_vector_v[12] = 1.0;
        knot_vector_v[13] = 1.0;

        int p = 5;
        int q = 5;

        Vector weights = ZeroVector(100);
        weights[0] = 1.0;
        weights[1] = 9.765685424949238e-01;
        weights[2] = 9.343919189857868e-01;
        weights[3] = 8.851858582251269e-01;
        weights[4] = 8.476955262170048e-01;
        weights[5] = 8.476955262170048e-01;
        weights[6] = 8.851858582251269e-01;
        weights[7] = 9.343919189857868e-01;
        weights[8] = 9.765685424949239e-01;
        weights[9] = 1.0;

        weights[10] = 9.765685424949238e-01;
        weights[11] = 9.536861181906597e-01;
        weights[12] = 9.124977544429849e-01;
        weights[13] = 8.644446634040304e-01;
        weights[14] = 8.278327845172080e-01;
        weights[15] = 8.278327845172080e-01;
        weights[16] = 8.644446634040304e-01;
        weights[17] = 9.124977544429849e-01;
        weights[18] = 9.536861181906600e-01;
        weights[19] = 9.765685424949238e-01;

        weights[20] = 9.343919189857868e-01;
        weights[21] = 9.124977544429846e-01;
        weights[22] = 8.730882582659412e-01;
        weights[23] = 8.271105127260567e-01;
        weights[24] = 7.920798494575735e-01;
        weights[25] = 7.920798494575735e-01;
        weights[26] = 8.271105127260567e-01;
        weights[27] = 8.730882582659412e-01;
        weights[28] = 9.124977544429846e-01;
        weights[29] = 9.343919189857868e-01;

        weights[30] = 8.851858582251269e-01;
        weights[31] = 8.644446634040303e-01;
        weights[32] = 8.271105127260568e-01;
        weights[33] = 7.835540036017543e-01;
        weights[34] = 7.503680918879998e-01;
        weights[35] = 7.503680918879998e-01;
        weights[36] = 7.835540036017543e-01;
        weights[37] = 8.271105127260568e-01;
        weights[38] = 8.644446634040303e-01;
        weights[39] = 8.851858582251269e-01;

        weights[40] = 8.476955262170048e-01;
        weights[41] = 8.278327845172080e-01;
        weights[42] = 7.920798494575736e-01;
        weights[43] = 7.503680918880000e-01;
        weights[44] = 7.185877051683247e-01;
        weights[45] = 7.185877051683247e-01;
        weights[46] = 7.503680918880000e-01;
        weights[47] = 7.920798494575736e-01;
        weights[48] = 8.278327845172079e-01;
        weights[49] = 8.476955262170048e-01;

        weights[50] = 8.476955262170048e-01;
        weights[51] = 8.278327845172080e-01;
        weights[52] = 7.920798494575735e-01;
        weights[53] = 7.503680918880000e-01;
        weights[54] = 7.185877051683246e-01;
        weights[55] = 7.185877051683246e-01;
        weights[56] = 7.503680918880000e-01;
        weights[57] = 7.920798494575735e-01;
        weights[58] = 8.278327845172080e-01;
        weights[59] = 8.476955262170048e-01;

        weights[60] = 8.851858582251269e-01;
        weights[61] = 8.644446634040304e-01;
        weights[62] = 8.271105127260567e-01;
        weights[63] = 7.835540036017543e-01;
        weights[64] = 7.503680918879998e-01;
        weights[65] = 7.503680918879998e-01;
        weights[66] = 7.835540036017543e-01;
        weights[67] = 8.271105127260567e-01;
        weights[68] = 8.644446634040304e-01;
        weights[69] = 8.851858582251269e-01;

        weights[70] = 9.343919189857868e-01;
        weights[71] = 9.124977544429846e-01;
        weights[72] = 8.730882582659411e-01;
        weights[73] = 8.271105127260567e-01;
        weights[74] = 7.920798494575735e-01;
        weights[75] = 7.920798494575735e-01;
        weights[76] = 8.271105127260567e-01;
        weights[77] = 8.730882582659411e-01;
        weights[78] = 9.124977544429846e-01;
        weights[79] = 9.343919189857868e-01;

        weights[80] = 9.765685424949239e-01;
        weights[81] = 9.536861181906598e-01;
        weights[82] = 9.124977544429848e-01;
        weights[83] = 8.644446634040305e-01;
        weights[84] = 8.278327845172080e-01;
        weights[85] = 8.278327845172080e-01;
        weights[86] = 8.644446634040305e-01;
        weights[87] = 9.124977544429848e-01;
        weights[88] = 9.536861181906598e-01;
        weights[89] = 9.765685424949239e-01;

        weights[90] = 1.0;
        weights[91] = 9.765685424949238e-01;
        weights[92] = 9.343919189857868e-01;
        weights[93] = 8.851858582251269e-01;
        weights[94] = 8.476955262170048e-01;
        weights[95] = 8.476955262170048e-01;
        weights[96] = 8.851858582251269e-01;
        weights[97] = 9.343919189857868e-01;
        weights[98] = 9.765685424949239e-01;
        weights[99] = 1.0;

        return NurbsSurfaceGeometry<3, Point>(
            points, p, q, knot_vector_u, knot_vector_v, weights);
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfacePoint, KratosCoreNurbsGeometriesFastSuite) {
        auto surface = GenerateReferencePointSurface();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), true);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 2);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 1);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 5);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 3);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 4);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 3);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 12);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.0;
        parameter[1] = 0.0;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCylinderSurface, KratosCoreNurbsGeometriesFastSuite) {
        auto surface = GenerateReferencePieceOfCylinderNurbsSurface();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), true);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 2);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 1);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 6);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 2);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 5);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 2);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 10);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 10.0;
        parameter[1] = 3.5;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
        double length = sqrt(result[0] * result[0] + result[1] * result[1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], parameter[1], TOLERANCE);

        auto derivatives = surface.GlobalDerivatives(parameter, 3);
        array_1d<double, 3> cross(0.0);
        array_1d<double, 3> colinear_vector(0.0);
        derivatives[0][2] = 0.0;
        MathUtils<double>::CrossProduct(cross, derivatives[1], derivatives[2]);
        MathUtils<double>::CrossProduct(colinear_vector, cross, derivatives[0]);
        KRATOS_CHECK_NEAR(norm_2(colinear_vector), 0, TOLERANCE);

        parameter[0] = 6.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        length = sqrt(result[0] * result[0] + result[1] * result[1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], parameter[1], TOLERANCE);

        parameter[0] = 0.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        length = sqrt(result[0] * result[0] + result[1] * result[1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], parameter[1], TOLERANCE);

        parameter[0] = 0;
        parameter[1] = 1.0;
        auto derivatives_2 = surface.GlobalDerivatives(parameter, 3);
        length = sqrt(derivatives_2[0][0] * derivatives_2[0][0]
            + derivatives_2[0][1] * derivatives_2[0][1]);
        KRATOS_CHECK_NEAR(length, 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[0][2], parameter[1], TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][0]/norm_2(derivatives_2[1]), 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][2], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][2], 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsSurfaceNode, KratosCoreNurbsGeometriesFastSuite) {
        auto surface = GenerateReferenceNodeSurface();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), false);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 2);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 1);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 4);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 2);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 3);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 2);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 6);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 10.0;
        parameter[1] = 3.5;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 10.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 1.5, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], -4.0, TOLERANCE);

        parameter[0] = 6.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 6.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 4.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], - 1.44, TOLERANCE);

        parameter[0] = 0.0;
        parameter[1] = 1.0;

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 4.0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsQuarterSphereSurface, KratosCoreNurbsGeometriesFastSuite) {
        auto surface = GenerateReferenceQuarterSphereGeometry();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(surface.Dimension(), 2);
        KRATOS_CHECK_EQUAL(surface.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(surface.LocalSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(surface.IsRational(), true);

        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeU(), 5);
        KRATOS_CHECK_EQUAL(surface.PolynomialDegreeV(), 5);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsU(), 14);
        KRATOS_CHECK_EQUAL(surface.NumberOfKnotsV(), 14);

        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsU(), 10);
        KRATOS_CHECK_EQUAL(surface.NumberOfControlPointsV(), 10);
        KRATOS_CHECK_EQUAL(surface.PointsNumber(), 100);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = .25;
        parameter[1] = .75;
        array_1d<double, 3> result(0.0);

        surface.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 0.027607103217140, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], -0.025668761597520, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 0.064837971359442, TOLERANCE);

        auto derivatives = surface.GlobalDerivatives(parameter, 6);

        // Compare the position vectors and the gradients up to 5th order
        double positionVct[3] = {0.027607103217140, -0.025668761597520, 0.064837971359442};
        double gradient1[3] = {0.110787255345493, 0.016144510322666, -0.040780202579556};
        double gradient2[3] = {-0.000000000000000, 0.103008693927056, 0.040780202579556};
        double gradient3[3] = {-0.033227653950941, 0.070099959163148, -0.177068890809183};
        double gradient4[3] = {-0.000000000000000, -0.064787890764092, -0.025648935146511};
        double gradient5[3] = {-0.000000000000000, 0.030894683915335, -0.177068890809183};
        double gradient6[3] = {-0.470230208134118, 0.005450671128399, -0.013768114880420};
        double gradient7[3] = {-0.000000000000003, -0.281311009504776, -0.111368463360849};
        double gradient8[3] = {-0.000000000000001, -0.019431383220095, 0.111368463360851};
        double gradient9[3] = {-0.000000000000002, -0.437214546329255, 0.013768114880415};
        double gradient10[3] = {-0.356492597666529, -0.546390143805925, 1.380153396204004};
        double gradient11[3] = {-0.000000000000008, -0.021873533393088, -0.008659532403435};
        double gradient12[3] = {0.000000000000013, -0.084371662130872, 0.483565284894256};
        double gradient13[3] = {-0.000000000000012, 0.274988519785633, -0.008659532403469};
        double gradient14[3] = {-0.000000000000040, 0.331462646725662, 1.380153396204038};

        for (IndexType i = 0; i < 3; ++i) {
            KRATOS_CHECK_NEAR(derivatives[0][i], positionVct[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][i], gradient1[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][i], gradient2[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][i], gradient3[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[4][i], gradient4[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[5][i], gradient5[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[6][i], gradient6[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[7][i], gradient7[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[8][i], gradient8[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[9][i], gradient9[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[10][i], gradient10[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[11][i], gradient11[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[12][i], gradient12[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[13][i], gradient13[i], TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[14][i], gradient14[i], TOLERANCE);
        }
    }
} // namespace Testing.
} // namespace Kratos.
