//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

namespace
{

class TestLaplacianElement1D : public Element
{
public:

    using IndexType = Element::IndexType;

    using DofsArrayType = Element::DofsArrayType;

    using DofsVectorType = Element::DofsVectorType;

    using EquationIdVectorType = Element::EquationIdVectorType;

    TestLaplacianElement1D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(
        MatrixType &rLeftHandSideMatrix,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rLeftHandSideMatrix.size1() != n_nodes || rLeftHandSideMatrix.size2() != n_nodes) {
            rLeftHandSideMatrix.resize(n_nodes, n_nodes, false);
        }
        const double aux = 1.0 / r_geom.Length();
        rLeftHandSideMatrix(0,0) = aux;
        rLeftHandSideMatrix(0,1) = -aux;
        rLeftHandSideMatrix(1,0) = -aux;
        rLeftHandSideMatrix(1,1) = aux;
    }

    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rRightHandSideVector.size() != n_nodes) {
            rRightHandSideVector.resize(n_nodes);
        }
        const double f = 1.0;
        const double h = r_geom.Length();
        const double aux_f = f * h / 2.0;
        const double aux_phi =  (r_geom[0].GetSolutionStepValue(DISTANCE) - r_geom[1].GetSolutionStepValue(DISTANCE)) / h;
        rRightHandSideVector(0) = aux_f - aux_phi;
        rRightHandSideVector(1) = aux_f + aux_phi;
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rResult.size() != n_nodes) {
            rResult.resize(n_nodes);
        }
        for (unsigned int i = 0; i < n_nodes; ++i) {
            rResult[i] = r_geom[i].GetDof(DISTANCE).EquationId();
        }
    }

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        if (rElementalDofList.size() != n_nodes) {
            rElementalDofList.resize(n_nodes);
        }
        for (unsigned int i = 0; i < n_nodes; i++) {
            rElementalDofList[i] = r_geom[i].pGetDof(DISTANCE);
        }
    }

};

class TestElasticityElement2D : public Element
{
public:

    using IndexType = Element::IndexType;

    using DofsArrayType = Element::DofsArrayType;

    using DofsVectorType = Element::DofsVectorType;

    using EquationIdVectorType = Element::EquationIdVectorType;

    TestElasticityElement2D(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {}

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
        CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(
        MatrixType &rLeftHandSideMatrix,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto& r_geom = GetGeometry();
        const std::size_t local_size = 2 * r_geom.PointsNumber();
        if (rLeftHandSideMatrix.size1() != local_size || rLeftHandSideMatrix.size2() != local_size) {
            rLeftHandSideMatrix.resize(local_size, local_size, false);
        }

        const double C = E / (1.0 - std::pow(nu,2));
        const double a = r_geom[1].X0() - r_geom[0].X0();
        const double b = r_geom[3].Y0() - r_geom[0].Y0();
        const double k1 = C *((b/(3.0*a))*(1+nu) + (a/(3.0*b))*(1-nu));
        const double k2 = C *((b/(6.0*a))*(1+nu) - (a/(6.0*b))*(1-nu));
        const double k3 = C *((b/(6.0*a))*(1+nu) + (a/(6.0*b))*(1-nu));
        const double k4 = C *((b/(3.0*a))*(1+nu) - (a/(3.0*b))*(1-nu));
        const double k5 = C *((b/(4.0*a))*(1-nu) + (a/(4.0*b))*(1-nu));
        const double k6 = C *((b/(4.0*a))*(1-nu) - (a/(4.0*b))*(1-nu));

        rLeftHandSideMatrix(0,0) = k1;
        rLeftHandSideMatrix(1,1) = k5;
        rLeftHandSideMatrix(2,2) = k1;
        rLeftHandSideMatrix(3,3) = k5;
        rLeftHandSideMatrix(4,4) = k1;
        rLeftHandSideMatrix(5,5) = k5;
        rLeftHandSideMatrix(6,6) = k1;
        rLeftHandSideMatrix(7,7) = k5;

        rLeftHandSideMatrix(0,2) = k2;
        rLeftHandSideMatrix(0,4) = -k3;
        rLeftHandSideMatrix(0,6) = -k4;
        rLeftHandSideMatrix(2,0) = rLeftHandSideMatrix(0,2);
        rLeftHandSideMatrix(4,0) = rLeftHandSideMatrix(0,4);
        rLeftHandSideMatrix(6,0) = rLeftHandSideMatrix(0,6);

        rLeftHandSideMatrix(1,3) = k6;
        rLeftHandSideMatrix(1,5) = -k5;
        rLeftHandSideMatrix(1,7) = -k6;
        rLeftHandSideMatrix(3,1) = rLeftHandSideMatrix(1,3);
        rLeftHandSideMatrix(5,1) = rLeftHandSideMatrix(1,5);
        rLeftHandSideMatrix(7,1) = rLeftHandSideMatrix(1,7);

        rLeftHandSideMatrix(2,4) = -k4;
        rLeftHandSideMatrix(2,6) = -k3;
        rLeftHandSideMatrix(4,2) = rLeftHandSideMatrix(2,4);
        rLeftHandSideMatrix(6,2) = rLeftHandSideMatrix(2,6);

        rLeftHandSideMatrix(3,5) = -k6;
        rLeftHandSideMatrix(3,7) = -k5;
        rLeftHandSideMatrix(5,3) = rLeftHandSideMatrix(3,5);
        rLeftHandSideMatrix(7,3) = rLeftHandSideMatrix(3,7);

        rLeftHandSideMatrix(4,6) = k2;
        rLeftHandSideMatrix(6,4) = rLeftHandSideMatrix(4,6);

        rLeftHandSideMatrix(5,7) = k6;
        rLeftHandSideMatrix(7,5) = rLeftHandSideMatrix(5,7);
    }

    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override
    {
        const auto &r_geom = GetGeometry();
        const std::size_t local_size = 2 * r_geom.PointsNumber();
        if (rRightHandSideVector.size() != local_size) {
            rRightHandSideVector.resize(local_size);
        }
        const double f_x = 0.0;
        const double f_y = -1.0;
        const double a = r_geom[1].X0() - r_geom[0].X0();
        const double b = r_geom[3].Y0() - r_geom[0].Y0();
        const double aux = a * b / 4.0;

        const double C = E / (1.0 - std::pow(nu,2));
        const double k1 = C *((b/(3.0*a))*(1+nu) + (a/(3.0*b))*(1-nu));
        const double k2 = C *((b/(6.0*a))*(1+nu) - (a/(6.0*b))*(1-nu));
        const double k3 = C *((b/(6.0*a))*(1+nu) + (a/(6.0*b))*(1-nu));
        const double k4 = C *((b/(3.0*a))*(1+nu) - (a/(3.0*b))*(1-nu));
        const double k5 = C *((b/(4.0*a))*(1-nu) + (a/(4.0*b))*(1-nu));
        const double k6 = C *((b/(4.0*a))*(1-nu) - (a/(4.0*b))*(1-nu));

        const double u1 = r_geom[0].GetSolutionStepValue(DISPLACEMENT_X);
        const double v1 = r_geom[0].GetSolutionStepValue(DISPLACEMENT_Y);
        const double u2 = r_geom[1].GetSolutionStepValue(DISPLACEMENT_X);
        const double v2 = r_geom[1].GetSolutionStepValue(DISPLACEMENT_Y);
        const double u3 = r_geom[2].GetSolutionStepValue(DISPLACEMENT_X);
        const double v3 = r_geom[2].GetSolutionStepValue(DISPLACEMENT_Y);
        const double u4 = r_geom[3].GetSolutionStepValue(DISPLACEMENT_X);
        const double v4 = r_geom[3].GetSolutionStepValue(DISPLACEMENT_Y);

        rRightHandSideVector(0) = aux * f_x - (k1*u1 + k2*u2 - k3*u3 - k4*u4);
        rRightHandSideVector(1) = aux * f_y - (k5*v1 + k6*v2 - k5*v3 - k6*v4);
        rRightHandSideVector(2) = aux * f_x - (k2*u1 + k1*u2 - k4*u3 - k3*u4);
        rRightHandSideVector(3) = aux * f_y - (k6*v1 + k5*v2 - k6*v3 - k5*v4);
        rRightHandSideVector(4) = aux * f_x - (-k3*u1 - k4*u2 + k1*u3 + k2*u4);
        rRightHandSideVector(5) = aux * f_y - (-k5*v1 - k6*v2 + k5*v3 + k6*v4);
        rRightHandSideVector(6) = aux * f_x - (-k4*u1 - k3*u2 + k2*u3 + k1*u4);
        rRightHandSideVector(7) = aux * f_y - (-k6*v1 - k5*v2 + k6*v3 + k5*v4);
    }

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const std::size_t n_nodes = r_geom.PointsNumber();
        const std::size_t local_size = 2 * n_nodes;
        if (rResult.size() != local_size) {
            rResult.resize(local_size);
        }
        for (unsigned int i = 0; i < n_nodes; ++i) {
            rResult[2*i] = r_geom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[2*i + 1] = r_geom[i].GetDof(DISPLACEMENT_Y).EquationId();
        }
    }

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& r_geom = GetGeometry();
        const SizeType n_nodes = r_geom.PointsNumber();
        const std::size_t local_size = 2 * n_nodes;
        if (rElementalDofList.size() != local_size) {
            rElementalDofList.resize(local_size);
        }
        for (unsigned int i = 0; i < n_nodes; i++) {
            rElementalDofList[2*i] = r_geom[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[2*i + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        }
    }

private:

    static constexpr double E = 1e6;

    static constexpr double nu = 0.3;

};

}

class SolvingStrategiesTestUtilities
{
public:

    static void SetUpTestModelPart1D(
        const std::size_t NumElems,
        const double ElemSize,
        ModelPart& rModelPart)
    {
        const int domain_size = 2;
        const std::size_t buffer_size = 3;
        rModelPart.SetBufferSize(buffer_size);
        rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

        rModelPart.AddNodalSolutionStepVariable(DISTANCE);

        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        for (std::size_t i_elem = 0; i_elem < NumElems; ++i_elem) {
            auto p_node_a = rModelPart.pGetNode(i_elem + 1);
            auto p_node_b = rModelPart.CreateNewNode(i_elem + 2, (i_elem + 1)*ElemSize, 0.0, 0.0);
            auto p_geom = Kratos::make_shared<Line2D2<Node>>(p_node_a, p_node_b);
            rModelPart.AddElement(Kratos::make_intrusive<TestLaplacianElement1D>(i_elem + 1, p_geom));
        }

        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISTANCE);
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISTANCE);
        }
    }

    static void SetUpTestModelPart2D(
        const std::size_t NumElemsX,
        const std::size_t NumElemsY,
        const double ElemSizeX,
        const double ElemSizeY,
        ModelPart& rModelPart)
    {
        const int domain_size = 2;
        const std::size_t buffer_size = 3;
        rModelPart.SetBufferSize(buffer_size);
        rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, domain_size);

        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        for (std::size_t i = 0; i < NumElemsY; ++i) {
            for (std::size_t j = 0; j < NumElemsX; ++j) {
                auto p_node_a = rModelPart.pGetNode(i * (NumElemsX + 1) + j + 1);
                auto p_node_b = rModelPart.CreateNewNode(i * (NumElemsX + 1) + (j + 1) + 1, p_node_a->X0() + ElemSizeX, p_node_a->Y0(), 0.0);
                auto p_node_c = rModelPart.CreateNewNode((i + 1) * (NumElemsX + 1) + (j + 1) + 1, p_node_a->X0() + ElemSizeX, p_node_a->Y0()+ ElemSizeY, 0.0);
                auto p_node_d = j == 0 ?
                    rModelPart.CreateNewNode((i + 1) * (NumElemsX + 1) + j + 1, p_node_a->X0(), p_node_a->Y0() + ElemSizeY, 0.0) :
                    rModelPart.pGetNode((i + 1) * (NumElemsX + 1) + j + 1);
                auto p_geom = Kratos::make_shared<Quadrilateral2D4<Node>>(p_node_a, p_node_b, p_node_c, p_node_d);
                rModelPart.AddElement(Kratos::make_intrusive<TestElasticityElement2D>(i * NumElemsX + j + 1, p_geom));
            }
        }

        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISPLACEMENT_X);
        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&DISPLACEMENT_Y);
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
        }
    }

};

}  // namespace Kratos::Testing.

