//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_FEMDEM3D_HEXAHEDRON_ELEMENT_H_INCLUDED)
#define KRATOS_FEMDEM3D_HEXAHEDRON_ELEMENT_H_INCLUDED

#include "custom_elements/femdem3d_element.hpp"

namespace Kratos
{
class FemDem3DHexahedronElement : public FemDem3DElement
{

  public:
    /// Default constructors
    FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry);

    FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    FemDem3DHexahedronElement(FemDem3DHexahedronElement const &rOther);

    /// Destructor.
    virtual ~FemDem3DHexahedronElement();

    /// Assignment operator.
    FemDem3DHexahedronElement &operator=(FemDem3DHexahedronElement const &rOther);

    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const;

    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const;

    FemDem3DHexahedronElement()
    {
    }

    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override;

    // void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rValues,
        const ProcessInfo &rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(
        const Variable<Vector> &rVariable,
        std::vector<Vector> &rValues,
        const ProcessInfo &rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(
        const Variable<Matrix> &rVariable,
        std::vector<Matrix> &rValues,
        const ProcessInfo &rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(
        const Variable<double> &rVariable,
        std::vector<double> &rOutput,
        const ProcessInfo &rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(
        const Variable<Vector> &rVariable,
        std::vector<Vector> &rOutput,
        const ProcessInfo &rCurrentProcessInfo);

    void CalculateOnIntegrationPoints(
        const Variable<Matrix> &rVariable,
        std::vector<Matrix> &rOutput,
        const ProcessInfo &rCurrentProcessInfo);

    Matrix GetEdgeNodeNumbering() 
    {
        Matrix M = ZeroMatrix(mNumberOfEdges, 2);
        M(0, 0)  = 0; M(0, 1)  = 1;
        M(1, 0)  = 1; M(1, 1)  = 2;
        M(2, 0)  = 2; M(2, 1)  = 3;
        M(3, 0)  = 3; M(3, 1)  = 0;
        M(4, 0)  = 4; M(4, 1)  = 5;
        M(5, 0)  = 5; M(5, 1)  = 6;
        M(6, 0)  = 6; M(6, 1)  = 7;
        M(7, 0)  = 7; M(7, 1)  = 4;
        M(8, 0)  = 0; M(8, 1)  = 4;
        M(9, 0)  = 5; M(9, 1)  = 1;
        M(10, 0) = 6; M(10, 1) = 2;
        M(11, 0) = 7; M(11, 1) = 3;
        return M;
    }

    void CalculateAverageStressOnEdge(
        Vector& rEdgeStressVector,
        const unsigned int edge);

    void CalculateAverageStrainOnEdge(
        Vector& rEdgeStrainVector,
        const unsigned int edge);

    void CalculateCharacteristicLength(
        double& rcharacteristic_length, 
        const int Edge);

    double CalculateElementalDamage(const Vector &EdgeDamages);

    
  private:

    
}; // Class FemDem3DHexahedronElement

} // Namespace Kratos
#endif // KRATOS_FEMDEM3D_HEXAHEDRON_ELEMENT_H_INCLUDED  defined