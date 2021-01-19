// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

#if !defined(KRATOS_BEAM_ELEMENT_2D2N_H_INCLUDED )
#define  KRATOS_BEAM_ELEMENT_2D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{
/**
 * @class BeamElement2D2N
 *
 * @brief This is a 2D-2node beam element with 3 translational dofs and 3 rotational dof per node
 *       from "Co-rotational beam elements in instability problems - Jean-Marc Battini"
 *       and "2D Corotational Beam Formulation - Louie L. Yaw"
 *
 * @author Klaus B Sautter
 */

class BeamElement2D2N : public Element
{
protected:
    //const values
    static constexpr int msNumberOfNodes = 2;
    static constexpr int msDimension = 2;
    static constexpr unsigned int msLocalSize = 3;
    static constexpr unsigned int msElementSize = msLocalSize * msNumberOfNodes;

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BeamElement2D2N);


    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::DofsVectorType DofsVectorType;

    enum class ConfigurationType {
      Current,
      Reference
    };

    BeamElement2D2N() {};
    BeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry);
    BeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);


    ~BeamElement2D2N() override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param pGeom The pointer to the geometry of the element
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param ThisNodes The array containing nodes
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    Matrix CalculateInitialLocalCS() const;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void ConstCalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) const;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    void ConstCalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) const;

    void ConstCalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) const;


    void GetValuesVector(
        Vector& rValues,
        int Step = 0) const override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;


    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    double RigidBodyOrientation(const ConfigurationType& rConfiguration) const;

    void CalculateBMatrix(Matrix& rBMatrix) const;

    void LocalMaterialMatrix(Matrix& rMaterialMatrix) const;

    void MaterialTangentMatrix(Matrix& rMaterialMatrix) const;

    void GeometricTangentMatrix(Matrix& rGeometricMatrix) const;

    void GetLocalDeformations(Vector& rLocalDeformation) const;

    void GetLocalInternalForces(Vector& rLocalForces) const;

    void GetGlobalInternalForces(Vector& rGlobalForces) const;

    void CalculateMassMatrix(MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void AddExplicitContribution(
        const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo) override;

    BoundedVector<double, msElementSize> CalculateBodyForces() const;

    void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
        const BoundedVector<double, 3> ForceInput,
        BoundedVector<double, msElementSize>
        & rRightHandSideVector,
        const double GeometryLength) const;

private:


    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

};


}

#endif