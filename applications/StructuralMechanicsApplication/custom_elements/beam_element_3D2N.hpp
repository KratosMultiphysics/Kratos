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

#if !defined(KRATOS_BEAM_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_BEAM_ELEMENT_3D2N_H_INCLUDED

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
 * @class BeamElement3D2N
 *
 * @brief This is a 3D-2node beam element with 3 translational dofs and 3 rotational dof per node
 *       from "Co-rotational beam elements in instability problems - Jean-Marc Battini"
 *
 * @author Klaus B Sautter
 */

class BeamElement3D2N : public Element
{
protected:
    //const values
    static constexpr int msNumberOfNodes = 2;
    static constexpr int msDimension = 3;
    static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
    static constexpr unsigned int msElementSize = msLocalSize * 2;

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BeamElement3D2N);


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

    BeamElement3D2N() {};
    BeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
    BeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);


    ~BeamElement3D2N() override;

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

    /**
     * @brief This function calculates the elastic part of the total stiffness matrix
     */
    Matrix CreateElementStiffnessMatrix_Material() const;
    virtual Matrix CreateElementStiffnessMatrixIntermediate() const;
    virtual Matrix GlobalTangentStiffnessMatrix() const;

    /**
     * @brief This function calculates the current nodal position
     */
    virtual BoundedVector<double,msLocalSize> GetCurrentNodalPosition() const;


    /**
     * @brief This function calculates the initial transformation matrix to globalize/localize vectors and/or matrices
     */
    Matrix CalculateInitialLocalCS() const;


    Vector CurrentLocalAxis1() const;
    virtual Matrix CoRotatingCS() const;
    Vector LocalDeformations() const;
    Matrix LogRotationMatrix(const Matrix& rRotationMatrix) const;
    Matrix SkewSymmetricMatrix(const Vector& rinput_vec) const;
    void UpdateGlobalNodalRotations();
    Vector LocalInternalForces() const;
    Matrix DMatrix() const;
    Matrix GMatrix() const;
    Vector RVector() const;
    Vector AVector() const;
    Matrix PMatrix() const;
    Matrix EMatrix() const;
    Matrix QMatrix() const;
    Matrix BMatrixGlobal() const;
    Matrix BMatrixIntermediateA() const;
    Matrix InverseLocalRotation(const int node_nr) const;
    Vector LocalInternalIntermediateForces() const;


    void CalculateConsistentMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const;
    void CalculateLumpedMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) const;
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;
    void AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,const ProcessInfo& rCurrentProcessInfo) override;


    void GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo) override;

    IntegrationMethod GetIntegrationMethod() const override;

    void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
    const BoundedVector<double, msDimension> ForceInput,
    BoundedVector<double, msElementSize>
    & rRightHandSideVector,const double GeometryLength) const;

    BoundedVector<double, msElementSize> CalculateBodyForces() const;

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


    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This function calculates shear modulus from user input values
     */
    double CalculateShearModulus() const;

    virtual Vector CalculateGlobalNodalForces() const;

    Vector GetIncrementDeformation() const;

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


private:


    Vector mDeformationCurrentIteration = ZeroVector(msElementSize);
    Vector mDeformationPreviousIteration = ZeroVector(msElementSize);
    Matrix mGlobalRotationNode1 = IdentityMatrix(msDimension);
    Matrix mGlobalRotationNode2 = IdentityMatrix(msDimension);


    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

};


}

#endif