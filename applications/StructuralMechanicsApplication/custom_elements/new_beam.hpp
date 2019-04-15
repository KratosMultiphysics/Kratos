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

#if !defined(KRATOS_NEW_BEAM_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_NEW_BEAM_ELEMENT_3D2N_H_INCLUDED

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
 * @class NewBeamElement3D2N
 *
 * @brief This is a 3D-2node beam element with 3 translational dofs and 3 rotational dof per node
 *
 * @author Klaus B Sautter
 */

class NewBeamElement3D2N : public Element
{
protected:
    //const values
    static constexpr int msNumberOfNodes = 2;
    static constexpr int msDimension = 3;
    static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
    static constexpr unsigned int msElementSize = msLocalSize * 2;

public:
    KRATOS_CLASS_POINTER_DEFINITION(NewBeamElement3D2N);


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

    NewBeamElement3D2N() {};
    NewBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
    NewBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);


    ~NewBeamElement3D2N() override;

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

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo) override;



    void GetValuesVector(
        Vector& rValues,
        int Step = 0) override;

    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) override;

    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) override;


    int Check(const ProcessInfo& rCurrentProcessInfo) override;


    double CalculateShearModulus() const;
    double CalculatePsi(const double I, const double A_eff) const;

    BoundedMatrix<double, msLocalSize,msLocalSize>
        CalculateDeformationStiffness() const;

    BoundedMatrix<double, msElementSize,msElementSize>
        CreateElementStiffnessMatrix_Material() const;

    BoundedMatrix<double, msElementSize,msElementSize>
        CalculateInitialLocalCS() const;

    void AssembleSmallInBigMatrix(
        Matrix SmallMatrix,
        BoundedMatrix<double, msElementSize,
        msElementSize>& BigMatrix) const;

    BoundedMatrix<double, msElementSize,msLocalSize>
        CalculateTransformationS() const;

    Vector UpdateIncrementDeformation();

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    BoundedVector<double,msLocalSize> GetCurrentNodalPosition() const;

    BoundedMatrix<double, NewBeamElement3D2N::msDimension, msDimension>
        UpdateRotationMatrixLocal(Vector& Bisectrix,
            Vector& VectorDifference,const bool& rFinalize);

    Vector CalculateAntiSymmetricDeformationMode();
    Vector CalculateSymmetricDeformationMode();
    BoundedVector<double, msLocalSize> CalculateElementForces();
    Vector CalculateLocalNodalForces();
    Vector CalculateGlobalNodalForces();

private:

    Vector mTotalNodalDeformation = ZeroVector(msElementSize);
    Vector mQuaternionVEC_A = ZeroVector(msDimension);
    Vector mQuaternionVEC_B = ZeroVector(msDimension);
    double mQuaternionSCA_A = 1.00;
    double mQuaternionSCA_B = 1.00;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}

#endif
