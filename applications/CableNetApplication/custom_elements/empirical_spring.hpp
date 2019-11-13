//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

#if !defined(KRATOS_EMPIRICAL_SPRING_H_INCLUDED )
#define  KRATOS_EMPIRICAL_SPRING_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"

namespace Kratos
{
    /**
     * @class EmpiricalSpringElement3D2N
     *
     * @brief This spring reads a fitted polynomial as displacement-load curve
     *          takes u and returns f
     *
     * @author Klaus B Sautter
     */

    class EmpiricalSpringElement3D2N : public Element
    {
    protected:
        //const values
        static constexpr int msNumberOfNodes = 2;
        static constexpr int msDimension = 3;
        static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;

    public:
        KRATOS_CLASS_POINTER_DEFINITION(EmpiricalSpringElement3D2N);


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


        EmpiricalSpringElement3D2N() {};
        EmpiricalSpringElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        EmpiricalSpringElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~EmpiricalSpringElement3D2N() override;

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

        /**
         * @brief This function calculates the total stiffness matrix for the element
         */
        virtual BoundedMatrix<double,msLocalSize,msLocalSize>
         CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo);


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

        void AddExplicitContribution(
            const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
            Variable<array_1d<double, 3>>& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo) override;

        void AddExplicitContribution(
            const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            Variable<double >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

        void CalculateLumpedMassVector(VectorType& rMassVector);

        void CalculateDampingMatrix(MatrixType& rDampingMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function evaluates a given polynomial
         * p(x) = p[0] * x**deg + ... + p[deg]
         */
        double EvaluatePolynomial(const Vector& rPolynomial) const;
        double EvaluatePolynomialFirstDerivative(const Vector& rPolynomial) const;

        void LocalizeVector(BoundedVector<double,msLocalSize>& rInputVector);
        void GlobalizeVector(BoundedVector<double,msLocalSize>& rInputVector);
        void GlobalizeMatrix(BoundedMatrix<double,msLocalSize,msLocalSize>& rInputMatrix);

        void WriteTransformationCoordinates(BoundedVector<double,msLocalSize>
            & rReferenceCoordinates);

        void CreateTransformationMatrix(BoundedMatrix<double,msLocalSize,
            msLocalSize>& rRotationMatrix);

        double GetElementElongation() const;


        void GetValuesVector(
            Vector& rValues,
            int Step = 0) override;

        void GetSecondDerivativesVector(
            Vector& rValues,
            int Step = 0) override;

        void GetFirstDerivativesVector(
            Vector& rValues,
            int Step = 0) override;

        int  Check(
            const ProcessInfo& rCurrentProcessInfo) override;

private:

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
    };


}


#endif
