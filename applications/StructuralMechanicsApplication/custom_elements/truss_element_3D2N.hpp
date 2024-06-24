// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
    /**
     * @class TrussElement3D2N
     *
     * @brief This is a 3D-2node truss element with 3 translational dofs per node
     *
     * @author Klaus B Sautter
     */

    class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TrussElement3D2N : public Element
    {
    protected:
        //const values
        static constexpr int msNumberOfNodes = 2;
        static constexpr int msDimension = 3;
        static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
        ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TrussElement3D2N);


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


        TrussElement3D2N() {};
        TrussElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        TrussElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~TrussElement3D2N() override;

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
            const ProcessInfo& rCurrentProcessInfo) const override;

        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo) const override;

        void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function calculates the total stiffness matrix for the element
         */
        virtual BoundedMatrix<double,msLocalSize,msLocalSize>
         CreateElementStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo);

        void Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

        void Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function updates the internal normal force w.r.t. the current deformations
         * @param rinternalForces The current updated internal forces
         */
        virtual void UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces, const ProcessInfo& rCurrentProcessInfo);

        /**
         * @brief This function calculates the transformation matrix to globalize vectors and/or matrices
         * @param rRotationMatrix The transformation matrix
         */
        void CreateTransformationMatrix(BoundedMatrix<double,msLocalSize,msLocalSize>& rRotationMatrix);


        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;


        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateConsistentMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo) const;

        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable rDestinationVariable (double version)
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(
        const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<double >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function is designed to make the element to assemble an rRHS vector identified by a variable rRHSVariable by assembling it to the nodes on the variable (array_1d<double, 3>) version rDestinationVariable.
     * @details The "AddEXplicit" FUNCTIONS THE ONLY FUNCTIONS IN WHICH AN ELEMENT IS ALLOWED TO WRITE ON ITS NODES.
     * The caller is expected to ensure thread safety hence SET/UNSETLOCK MUST BE PERFORMED IN THE STRATEGY BEFORE CALLING THIS FUNCTION
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSVector will be assembled
     * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHSVector,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double, 3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


        void GetValuesVector(
            Vector& rValues,
            int Step = 0) const override;

        void GetSecondDerivativesVector(
            Vector& rValues,
            int Step = 0) const override;

        void GetFirstDerivativesVector(
            Vector& rValues,
            int Step = 0) const override;

        int  Check(
            const ProcessInfo& rCurrentProcessInfo) const override;

        /**
         * @brief This function calculates the current Green-Lagrange strain
         */
        double CalculateGreenLagrangeStrain()const;

        /**
         * @brief This function calculates self-weight forces
         */
        BoundedVector<double,msLocalSize> CalculateBodyForces();

        /**
         * @brief This function assembles the geometric stiffness part of the total stiffness matrix
         * @param rGeometricStiffnessMatrix The geometric stiffness matrix
         * @param rCurrentProcessInfo The current process information
         */
        void CalculateGeometricStiffnessMatrix(BoundedMatrix<double,msLocalSize,msLocalSize>& rGeometricStiffnessMatrix,
            const ProcessInfo& rCurrentProcessInfo);

        /**
         * @brief This function assembles the elastic stiffness part of the total stiffness matrix
         * @param rElasticStiffnessMatrix The elastic stiffness matrix
         * @param rCurrentProcessInfo The current process information
         */
        void CalculateElasticStiffnessMatrix(BoundedMatrix<double,msLocalSize,msLocalSize>& rElasticStiffnessMatrix,
            const ProcessInfo& rCurrentProcessInfo);

        /**
         * @brief This function calculates the current nodal postion for the transformation matrix
         * @param rReferenceCoordinates The current coordinates
         */
        virtual void WriteTransformationCoordinates(
            BoundedVector<double,msLocalSize>& rReferenceCoordinates);

        virtual double ReturnTangentModulus1D(const ProcessInfo& rCurrentProcessInfo);

        void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function checks if self weight is present
         */
        bool HasSelfWeight() const;

        const Parameters GetSpecifications() const override;

protected:
    double ReturnTangentModulus1D(double Strain, const ProcessInfo& rCurrentProcessInfo) const;

private:
    /**
     * @brief This method computes directly the lumped mass vector
     * @param rMassVector The lumped mass vector
     */
    void CalculateLumpedMassVector(
        VectorType& rMassVector,
        const ProcessInfo& rCurrentProcessInfo) const override;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
    };


}
