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

#if !defined(KRATOS_CR_BEAM_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_3D2N_H_INCLUDED

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
     * @class CrBeamElement3D2N
     *
     * @brief This is a 3D-2node beam element with 3 translational dofs and 3 rotational dof per node
     *
     * @author Klaus B Sautter
     */

    class CrBeamElement3D2N : public Element
    {
    protected:
        //const values
        static constexpr int msNumberOfNodes = 2;
        static constexpr int msDimension = 3;
        static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
        static constexpr unsigned int msElementSize = msLocalSize * 2;

    public:
        KRATOS_CLASS_POINTER_DEFINITION(CrBeamElement3D2N);


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

        CrBeamElement3D2N() {};
        CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
        CrBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~CrBeamElement3D2N() override;

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
        BoundedMatrix<double,msElementSize,msElementSize> CreateElementStiffnessMatrix_Material();

        /**
         * @brief This function calculates the geometric part of the total stiffness matrix
         */
        BoundedMatrix<double,msElementSize,msElementSize>  CreateElementStiffnessMatrix_Geometry();

        /**
         * @brief This function calculates the element stiffness w.r.t. deformation modes
         */
        virtual BoundedMatrix<double,msLocalSize,msLocalSize> CalculateDeformationStiffness();

        /**
         * @brief This function calculates a transformation matrix from deformation modes to real deformations
         */
        BoundedMatrix<double,msElementSize,msLocalSize> CalculateTransformationS();

        /**
         * @brief This function calculates the current nodal position
         */
        BoundedVector<double,msLocalSize> GetCurrentNodalPosition();

        /**
         * @brief This function calculates the internal element forces
         */
        BoundedVector<double,msLocalSize> CalculateElementForces(const Vector& Bisectrix,const Vector& VectorDifference);


        /**
         * @brief This function calculates the transformation matrix to globalize/localize vectors and/or matrices
         * @param rRotationMatrix The current transformation matrix
         */
        void CalculateTransformationMatrix(
            BoundedMatrix<double,msElementSize,msElementSize>& rRotationMatrix,
            Vector& Bisectrix, Vector& VectorDifference);


        /**
         * @brief This function calculates the initial transformation matrix to globalize/localize vectors and/or matrices
         */
        BoundedMatrix<double,msElementSize,msElementSize> CalculateInitialLocalCS();


        /**
         * @brief This function updates constantly the transformation matrix
         */
        BoundedMatrix<double,msDimension,msDimension> UpdateRotationMatrixLocal(Vector& Bisectrix, Vector& VectorDifference);

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

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) override;


        /**
         * @brief This function calculates the lumped mass matrix
         * @param rMassMatrix The current mass matrix
         * @param rCurrentProcessInfo The current Process information
         */
        void CalculateLumpedMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo);


        /**
         * @brief This function calculates the consistent mass matrix
         * @param rMassMatrix The current mass matrix
         * @param rCurrentProcessInfo The current Process information
         */
        void CalculateConsistentMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo);


        /**
         * @brief This function calculates parts of the total consistent mass matrix to simplify the code
         * @param rMassMatrix The current mass matrix
         * @param Phi The reduction value in case of shear-deformable structures
         * @param CT A scaling factor
         * @param CR A scaling factor
         * @param L The element length
         * @param dir The direction of the current cs
         */
        void BuildSingleMassMatrix(
            MatrixType& rMassMatrix,
            const double Phi, const double CT, const double CR, const double L, const double dir);

        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

        void AddExplicitContribution(const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            Variable<array_1d<double, 3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo) override;

        void GetValuesVector(
            Vector& rValues,
            int Step = 0) override;

        void GetSecondDerivativesVector(
            Vector& rValues,
            int Step = 0) override;

        void GetFirstDerivativesVector(
            Vector& rValues,
            int Step = 0) override;

        /**
         * @brief This function is used to assemble single transformation matrix in the big global rotation matrix
         * @param SmallMatrix The local transformation matrix
         * @param BigMatrix The total global rotation matrix
         */
        void AssembleSmallInBigMatrix(Matrix SmallMatrix, BoundedMatrix<double,
            msElementSize,msElementSize>& BigMatrix);

        int Check(const ProcessInfo& rCurrentProcessInfo) override;


        /**
         * @brief This function calculates reduction values in case of shear-deformable structures
         * @param I The second moment of area
         * @param A_eff The shear-effective area
         */
        double CalculatePsi(const double I, const double A_eff);

        /**
         * @brief This function calculates shear modulus from user input values
         */
        double CalculateShearModulus();

        /**
         * @brief This function calculates the reference length
         */
        double CalculateReferenceLength();

        /**
         * @brief This function calculates the current length
         */
        double CalculateCurrentLength();

        /**
         * @brief This function updates incremental deformation w.r.t. to current and previous deformations
         */
        Vector UpdateIncrementDeformation();


        /**
         * @brief This function calculates self-weight forces
         */
        BoundedVector<double,msElementSize> CalculateBodyForces();

        void CalculateOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void GetValueOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<Vector >& rVariable,
            std::vector< Vector >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void GetValueOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;


        IntegrationMethod GetIntegrationMethod() const override;


        /**
         * @brief This function calculates nodal moments due to self-weight
         * @param ForceInput The self-weight line load vector
         * @param rRightHandSideVector The right hand side of the problem
         * @param GeometryLength The element length
         */
        void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
            const BoundedVector<double,msDimension> ForceInput,
            BoundedVector<double,msElementSize>& rRightHandSideVector,
            const double GeometryLength);


        /**
         * @brief This function calculates the symmetric deformation modes
         * @param VectorDifference The vector differences of the quaternions
         */
        Vector CalculateSymmetricDeformationMode(const Vector& VectorDifference);

        /**
         * @brief This function calculates the antisymmetric deformation modes
         * @param Bisectrix The bisectrix between the local axis1 from the last iter. step and the updated axis 1
         */
        Vector CalculateAntiSymmetricDeformationMode(const Vector& Bisectrix);

        /**
         * @brief This function calculates the local nodal forces
         * @param Bisectrix The bisectrix between the local axis1 from the last iter. step and the updated axis 1
         * @param VectorDifference The vector differences of the quaternions
         */
        void CalculateLocalNodalForces(const Vector& Bisectrix,const Vector& VectorDifference);

    private:

        int mIterationCount = 0;
        Vector mTotalNodalDeformation = ZeroVector(msElementSize); // save as the displacement from the last iteration step is needed
        Matrix mLocalRotationMatrix  = ZeroMatrix(msDimension); // save this as updating the matrix takes rather long
        Vector mQuaternionVEC_A = ZeroVector(msDimension);
        Vector mQuaternionVEC_B = ZeroVector(msDimension);
        double mQuaternionSCA_A = 1.00;
        double mQuaternionSCA_B = 1.00;
        Vector mNodalForces = ZeroVector(msElementSize);




        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;


    public:
        void IncrementIterationCounter() {this->mIterationCount += 1;};
    };


}

#endif
