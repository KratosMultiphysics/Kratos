// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//  license:      structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

#if !defined(KRATOS_CR_BEAM_ELEMENT_2D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_2D2N_H_INCLUDED

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
     * @class CrBeamElement2D2N
     *
     * @brief This is a 2D-2node beam element with 2 translational dofs and 1 rotational dof per node
     *
     * @author Klaus B Sautter
     */

    class CrBeamElement2D2N : public Element
    {
    protected:
        //const values
        static constexpr int msNumberOfNodes = 2;
        static constexpr int msDimension = 2;
        static constexpr unsigned int msLocalSize = 3;
        static constexpr unsigned int msElementSize = msLocalSize * 2;

    public:
        KRATOS_CLASS_POINTER_DEFINITION(CrBeamElement2D2N);


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

        CrBeamElement2D2N() {};
        CrBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry);
        CrBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~CrBeamElement2D2N() override;


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

        void GetValuesVector(
            Vector& rValues,
            int Step = 0) override;

        void GetSecondDerivativesVector(
            Vector& rValues,
            int Step = 0) override;

        void GetFirstDerivativesVector(
            Vector& rValues,
            int Step = 0) override;

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

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

        void AddExplicitContribution(const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            Variable<array_1d<double, 3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo) override;

        int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /////////////////////////////////////////////////
    ///////////// CUSTOM FUNCTIONS --->>
    /////////////////////////////////////////////////

        /**
         * @brief This function calculates shear modulus from user input values
         */
        double CalculateShearModulus();


        /**
         * @brief This function calculates reduction values in case of shear-deformable structures
         * @param I The second moment of area
         * @param A_eff The shear-effective area
         */
        double CalculatePsi(const double I, const double A_eff);

        /**
         * @brief This function calculates the initial angle
         */
        double CalculateInitialElementAngle();

        /**
         * @brief This function calculates the current angle
         */
        double CalculateDeformedElementAngle();

        /**
         * @brief This function calculates self-weight forces
         */
        BoundedVector<double,msElementSize> CalculateBodyForces();

        /**
         * @brief This function calculates nodal moments due to self-weight
         * @param ForceInput The self-weight line load vector
         * @param rRightHandSideVector The right hand side of the problem
         * @param GeometryLength The element length
         */
        void CalculateAndAddWorkEquivalentNodalForcesLineLoad(
            const BoundedVector<double,3> ForceInput,
            BoundedVector<double,msElementSize>& rRightHandSideVector,
            const double GeometryLength);

        IntegrationMethod GetIntegrationMethod() const override;

        /**
         * @brief This function calculates a transformation matrix from deformation modes to real deformations
         */
        BoundedMatrix<double,msElementSize,msLocalSize> CalculateTransformationS();

        /**
         * @brief This function calculates the current length
         */
        virtual double CalculateLength();

        /**
         * @brief This function calculates the reference length
         */
        double CalculateReferenceLength();


        /**
         * @brief This function calculates the elastic part of the total stiffness matrix
         */
        BoundedMatrix<double,msLocalSize,msLocalSize> CreateElementStiffnessMatrix_Kd_mat();

        /**
         * @brief This function calculates the geometric part of the total stiffness matrix
         */
        BoundedMatrix<double,msLocalSize,msLocalSize> CreateElementStiffnessMatrix_Kd_geo();

        /**
         * @brief This function calculates the co-rotating part of the total stiffness matrix
         */
        BoundedMatrix<double,msElementSize,msElementSize> CreateElementStiffnessMatrix_Kr();

        /**
         * @brief This function assembles the total stiffness matrix
         */
        BoundedMatrix<double,msElementSize,msElementSize> CreateElementStiffnessMatrix_Total();


        /**
         * @brief This function globalizes matrices
         * @param A The matrix to be globalized
         */
        void GlobalizeMatrix(Matrix &A);

        /**
         * @brief This function globalizes vectors
         * @param A The vector to be globalized
         */
        void GlobalizeVector(Vector &A);

        /**
         * @brief This function calculates the modulus to 2 PI to keep the rotation angle between 0 and 2PI
         * @param A The current angle
         */
        double Modulus2Pi(double A);

        /**
         * @brief This function calculates the transformation matrix to globalize/localize vectors and/or matrices
         */
        virtual BoundedMatrix<double,msElementSize,msElementSize> CreateRotationMatrix();


        /**
         * @brief This function calculates the deformation parameters
         */
        BoundedVector<double,msLocalSize> CalculateDeformationParameters();

        /**
         * @brief This function calculates the internal forces w.r.t. the deformation parameters
         */
        BoundedVector<double,msLocalSize> CalculateInternalStresses_DeformationModes();

        /**
         * @brief This function calculates the "real" internal forces in a local reference frame
         */
        BoundedVector<double,msElementSize> ReturnElementForces_Local();


        /**
         * @brief This function calculates the element contributions to an explicit time integration
         */

        void GetValueOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

    private:

        // stores the deformation modes
        BoundedVector<double,msLocalSize> mDeformationForces = ZeroVector(msLocalSize);

        // stores the globalized internal forces for calculation of the residual
        Vector mInternalGlobalForces = ZeroVector(msElementSize);


        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;

    };

}

#endif
