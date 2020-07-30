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

#if !defined(KRATOS_SLIDING_CABLE_ELEMENT_H_INCLUDED )
#define  KRATOS_SLIDING_CABLE_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
    /**
     * @class SlidingCableElement3D
     *
     * @brief This is a sliding node element with 3 translational dofs per node
     *
     * @author Klaus B Sautter
     */

    class SlidingCableElement3D : public Element
    {
    protected:

        ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    public:
        KRATOS_CLASS_POINTER_DEFINITION(SlidingCableElement3D);


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


        SlidingCableElement3D() {};
        SlidingCableElement3D(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        SlidingCableElement3D(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~SlidingCableElement3D() override;

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

        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

        void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
            VectorType &rRightHandSideVector,
            ProcessInfo &rCurrentProcessInfo) override;

        void CalculateRightHandSide(VectorType &rRightHandSideVector,
            ProcessInfo &rCurrentProcessInfo) override;

        void GetValuesVector(Vector& rValues,int Step = 0) override;
        void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;
        void GetFirstDerivativesVector(Vector& rValues,int Step = 0) override;

        Vector GetCurrentLengthArray(int Step = 0) const;
        Vector GetRefLengthArray() const;
        Vector GetDeltaPositions(const int& rDirection) const;
        Vector GetDirectionVectorNt() const;
        Vector GetInternalForces();

        Matrix ElasticStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo) const;
        Matrix GeometricStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo) const;
        inline Matrix TotalStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo) const;

        double GetCurrentLength() const;
        double GetRefLength() const;
        double CalculateGreenLagrangeStrain() const;
        inline double LinearStiffness() const
        {
            return (this->GetProperties()[CROSS_AREA] * this->GetProperties()[YOUNG_MODULUS] / this->GetRefLength());
        };

        double GetPK2PrestressValue() const
        {
            double prestress = 0.00;
            if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
                prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
            }
            return prestress;
        };

        void CalculateLumpedMassVector(VectorType &rMassVector);

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) override;

        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            ProcessInfo& rCurrentProcessInfo) override;


        void AddExplicitContribution(
            const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            Variable<double >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        void AddExplicitContribution(const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            Variable<array_1d<double, 3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        int Check(const ProcessInfo& rCurrentProcessInfo) override;


        void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
        void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;
        void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


        void GetConstitutiveLawTrialResponse(
            const ProcessInfo& rCurrentProcessInfo);

        /**
         * @brief This function checks if self weight is present
         */
        bool HasSelfWeight() const;

        /**
         * @brief This function calculates self-weight forces
         */
        Vector CalculateBodyForces();

        Vector GetCustomInternalForceWithFriction(const Vector& rNormalForces);

        Vector CalculateProjectionLengths();

        double ReturnTangentModulus1D(const ProcessInfo& rCurrentProcessInfo) const;

    private:

        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };


}


#endif
