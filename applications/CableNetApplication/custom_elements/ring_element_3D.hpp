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

#if !defined(KRATOS_RING_ELEMENT_3D_H_INCLUDED )
#define  KRATOS_RING_ELEMENT_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"

namespace Kratos
{
    /**
     * @class RingElement3D
     *
     * @brief This is a ring elemen with 3 translational dofs per node
     *
     * @author Klaus B Sautter
     */

    class KRATOS_API(CABLE_NET_APPLICATION) RingElement3D : public Element
    {
    protected:

    public:
        KRATOS_CLASS_POINTER_DEFINITION(RingElement3D);


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


        RingElement3D() {};
        RingElement3D(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        RingElement3D(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~RingElement3D() override;

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

        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
            VectorType &rRightHandSideVector,
            const ProcessInfo &rCurrentProcessInfo) override;

        void CalculateRightHandSide(VectorType &rRightHandSideVector,
            const ProcessInfo &rCurrentProcessInfo) override;

        void GetValuesVector(Vector& rValues,int Step = 0) const override;
        void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;
        void GetFirstDerivativesVector(Vector& rValues,int Step = 0) const override;

        Matrix ElasticStiffnessMatrix() const;
        Matrix GeometricStiffnessMatrix() const;
        inline Matrix TotalStiffnessMatrix() const;

        double GetCurrentLength() const;
        double GetRefLength() const;
        double CalculateGreenLagrangeStrain() const;
        double LinearStiffness() const;

        void CalculateLumpedMassVector(VectorType &rMassVector);

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;


        void AddExplicitContribution(
            const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            const Variable<double >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        void AddExplicitContribution(const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            const Variable<array_1d<double, 3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        int Check(const ProcessInfo& rCurrentProcessInfo) const override;

        /**
         * @brief This function checks if self weight is present
         */
        bool HasSelfWeight() const;

        /**
         * @brief This function calculates self-weight forces
         */
        Vector CalculateBodyForces();

    private:

        Vector GetCurrentLengthArray() const;
        Vector GetRefLengthArray() const;
        Vector GetDeltaPositions(const int& rDirection) const;
        Vector GetDirectionVectorNt() const;
        Vector GetInternalForces() const;

        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };


}


#endif
