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

#if !defined(KRATOS_RING_ELEMENT_ROCCO_3D_H_INCLUDED )
#define  KRATOS_RING_ELEMENT_ROCCO_AV_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "cable_net_application_variables.h"
#define PI 3.1415926535898

namespace Kratos
{
    /**
     * @class RingElementROCCO3D
     *
     * @brief This is a ring elemen with 3 translational dofs per node. Volkwein 2004
     *
     * @author Klaus B Sautter
     */

    class KRATOS_API(CABLE_NET_APPLICATION) RingElementROCCO3D : public Element
    {
    protected:

    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RingElementROCCO3D);


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


        RingElementROCCO3D() {};
        RingElementROCCO3D(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        RingElementROCCO3D(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~RingElementROCCO3D() override;

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


        void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

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

        void CalculateLumpedMassVector(
            VectorType &rLumpedMassVector,
            const ProcessInfo& rCurrentProcessInfo) const override;

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


        void InternalForcesCircumference(VectorType &rRightHandSideVector);

        void InternalForcesDiagonal(VectorType &rRightHandSideVector);

        Vector GetCurrentLengthCircumferenceArray() const;
        Vector GetRefLengthCircumferenceArray() const;
        Vector GetDiagonalLengthArray(const int step = 0) const;

        Vector GetDirectionVectorCircumference() const;
        Vector GetDirectionVectorDiagonal() const;

        Vector GetDeltaPositions(const int& rDirection) const;

        Vector DistanceVectorNodes(const int node_a, const int node_b) const;

        int Check(const ProcessInfo& rCurrentProcessInfo) const override;

        /**
         * @brief This function checks if self weight is present
         */
        bool HasSelfWeight() const;

        /**
         * @brief This function calculates self-weight forces
         */
        Vector CalculateBodyForces();

        void CalculateOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rOutput,
            const ProcessInfo& rCurrentProcessInfo);

    private:

        double mDeformedLength = 0.0;
        double mMaxDeformedLength = 0.0;
        Vector mDeformedDiagonal = ZeroVector(2);
        Vector mMaxDeformedDiagonal = ZeroVector(2);
        double mPerimeterForce = 0.0;
        Vector mDiagonalForces = ZeroVector(2);

        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };


}


#endif
