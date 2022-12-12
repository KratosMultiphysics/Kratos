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

#if !defined(KRATOS_BRAKING_ELEMENT_H_INCLUDED )
#define  KRATOS_BRAKING_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"

namespace Kratos
{
    /**
     * @class BrakingElement3D1N
     *
     * @author Klaus B Sautter
     */

    class KRATOS_API(CABLE_NET_APPLICATION) BrakingElement3D1N : public Element
    {
    protected:
        //const values
        static constexpr int msNumberOfNodes = 1;
        static constexpr int msDimension = 3;
        static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;

    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BrakingElement3D1N);


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


        BrakingElement3D1N() {};
        BrakingElement3D1N(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        BrakingElement3D1N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~BrakingElement3D1N() override;

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

        void AddExplicitContribution(
            const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,
            const Variable<array_1d<double, 3>>& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo) override;

        void AddExplicitContribution(
            const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            const Variable<double >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateLumpedMassVector(
            VectorType &rLumpedMassVector,
            const ProcessInfo &rCurrentProcessInfo) const override;

        void CalculateDampingMatrix(MatrixType& rDampingMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        void GetValuesVector(
            Vector& rValues,
            int Step = 0) const override;

        void GetSecondDerivativesVector(
            Vector& rValues,
            int Step = 0) const override;

        void GetFirstDerivativesVector(
            Vector& rValues,
            int Step = 0) const override;

        void CalculateInternalForces(
            Vector& rInternalForces
        );

        int  Check(
            const ProcessInfo& rCurrentProcessInfo) const override;

private:

    double mAccumulatedPlasticDisplacement = 0.0;
    double mAccumulatedPlasticAlpha = 0.0;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
    };


}


#endif
