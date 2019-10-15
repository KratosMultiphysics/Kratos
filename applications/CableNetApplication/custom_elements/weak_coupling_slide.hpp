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

#if !defined(KRATOS_WEAK_SLIDING_ELEMENT_H_INCLUDED )
#define  KRATOS_WEAK_SLIDING_ELEMENT_H_INCLUDED

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
     * @class WeakSlidingElement3D3N
     *
     * @brief This is a penalty element to realize sliding nodes element with 3 translational dofs per node
     *
     * @author Klaus B Sautter
     */

    class WeakSlidingElement3D3N : public Element
    {
    protected:
        //const values
        static constexpr int msNumberOfNodes = 3;
        static constexpr int msDimension = 3;
        static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;

    public:
        KRATOS_CLASS_POINTER_DEFINITION(WeakSlidingElement3D3N);


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


        WeakSlidingElement3D3N() {};
        WeakSlidingElement3D3N(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        WeakSlidingElement3D3N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~WeakSlidingElement3D3N() override;

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
