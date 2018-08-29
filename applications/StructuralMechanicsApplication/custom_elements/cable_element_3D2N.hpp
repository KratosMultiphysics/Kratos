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

#if !defined(KRATOS_CABLE_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_CABLE_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
    /**
     * @class CableElement3D2N
     *
     * @brief This is a 3D-2node cable element with 3 translational dofs per node inheriting from the TrussElement3D2N
     *
     * @author Klaus B Sautter
     */
    class CableElement3D2N : public TrussElement3D2N
    {

    public:
        KRATOS_CLASS_POINTER_DEFINITION(CableElement3D2N);



        CableElement3D2N() {};
        CableElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry);
        CableElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);


        ~CableElement3D2N() override;

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

        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override;

        BoundedMatrix<double,msLocalSize,msLocalSize>
         CreateElementStiffnessMatrix(ProcessInfo& rCurrentProcessInfo) override;

        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function updates the internal normal force w.r.t. the current deformations
         * @param rinternalForces The current updated internal forces
         */
        void UpdateInternalForces(BoundedVector<double,msLocalSize>& rinternalForces) override;

    private:
        // boolean for the cable --> does not resist to compression
        bool mIsCompressed;

        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };


}


#endif
