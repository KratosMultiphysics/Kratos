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

#if !defined(KRATOS_CR_BEAM_ELEMENT_LINEAR_3D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_LINEAR_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/cr_beam_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{
    /**
     * @class CrBeamElementLinear3D2N
     *
     * @brief This is a linear 3D-2node beam element with 3 translational dofs and 3 rotational dof per node inheriting from CrBeamElement3D2N
     *
     * @author Klaus B Sautter
     */

    class CrBeamElementLinear3D2N : public CrBeamElement3D2N
    {

        public:
            KRATOS_CLASS_POINTER_DEFINITION(CrBeamElementLinear3D2N);

            CrBeamElementLinear3D2N() {};
            CrBeamElementLinear3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
            CrBeamElementLinear3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties);


            ~CrBeamElementLinear3D2N() override;

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
             * @brief This function calculates the element stiffness w.r.t. deformation modes
             */
            BoundedMatrix<double,msLocalSize,msLocalSize> CalculateDeformationStiffness() override;

            void CalculateOnIntegrationPoints(
                const Variable<array_1d<double, 3 > >& rVariable,
                std::vector< array_1d<double, 3 > >& rOutput,
                const ProcessInfo& rCurrentProcessInfo) override;

            void CalculateOnIntegrationPoints(
                const Variable<Vector >& rVariable,
                std::vector< Vector >& rOutput,
                const ProcessInfo& rCurrentProcessInfo) override;

        private:

            friend class Serializer;
            void save(Serializer& rSerializer) const override;
            void load(Serializer& rSerializer) override;
    };

}

#endif
