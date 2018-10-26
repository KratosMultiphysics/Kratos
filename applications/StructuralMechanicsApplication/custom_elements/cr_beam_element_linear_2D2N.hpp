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

#if !defined(KRATOS_CR_BEAM_ELEMENT_LINEAR_2D2N_H_INCLUDED )
#define  KRATOS_CR_BEAM_ELEMENT_LINEAR_2D2N_H_INCLUDED

// System includes

// External includes


// Project includes
#include "custom_elements/cr_beam_element_2D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{
    /**
     * @class CrBeamElementLinear2D2N
     *
     * @brief This is a linear 2D-2node beam element with 2 translational dofs and 1 rotational dof per node inheriting from CrBeamElement2D2N
     *
     * @author Klaus B Sautter
     */

    class CrBeamElementLinear2D2N : public CrBeamElement2D2N
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(CrBeamElementLinear2D2N);


        CrBeamElementLinear2D2N(IndexType NewId, GeometryType::Pointer pGeometry);
        CrBeamElementLinear2D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);

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


        ~CrBeamElementLinear2D2N() override;

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


    /////////////////////////////////////////////////
    ///////////// CUSTOM FUNCTIONS --->>
    /////////////////////////////////////////////////

        /**
         * @brief This function calculates the transformation matrix to globalize/localize vectors and/or matrices
         */
        BoundedMatrix<double,msElementSize,msElementSize> CreateRotationMatrix() override;

        /**
         * @brief This function calculates the reference length
         */
        double CalculateLength() override;


        void GetValueOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(
            const Variable<array_1d<double, 3 > >& rVariable,
            std::vector< array_1d<double, 3 > >& rOutput,
            const ProcessInfo& rCurrentProcessInfo) override;

    private:
        CrBeamElementLinear2D2N() {};

        // stores the globalized master stiffness matrix for reaction forces
        Matrix mK_Master = ZeroMatrix(msElementSize,msElementSize);


        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };

}

#endif
