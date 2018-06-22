// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: L. Chen
//
//
//

#if !defined( KRATOS_FIBER_BEAM_COLUMN_ELEMENT_3D2N_INCLUDED )
#define KRATOS_FIBER_BEAM_COLUMN_ELEMENT_3D2N_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_element/fiber_beam_column_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variable.h"
#include "includes/serializer.h"

namespace KRATOS
{
    /**
     * @class FiberBeamColumnElement3D2N
     *
     * @brief This is a Fiber Beam-Column Element for reinforced concrete structures
     *
     * @author Long Chen
     */

    class FiberBeamColumnElement3D2N:public FiberBeamColumnElement3D2N
    {
        protected:
        //const values
        static constexpr int msNumberOfNodes = 2;
        static constexpr int msDimension = 3;
        static constexpr unsigned int msLocalSize = 5;
        static constexpr unsigned int msElementSize = 10;




        public:
        KRATOS_CLASS_POINTER_DEFINITION(FiberBeamColumnElement3D2N);

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


        FiberBeamColumnElement3D2N() {};
        FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
        FiberBeamColumnElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties);

        ~FiberBeamColumnElement3D2N() override;
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
            DofsVectorType& rElementalDofList,
            ProcessInfo& rCurrentProcessInfo) override;
        )

        void Initialize () override;


        /**
         * @brief This function calculates the transformation matrix to globalize/localize vectors and/or matrices
         * @param rRotationMatrix The current transformation matrix
         */
        void CalculateTransformationMatrix(
            BoundedMatrix<double, msElementSize, msElementSize>& rRotationMatrix,
            Vector& Bisectrix, Vector& VectorDifference);
        )
        /**
         * @brief This function calculates the initial transformation matrix to globalize/localize vectors and/or matrices
         */
        BoundedMatrix<double, msElementSize, msElementSize> CalculateInitialLocalCS();

        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override;
        )

        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) override;
        )

        void CalculateLeftHandSide(
            MatrixType& rLeftHandsideMatrix,
            ProcessInfo& rCurrentProcessInfo) override;
        )

        void CalclateMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) override;
        )

        void getValuesVector(
            Vector& rValues,
            int Step = 0 ) override;


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
        void CalculateLocalNodalForces(const Vector& Bisectrix, const Vector& VectorDifference);
            

    private:

    int mIterationCount = 0;        // iteration counter j used in element state determination
    Vector mTotalNodalDeformation = ZeroVector(msLocalSize); // save as the nodal displacement from the last iteration step j-1
    vector mNodalForces = ZeroVector(msLocalSize); // save as the nodal displacement from the last iteration step j-1

    Vector mChangeInElementDeformationIncr = ZeroVector(msLocalSize);  // Change in the element deformation increments ddq_i
    Vector mChangeInElementForceIncr = ZeroVector(msLocalSize); // Change in the element force increments ddQ_j
    Vector mElementForceIncr = ZeroVector(msLocalSize);  // element force increments dQ_j
    Vecotr mElementResistingForces = ZeroVector(msLocalSize); // element resisting forces

    /**
     * Element state determination process consists of the following functions
     * 
     */ 

    /**
     * (4) Compute the element deformation increments
     * 
     * using the compatibility matrix L_ele, the element deformation increments ddq_i
     * is computed from the structure displacement increments ddp_i
     */
    void ComputeElementDeformationIncrements();


    /**
     * (5) Start the element state deternimation
     * 
     * Set j = 1
     * 
     */    

    /**
     * (6) Compute the change in the element force increments
     * 
     */    
    void ComputeChangeInElementForceIncr();


    /**
     * (7) Update the element force increments and the element resisting forces
     * 
     */  
    void UpdateElementForceIncr(Vector mElementForceIncr);
    void UpdateElementResistingForces(Vector mElementResistingForces);

    /**
     * (8) Compute the section force increments
     * 
     */
    void ComputeSectionForceIncr();


    




    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
    };

}



#endif