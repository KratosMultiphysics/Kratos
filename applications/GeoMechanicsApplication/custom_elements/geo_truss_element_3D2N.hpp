// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_TRUSS_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_GEO_TRUSS_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "../StructuralMechanicsApplication/custom_elements/truss_element_3D2N.hpp"

namespace Kratos
{
    /**
     * @class GeoTrussElement3D2N
     *
     * @brief This is a 3D-2node truss element with 3 translational dofs per node
     *
     * @author Vahid Galavi
     */

    class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoTrussElement3D2N : public TrussElement3D2N
    {
    protected:
        //const values
        static constexpr int mStressVectorSize = 1;
        Vector mInternalStresses = ZeroVector(mStressVectorSize);
        Vector mInternalStressesFinalized = ZeroVector(mStressVectorSize);
        Vector mInternalStressesFinalizedPrevious = ZeroVector(mStressVectorSize);

    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoTrussElement3D2N);

        GeoTrussElement3D2N() {};
        GeoTrussElement3D2N(IndexType NewId,
                            GeometryType::Pointer pGeometry);
        GeoTrussElement3D2N(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties);


        ~GeoTrussElement3D2N() override;

        /**
         * @brief Creates a new element
         * @param NewId The Id of the new created element
         * @param pGeom The pointer to the geometry of the element
         * @param pProperties The pointer to property
         * @return The pointer to the created element
         */
        Element::Pointer Create( IndexType NewId,
                                 GeometryType::Pointer pGeom,
                                 PropertiesType::Pointer pProperties ) const override;

        /**
         * @brief Creates a new element
         * @param NewId The Id of the new created element
         * @param ThisNodes The array containing nodes
         * @param pProperties The pointer to property
         * @return The pointer to the created element
         */
        Element::Pointer Create( IndexType NewId,
                                 NodesArrayType const& ThisNodes,
                                 PropertiesType::Pointer pProperties ) const override;

        void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function updates the internal normal force w.r.t. the current deformations
         * @param rinternalForces The current updated internal forces
         */
        void UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces,
                                  const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                          std::vector<Vector>& rOutput,
                                          const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                          std::vector< array_1d<double, 3 > >& rOutput,
                                          const ProcessInfo& rCurrentProcessInfo) override;


        void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
        void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;


    private:
        bool mIsInitialization = false;


        friend class Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
    };
}


#endif
