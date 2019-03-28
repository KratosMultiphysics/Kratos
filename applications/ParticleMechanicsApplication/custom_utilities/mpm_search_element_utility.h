//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_SEARCH_ELEMENT_UTILITY
#define KRATOS_MPM_SEARCH_ELEMENT_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/binbased_fast_point_locator.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
namespace MPMSearchElementUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    /**
     * @brief Search element connectivity for each particle
     * @details A search is performed to know in which grid element the material point falls.
     * If one or more material points fall in the grid element, the grid element is
     * set to be active and its connectivity is associated to the material point
     * element.
     * STEPS:
     * 1) All the elements are set to be INACTIVE
     * 2) A searching is performed and the grid elements which contain at least a MP are set to be ACTIVE
     *
     */
    template<std::size_t TDim>
    void SearchElement(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart, const std::size_t MaxNumberOfResults,
        const double Tolerance)
    {
        // Reset elements to inactive
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i){
                auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
                auto& rGeom = element_itr->GetGeometry();
                element_itr->Reset(ACTIVE);

                for (IndexType j=0; j < rGeom.PointsNumber(); ++j)
                    rGeom[j].Reset(ACTIVE);

        }

        // Search background grid and make element active
        Vector N;
        const int max_result = 1000;

        #pragma omp parallel
        {
            BinBasedFastPointLocator<TDim> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();

            typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_result);

            // Element search and assign background grid
            #pragma omp for
            for(int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i){

                auto element_itr = rMPMModelPart.Elements().begin() + i;

                const array_1d<double,3>& xg = element_itr->GetValue(MP_COORD);
                typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg, N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                if (is_found == true) {
                        pelem->Set(ACTIVE);
                        element_itr->GetGeometry() = pelem->GetGeometry();
                        auto& rGeom = element_itr->GetGeometry();

                        for (IndexType j=0; j < rGeom.PointsNumber(); ++j)
                            rGeom[j].Set(ACTIVE);
                }
                else{
                        KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: " << element_itr->Id()
                            << " is failed. Geometry is cleared." << std::endl;

                        element_itr->GetGeometry().clear();
                        element_itr->Reset(ACTIVE);
                        element_itr->Set(TO_ERASE);
                }
            }

            // Condition search and assign background grid
            #pragma omp for
            for(int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i){

                auto condition_itr = rMPMModelPart.Conditions().begin() + i;

                if (condition_itr->Has(MPC_COORD)){
                    const array_1d<double,3>& xg = condition_itr->GetValue(MPC_COORD);
                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();

                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg, N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                            pelem->Set(ACTIVE);
                            condition_itr->GetGeometry() = pelem->GetGeometry();
                            auto& rGeom = condition_itr->GetGeometry();

                            for (IndexType j=0; j < rGeom.PointsNumber(); ++j)
                                rGeom[j].Set(ACTIVE);
                    }
                    else{
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point Condition: " << condition_itr->Id()
                                << " is failed. Geometry is cleared." << std::endl;

                            condition_itr->GetGeometry().clear();
                            condition_itr->Reset(ACTIVE);
                            condition_itr->Set(TO_ERASE);
                    }
                }
            }
        }
    }

} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY


