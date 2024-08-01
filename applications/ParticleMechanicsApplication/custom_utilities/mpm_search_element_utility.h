//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
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
#include "utilities/quadrature_points_utility.h"
#include "includes/kratos_flags.h"

#include "particle_mechanics_application_variables.h"

#include "geometries/geometry.h"
#include "includes/model_part.h"

#include "pqmpm_partition_utilities.h"

namespace Kratos
{
namespace MPMSearchElementUtility
{
    // Standard types
    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    typedef Node<3> NodeType;
    typedef typename ModelPart::GeometryType GeometryType;

    inline double CrossProductDet2D(array_1d<double, 3> VectorA, array_1d<double, 3> VectorB)
    {
        return (VectorA[0] * VectorB[1] - VectorB[0] * VectorA[1]);
    }

    inline bool CheckIsInside(const GeometryType& rGeom, array_1d<double, 3>& LocalCoords, const array_1d<double, 3>& Coords, const double Tolerance, const bool IsCalcLocalCoords = true)
    {

        // TODO some optimisation for simple 2D shapes.

        return rGeom.IsInside(Coords, LocalCoords, Tolerance);
    }

    inline void ConstructNeighbourRelations(GeometryType& rGeom, const ModelPart& rBackgroundGridModelPart)
    {
        std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours;
        std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours_edge_aligned;
        std::vector<typename Geometry<Node<3>>::Pointer> geometry_neighbours_surface_aligned;
        for (IndexType j = 0; j < rBackgroundGridModelPart.NumberOfElements(); j++)
        {
            auto p_geometry_neighbour = (rBackgroundGridModelPart.ElementsBegin() + j)->pGetGeometry();
            if (p_geometry_neighbour->Id() != rGeom.Id()) // dont add the parent as its own neighbour
            {
                for (IndexType n = 0; n < p_geometry_neighbour->size(); n++)
                {
                    for (IndexType k = 0; k < rGeom.size(); k++)
                    {
                        if (rGeom[k].Id() == (*p_geometry_neighbour)[n].Id()) {
                            // Prevent duplicate additions
                            bool add_entry = true;
                            for (size_t i = 0; i < geometry_neighbours.size(); i++)
                            {
                                if (geometry_neighbours[i]->Id() == p_geometry_neighbour->Id())
                                {
                                    add_entry = false;
                                    break;
                                }
                            }
                            if (add_entry)
                            {
                                geometry_neighbours.push_back(p_geometry_neighbour);
                            }
                            break;
                        }
                    }
                }
            }

        }

        
            
        for (IndexType j = 0; j < geometry_neighbours.size(); j++)
        {
            int counter =0;
            for (IndexType n = 0; n < (*geometry_neighbours[j]).size(); n++)
            {
                for (IndexType k = 0; k < rGeom.size(); k++)
                {
                    if (rGeom[k].Id() == (*geometry_neighbours[j])[n].Id())
                        counter+=1;

                }
            }
            if (counter==2){
                
                geometry_neighbours_edge_aligned.push_back(geometry_neighbours[j]);
            }
            else if (counter>2)
                geometry_neighbours_surface_aligned.push_back(geometry_neighbours[j]);
        }

 
        #pragma omp critical
        {
            rGeom.SetValue(GEOMETRY_NEIGHBOURS, geometry_neighbours);
            rGeom.SetValue(GEOMETRY_NEIGHBOURS_EDGE_ALIGNED, geometry_neighbours_edge_aligned);
            rGeom.SetValue(GEOMETRY_NEIGHBOURS_SURFACE_ALIGNED, geometry_neighbours_surface_aligned);
        }
        
    }


    inline bool IsExplicitAndNeedsCorrection(GeometryType::Pointer pQuadraturePoint, const ProcessInfo& rProcessInfo)
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                if (pQuadraturePoint->IntegrationPointsNumber() == 1)
                {
                    for (size_t i = 0; i < pQuadraturePoint->ShapeFunctionsValues().size2(); ++i)
                    {
                        if (pQuadraturePoint->ShapeFunctionsValues()(0, i) < std::numeric_limits<double>::epsilon()) return true;
                    }
                }
            }
        }

        return false;
    }

    inline GeometryType& FindGridGeom(GeometryType& rParentGeom,
        const ModelPart& rBackgroundGridModelPart,
        const double Tolerance,
        const array_1d<double, 3>& xg,
        array_1d<double, 3>& rLocalCoords,
        const ProcessInfo& rProcessInfo,
        bool& IsFound)
    {
        IsFound = false;

        if (CheckIsInside(rParentGeom, rLocalCoords, xg, Tolerance)) {
            IsFound = true;
            return rParentGeom;
        }
        else
        {
            if (!rParentGeom.Has(GEOMETRY_NEIGHBOURS))
                ConstructNeighbourRelations(rParentGeom, rBackgroundGridModelPart);

            auto& geometry_neighbours = rParentGeom.GetValue(GEOMETRY_NEIGHBOURS);
            for (IndexType k = 0; k < geometry_neighbours.size(); ++k) {
                if (CheckIsInside(*geometry_neighbours[k], rLocalCoords, xg, Tolerance)) {
                    IsFound = true;
                    return *(geometry_neighbours[k].get());
                }
            }
        }

        return rParentGeom;
    }


    inline void UpdatePartitionedQuadraturePoint(const ModelPart& rBackgroundGridModelPart,
        const array_1d<double, 3>& rCoordinates,
        Element& rMasterMaterialPoint,
        typename GeometryType::Pointer pQuadraturePointGeometry,
        const double Tolerance)
    {
        KRATOS_TRY;

        array_1d<double, 3> local_coords;
        pQuadraturePointGeometry->IsInside(rCoordinates, local_coords, Tolerance);
        PQMPMPartitionUtilities::PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, rCoordinates,
            local_coords, rMasterMaterialPoint, pQuadraturePointGeometry, Tolerance);

        KRATOS_CATCH("");
    }


    inline void NeighbourSearchElements(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        const double Tolerance)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Elements().size()); ++i) {
            auto element_itr = (rMPMModelPart.ElementsBegin() + i);
            array_1d<double, 3> local_coordinates;
            bool is_found = false;
            std::vector<array_1d<double, 3>> xg;
            element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rBackgroundGridModelPart.GetProcessInfo());

            GeometryType& r_found_geom = FindGridGeom(element_itr->GetGeometry().GetGeometryParent(0),
                rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                rMPMModelPart.GetProcessInfo(), is_found);

            if (is_found)
            {
                const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                    ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                if (is_pqmpm)
                {
                    // Updates the quadrature point geometry.
                    (*element_itr).GetGeometry().SetGeometryParent(&r_found_geom);
                    PQMPMPartitionUtilities::PartitionMasterMaterialPointsIntoSubPoints(rBackgroundGridModelPart, xg[0],
                        local_coordinates, *element_itr, element_itr->pGetGeometry(), Tolerance);
                }
                else
                {
                    CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                        element_itr->pGetGeometry(), local_coordinates,
                        element_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);
                }

                if (IsExplicitAndNeedsCorrection(element_itr->pGetGeometry(), rBackgroundGridModelPart.GetProcessInfo()))
                    is_found = false;
                else {
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j)
                        r_found_geom.Points()[j].Set(ACTIVE);
                }
                int& mp_counter = r_found_geom.GetValue(MP_COUNTER) ; 
                    mp_counter +=1;
            }
            if(!is_found)
            {
                #pragma omp critical
                rMissingElements.push_back(&*element_itr);
            }
        }
    }


    //


    inline void NeighbourSearchConditions(const ModelPart& rMPMModelPart,
        const ModelPart& rBackgroundGridModelPart,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const double Tolerance)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rMPMModelPart.Conditions().size()); ++i) {
            auto condition_itr = rMPMModelPart.Conditions().begin() + i;

            std::vector<array_1d<double, 3>> xg;
            condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

            if (xg.size() > 0 && condition_itr->Is(BOUNDARY))
            {
                array_1d<double, 3> local_coordinates;
                bool is_found = false;

                GeometryType& r_found_geom = FindGridGeom(condition_itr->GetGeometry().GetGeometryParent(0),
                    rBackgroundGridModelPart, Tolerance, xg[0], local_coordinates,
                    rMPMModelPart.GetProcessInfo(), is_found);

                if (is_found)
                {
                    CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                        condition_itr->pGetGeometry(), local_coordinates,
                        condition_itr->GetGeometry().IntegrationPoints()[0].Weight(), r_found_geom);
                    
                    int& mpc_counter = r_found_geom.GetValue(MPC_COUNTER) ; 
                    mpc_counter +=1;
                    double& mpc_area_element = r_found_geom.GetValue(MPC_AREA_ELEMENT) ; 
                    mpc_area_element +=condition_itr->GetGeometry().IntegrationPoints()[0].Weight();
                    condition_itr->Set(ACTIVE);
                    condition_itr->Reset(MODIFIED);
                    
                    for (IndexType j = 0; j < r_found_geom.PointsNumber(); ++j)
                        r_found_geom[j].Set(ACTIVE);
                }
                else
                {
                    #pragma omp critical
                    rMissingConditions.push_back(&*condition_itr);
                }
            }
        }
    }


    inline bool IsFixExplicitAndOnElementEdge(const Vector& N, const ProcessInfo& rProcessInfo)
    {
        if (rProcessInfo.Has(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
            if (rProcessInfo.GetValue(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)) {
                // check if MP is exactly on the edge of the element, this gives spurious strains in explicit
                for (SizeType i = 0; i < N.size(); ++i) {
                    if (std::abs(N[i]) < std::numeric_limits<double>::epsilon()) {
                        return true;
                    }
                }
            }
        }

        return false;
    }


    template <std::size_t TDimension>
    void BinBasedSearchElementsAndConditions(ModelPart& rMPMModelPart,
        ModelPart& rBackgroundGridModelPart,
        std::vector<typename Element::Pointer>& rMissingElements,
        std::vector<typename Condition::Pointer>& rMissingConditions,
        const std::size_t MaxNumberOfResults, const double Tolerance)
    {
        const ProcessInfo& r_process_info = rBackgroundGridModelPart.GetProcessInfo();
        bool is_pqmpm = (r_process_info.Has(IS_PQMPM))
            ? r_process_info.GetValue(IS_PQMPM) : false;

        // Search background grid and make element active
        Vector N;
        const int max_result = 1000;

        #pragma omp parallel
        {
            BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();
            typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(max_result);

            // Element search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMissingElements.size()); ++i) {
                auto element_itr = *(rMissingElements.begin() + i);
                std::vector<array_1d<double, 3>> xg;
                element_itr->CalculateOnIntegrationPoints(MP_COORD, xg, rMPMModelPart.GetProcessInfo());
                typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                Element::Pointer pelem;

                // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                if (is_found == true) {
                    if (IsFixExplicitAndOnElementEdge(N, r_process_info) && !is_pqmpm) {
                        // MP is exactly on the edge. Now we give it a little 'nudge'
                        array_1d<double, 3> xg_nudged = array_1d<double, 3>(xg[0]);
                        std::vector<array_1d<double, 3>> mp_vel;
                        element_itr->CalculateOnIntegrationPoints(MP_VELOCITY, mp_vel, rMPMModelPart.GetProcessInfo());
                        xg_nudged += r_process_info[DELTA_TIME] / 1000.0 * mp_vel[0];
                        if (SearchStructure.FindPointOnMesh(xg_nudged, N, pelem, result_begin, MaxNumberOfResults, Tolerance)) {
                            element_itr->SetValuesOnIntegrationPoints(MP_COORD, { xg_nudged }, rMPMModelPart.GetProcessInfo());
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: To prevent spurious explicit stresses, Material Point "
                                << element_itr->Id() << " was nudged." << std::endl;
                        } else {
                            is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);
                            KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Material Point " << element_itr->Id()
                                << " lies exactly on an element edge and may give spurious results." << std::endl;
                        }
                    }
                    pelem->Set(ACTIVE);

                    const bool is_pqmpm = (rBackgroundGridModelPart.GetProcessInfo().Has(IS_PQMPM))
                        ? rBackgroundGridModelPart.GetProcessInfo().GetValue(IS_PQMPM) : false;
                    if (is_pqmpm)
                    {
                        // Updates the quadrature point geometry.
                        (*element_itr).GetGeometry().SetGeometryParent((pelem->pGetGeometry().get()));
                        UpdatePartitionedQuadraturePoint(rBackgroundGridModelPart, xg[0],
                            *element_itr, pelem->pGetGeometry(), Tolerance);
                    }
                    else
                    {
                        auto p_quadrature_point_geometry = element_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        p_quadrature_point_geometry->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());
                    }

                    auto& r_geometry = element_itr->GetGeometry();
                    for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                        r_geometry[j].Set(ACTIVE);

                    int& mp_counter = pelem->GetGeometry().GetValue(MP_COUNTER); 
                    mp_counter +=1;
                }
                else {
                    KRATOS_INFO("MPMSearchElementUtility") << "WARNING: Search Element for Material Point: " << element_itr->Id()
                        << " is failed. Geometry is cleared." << std::endl;

                    element_itr->GetGeometry().clear();
                    element_itr->Reset(ACTIVE);
                    element_itr->Set(TO_ERASE);
                }
            }

            // Condition search and assign background grid
            #pragma omp for
            for (int i = 0; i < static_cast<int>(rMissingConditions.size()); ++i) {
                auto condition_itr = *(rMissingConditions.begin() + i);
                std::vector<array_1d<double, 3>> xg;
                condition_itr->CalculateOnIntegrationPoints(MPC_COORD, xg, rMPMModelPart.GetProcessInfo());

                if (xg.size() > 0) {
                    // Only search for particle based BCs!
                    // Grid BCs are still applied on MP_model_part but we don't want to search for them.
                    typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
                    Element::Pointer pelem;

                    // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
                    bool is_found = SearchStructure.FindPointOnMesh(xg[0], N, pelem, result_begin, MaxNumberOfResults, Tolerance);

                    if (is_found == true) {
                        auto p_quadrature_point_geometry = condition_itr->pGetGeometry();
                        array_1d<double, 3> local_coordinates;
                        p_quadrature_point_geometry->PointLocalCoordinates(local_coordinates, xg[0]);
                        CreateQuadraturePointsUtility<Node<3>>::UpdateFromLocalCoordinates(
                            p_quadrature_point_geometry, local_coordinates,
                            p_quadrature_point_geometry->IntegrationPoints()[0].Weight(), pelem->GetGeometry());
                        
                        auto& r_geometry = condition_itr->GetGeometry();

                        int& mpc_counter = pelem->GetGeometry().GetValue(MPC_COUNTER) ; 
                        mpc_counter +=1;
                        double& mpc_area_element = pelem->GetGeometry().GetValue(MPC_AREA_ELEMENT) ; 
                        mpc_area_element +=p_quadrature_point_geometry->IntegrationPoints()[0].Weight();
                        condition_itr->Set(ACTIVE);
                        condition_itr->Reset(MODIFIED);

                        for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                            r_geometry[j].Set(ACTIVE);
                    } else {
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

    inline void ResetElementsAndNodes(ModelPart& rBackgroundGridModelPart)
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rBackgroundGridModelPart.Elements().size()); ++i) {
            auto element_itr = rBackgroundGridModelPart.Elements().begin() + i;
            auto& r_geometry = element_itr->GetGeometry();
            element_itr->Reset(ACTIVE);

            int& mpc_counter = element_itr->GetValue(MPC_COUNTER) ; 
            mpc_counter =0;
            double& mpc_area_element = element_itr->GetValue(MPC_AREA_ELEMENT) ; 
            mpc_area_element =0.0;
            int& mp_counter = element_itr->GetValue(MP_COUNTER); 
            mp_counter =0;
            

            for (IndexType j = 0; j < r_geometry.PointsNumber(); ++j)
                r_geometry[j].Reset(ACTIVE);

        }
    }

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
    template<std::size_t TDimension>
    void SearchElement(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart, const std::size_t MaxNumberOfResults,
        const double Tolerance)
    {
        ResetElementsAndNodes(rBackgroundGridModelPart);

        std::vector<typename Element::Pointer> missing_elements;
        std::vector<typename Condition::Pointer> missing_conditions;

        NeighbourSearchElements(rMPMModelPart, rBackgroundGridModelPart, missing_elements, Tolerance);
        NeighbourSearchConditions(rMPMModelPart, rBackgroundGridModelPart, missing_conditions, Tolerance);

        if (missing_conditions.size() > 0 || missing_elements.size() > 0)
            BinBasedSearchElementsAndConditions<TDimension>(rMPMModelPart,
                rBackgroundGridModelPart, missing_elements, missing_conditions,
                MaxNumberOfResults, Tolerance);

    }

    void SearchSuperfluousConstraints(ModelPart& rBackgroundGridModelPart, ModelPart& rMPMModelPart, ModelPart& rConditionModelPart){
    
        PointerVector<Condition> mpc_inactive_elements;
        PointerVector<Condition> mpc_active_elements;
        PointerVector<Condition> mpc_modify_elements;
        PointerVector<Condition> mpc_edge_elements;
        PointerVector<Condition> mpc_forced_active_elements;
        PointerVector<Condition> mpc_perspective_forced_active_elements;
        PointerVector<Condition> mpc_perspective_forced_active_elements2;
        PointerVector<Condition> mpc_checked_elements;
        PointerVector<Condition> mpc_not_inactive_elements;
        PointerVector<Condition> mpc_principal_inactive_elements;


        mpc_active_elements.clear();
        mpc_inactive_elements.clear();
        mpc_edge_elements.clear();
        mpc_forced_active_elements.clear();
        mpc_modify_elements.clear();
        mpc_perspective_forced_active_elements.clear();
        mpc_perspective_forced_active_elements2.clear();
        mpc_checked_elements.clear();
        mpc_not_inactive_elements.clear();

        #pragma omp critical
        {
            // create a list of elements which contain at least one boundary particle
            
            for (int i_cond = 0; i_cond < static_cast<int>(rConditionModelPart.Conditions().size()); ++i_cond) {
                auto condition_itr = rConditionModelPart.Conditions().begin() + i_cond;

                bool add = true;
                for(auto i=mpc_active_elements.begin(); i!=mpc_active_elements.end(); ++i)
                {
                    if (i->GetGeometry().GetGeometryParent(0).Id() == condition_itr->GetGeometry().GetGeometryParent(0).Id())
                        add = false;
                }

                // consider only conditions which are connected to the material and have at least one node with mass
                if (add){
                    int counter=0;
                    for (unsigned int i = 0; i < condition_itr->GetGeometry().size(); i++)
                    {
                        if (condition_itr->GetGeometry()[i].FastGetSolutionStepValue(NODAL_MASS, 0)< std::numeric_limits<double>::epsilon())
                            counter+=1;
                    }
                    // UNDO ONLY REQUIRED FOR CONTACT IF THIS IS ACTIVE SET ONLY IF PARTICLE IS INSIDE!!!
                    if(counter>0)
                        add = false;
                    
                    // if(counter== condition_itr->GetGeometry().size())
                    //     add = false;
                }
                
                // add only one condition per background grid element
                if (add){               
                    mpc_active_elements.push_back(*(condition_itr.base()));
                    // construct neighbour relations
                    if (!condition_itr->GetGeometry().GetGeometryParent(0).Has(GEOMETRY_NEIGHBOURS))
                        ConstructNeighbourRelations(condition_itr->GetGeometry().GetGeometryParent(0), rBackgroundGridModelPart);
                }
            }

            // Edge Elements should be constrained twice!
            int count_edge_element =0;
            for(auto i=mpc_active_elements.begin(); i!=mpc_active_elements.end(); ++i)
            {
                if (i->GetGeometry().WorkingSpaceDimension()==2){
                    auto& geometry_neighbours_edge_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_EDGE_ALIGNED);
                    int counter=0;

                    for (IndexType k = 0; k < geometry_neighbours_edge_aligned.size(); ++k) {
                        if (geometry_neighbours_edge_aligned[k]->GetValue(MPC_COUNTER) > 0)
                            counter +=1;     
                    }
                    if (counter ==1){
                        
                        count_edge_element+=1;
                        
                        // Fix just the start of the boundary condition
                        if (count_edge_element & 1)
                        {
                            mpc_edge_elements.push_back(*i.base());
                        }
                                
                            
                    }
                }
                else {
                    auto& geometry_neighbours_surface_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_SURFACE_ALIGNED);
                    int counter=0;

                    for (IndexType k = 0; k < geometry_neighbours_surface_aligned.size(); ++k) {

                        if (geometry_neighbours_surface_aligned[k]->GetValue(MPC_COUNTER) > 0)
                            counter +=1;  
                    }
                    const GeometryData::KratosGeometryType geo_type = i->GetGeometry().GetGeometryParent(0).GetGeometryType();
                    int min_surfaces=2;
                    if (geo_type == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8){
                        min_surfaces = 3;
                    }
                    if (counter < min_surfaces){
                        mpc_edge_elements.push_back(*i.base());       
                    }
                }
                    
            }
            

            // // SMALL CUT MODIFICATIONS
            // // ##############################################################################################################
            // // modify shape function values if elements are small cut parallel to one edge
            // for(auto i=mpc_active_elements.begin(); i!=mpc_active_elements.end(); ++i)
            // {
            //     int counter=0;
            //     for (int i_cond = 0; i_cond < static_cast<int>(rConditionModelPart.Conditions().size()); ++i_cond) {
            //         auto condition_itr = rConditionModelPart.Conditions().begin() + i_cond;
            
            //         if (i->GetGeometry().GetGeometryParent(0).Id() == condition_itr->GetGeometry().GetGeometryParent(0).Id()){
            //             auto rResult = row(condition_itr->GetGeometry().ShapeFunctionsValues(), 0);
            //             for (unsigned int n_nodes = 0; n_nodes < condition_itr->GetGeometry().size(); n_nodes++)
            //             {   
            //                 if (rResult[n_nodes]<0.01){
            //                     counter+=1;
            //                 }
            //             }        
            //         }

            //     }
            //     auto mpc_counter = i->GetGeometry().GetGeometryParent(0).GetValue(MPC_COUNTER);
            //     if (counter >= mpc_counter) {
            //         mpc_modify_elements.push_back(*i.base());
            //         KRATOS_WATCH("SMALL CUT MODIFICATION")
            //     }
                    
            // }

            // // deactivate small cut elements which are cut nearby a vertex
            // for(auto i=mpc_active_elements.begin(); i!=mpc_active_elements.end(); ++i)
            // {
            //     auto mpc_counter = i->GetGeometry().GetGeometryParent(0).GetValue(MPC_COUNTER);
            //     std::vector<double> area(1);
            //     i->CalculateOnIntegrationPoints(MPC_AREA, area, rMPMModelPart.GetProcessInfo());
            //     if (mpc_counter * area[0] < i->GetGeometry().GetGeometryParent(0).Length()/10){
            //         auto rResult = row(i->GetGeometry().ShapeFunctionsValues(), 0);
            //         for (unsigned int n_nodes = 0; n_nodes < i->GetGeometry().size(); n_nodes++)
            //         {   
            //             if (rResult[n_nodes]<0.01){
            //                 mpc_modify_elements.push_back(*i.base());
            //                 // mpc_inactive_elements.push_back(*i.base());
            //                 // auto& geometry_neighbours_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_ALIGNED);
            //                 KRATOS_WATCH("SMALL CUT")
            //                 // for (IndexType k = 0; k < geometry_neighbours_aligned.size(); ++k) {
            //                 //     checked_elements.push_back(geometry_neighbours_aligned[k]->Id());
            //                 // }
            //             }
            //         }
            //     }
                
            // }
            // // ##############################################################################################################
            
            // search elements which should be deactivated
            for(auto i=mpc_active_elements.begin(); i!=mpc_active_elements.end(); ++i){
                const GeometryData::KratosGeometryType geo_type = i->GetGeometry().GetGeometryParent(0).GetGeometryType();

                // ##############################################################################################################
                // 2D triangular elements
                // ##############################################################################################################
                if (geo_type == GeometryData::KratosGeometryType::Kratos_Triangle2D3){
                    const int max_aligned_edge_neighbors = 2;

                    auto& geometry_neighbours_edge_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_EDGE_ALIGNED);
                    int counter=0;

                    for (IndexType k = 0; k < geometry_neighbours_edge_aligned.size(); ++k) {
                        
                        // check if neighbor contains boundary particles
                        bool active_element=false;
                        for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                        {
                            if (geometry_neighbours_edge_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                active_element=true;
                                break;
                            }
                        }

                        if (active_element){
                            bool count = true;
                            // check if they are already deactivated
                            for(auto inactive=mpc_inactive_elements.begin(); inactive!=mpc_inactive_elements.end(); ++inactive)
                            {
                                if (geometry_neighbours_edge_aligned[k]->Id() == inactive->GetGeometry().GetGeometryParent(0).Id()){
                                    count=false;
                                    break;
                                }
                            }
                            // check if they are an edge element
                            for(auto edge=mpc_edge_elements.begin(); edge!=mpc_edge_elements.end(); ++edge)
                            {
                                if (geometry_neighbours_edge_aligned[k]->Id() == edge->GetGeometry().GetGeometryParent(0).Id()){
                                    count=false;
                                    break;
                                }
                            }
                            
                            if (count)
                                counter +=1; 
                            
                        }
                    }
                    // deactivate the element only if neighbouring elements are also constrained 
                    if (counter>=max_aligned_edge_neighbors)
                        mpc_inactive_elements.push_back(*i.base());
                }
                // ##############################################################################################################
                // 2D quadrilateral elements
                // ##############################################################################################################
                else if (geo_type == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4){
                    const int max_aligned_edge_neighbors = 3;

                    auto& geometry_neighbours_edge_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_EDGE_ALIGNED);
                    int counter=0;

                    for (IndexType k = 0; k < geometry_neighbours_edge_aligned.size(); ++k) {
                        
                        bool active_element=false;
                        for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                        {
                            if (geometry_neighbours_edge_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                active_element=true;
                                break;
                            }
                        }

                        if (active_element){
                            bool count = true;
                            // check if they are already deactivated
                            for(auto inactive=mpc_inactive_elements.begin(); inactive!=mpc_inactive_elements.end(); ++inactive)
                            {
                                if (geometry_neighbours_edge_aligned[k]->Id() == inactive->GetGeometry().GetGeometryParent(0).Id()){
                                    count=false;
                                    break;
                                }
                            }
                            // check if they are an edge element
                            for(auto edge=mpc_edge_elements.begin(); edge!=mpc_edge_elements.end(); ++edge)
                            {
                                if (geometry_neighbours_edge_aligned[k]->Id() == edge->GetGeometry().GetGeometryParent(0).Id()){
                                    count=false;
                                    break;
                                }
                            }
                            
                            if (count)
                                counter +=1; 
                                
                        }
                    }
                    // deactivate the element only if neighbouring elements are also constrained 
                    if (counter>=max_aligned_edge_neighbors)
                        mpc_inactive_elements.push_back(*i.base());
                }
                // ##############################################################################################################
                // 3D tetrahedra elements
                // ##############################################################################################################
                else if (geo_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4){
                    const int max_aligned_surface_neighbors = 2;

                    auto& geometry_neighbours_surface_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_SURFACE_ALIGNED);
                    int counter=0;
                    bool check = true;

                    for(auto checked=mpc_checked_elements.begin(); checked!=mpc_checked_elements.end(); ++checked)
                    {
                        if (checked->Id() == i->GetGeometry().GetGeometryParent(0).Id()){
                            check = false;
                            break;
                        }
                    }

                    if (check){
                        mpc_checked_elements.push_back(*i.base());
                        int count_forced_active=0;
                        mpc_perspective_forced_active_elements.clear();

                        for (IndexType k = 0; k < geometry_neighbours_surface_aligned.size(); ++k) {
                            // if neighbors contain boundary particles
                            bool active_element=false;
                            for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                            {
                                if (geometry_neighbours_surface_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                    active_element=true;
                                    break;
                                }
                            }

                            if (active_element){
                                bool count = true;
                                // check if they are already deactivated
                                for(auto inactive=mpc_inactive_elements.begin(); inactive!=mpc_inactive_elements.end(); ++inactive)
                                {
                                    if (geometry_neighbours_surface_aligned[k]->Id() == inactive->GetGeometry().GetGeometryParent(0).Id()){
                                        count=false;
                                        break;
                                    }
                                }
                                // check if they are an edge element
                                for(auto edge=mpc_edge_elements.begin(); edge!=mpc_edge_elements.end(); ++edge)
                                {
                                    if (geometry_neighbours_surface_aligned[k]->Id() == edge->GetGeometry().GetGeometryParent(0).Id()){
                                        count=false;
                                        break;
                                    }
                                }
                                
                                if (count){
                                    counter +=1; 

                                    // count neighbors which are forced active and search neighbors which could be set forced active, too such that the element is sufficiently constrained 
                                    bool deactivate = true;
                                    for(auto forced_active=mpc_forced_active_elements.begin(); forced_active!=mpc_forced_active_elements.end(); ++forced_active)
                                    {
                                        if (geometry_neighbours_surface_aligned[k]->Id() == forced_active->GetGeometry().GetGeometryParent(0).Id()){
                                            deactivate = false;
                                            count_forced_active+=1;
                                            break;
                                        }
                                    }
                                    if (deactivate){
                                        for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                                        {   
                                            if (geometry_neighbours_surface_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                                mpc_perspective_forced_active_elements.push_back(*active.base());
                                                break;
                                            }
                                        }
                                    }
                                }
                                    
                            }
                        }

                        // check if the condition within the element can be deactivated
                        if (counter>=max_aligned_surface_neighbors){
                            mpc_inactive_elements.push_back(*i.base());

                            for(auto perspective_forced_active=mpc_perspective_forced_active_elements.begin(); perspective_forced_active!=mpc_perspective_forced_active_elements.end(); ++perspective_forced_active){
                                if (count_forced_active<max_aligned_surface_neighbors){
                                    mpc_forced_active_elements.push_back(*perspective_forced_active.base());
                                    count_forced_active+=1;
                                    mpc_checked_elements.push_back(*perspective_forced_active.base());
                                }
                                
                            }
                        }
                        else{
                            mpc_forced_active_elements.push_back(*i.base());
                        }

                        
                    }
                }
                
                // ##############################################################################################################
                // 3D hexahedra elements
                // ##############################################################################################################
                else if (geo_type == GeometryData::KratosGeometryType::Kratos_Hexahedra3D8){
                    const int max_aligned_surface_neighbors = 5;
                    

                    auto& geometry_neighbours_surface_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_SURFACE_ALIGNED);
                    // auto& geometry_neighbours_edge_aligned = i->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_EDGE_ALIGNED);
                    int counter=0;
                    bool check = true;

                    for(auto checked=mpc_checked_elements.begin(); checked!=mpc_checked_elements.end(); ++checked)
                    {
                        if (checked->Id() == i->GetGeometry().GetGeometryParent(0).Id()){
                            check = false;
                            break;
                        }
                    }

                    if (check){
                        mpc_checked_elements.push_back(*i.base());
                        int count_forced_active=0;
                        mpc_perspective_forced_active_elements.clear();

                        for (IndexType k = 0; k < geometry_neighbours_surface_aligned.size(); ++k) {
                            // if neighbors contain boundary particles
                            bool active_element=false;
                            for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                            {
                                if (geometry_neighbours_surface_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                    active_element=true;
                                    break;
                                }
                            }

                            if (active_element){
                                bool count = true;
                                // check if they are already deactivated
                                for(auto inactive=mpc_inactive_elements.begin(); inactive!=mpc_inactive_elements.end(); ++inactive)
                                {
                                    if (geometry_neighbours_surface_aligned[k]->Id() == inactive->GetGeometry().GetGeometryParent(0).Id()){
                                        count=false;
                                        break;
                                    }
                                }
                                // check if they are an edge element
                                for(auto edge=mpc_edge_elements.begin(); edge!=mpc_edge_elements.end(); ++edge)
                                {
                                    if (geometry_neighbours_surface_aligned[k]->Id() == edge->GetGeometry().GetGeometryParent(0).Id()){
                                        count=false;
                                        break;
                                    }
                                }
                                
                                if (count){
                                    counter +=1; 

                                    // count neighbors which are forced active and search neighbors which could be set forced active, too such that the element is sufficiently constrained 
                                    bool deactivate = true;
                                    for(auto forced_active=mpc_forced_active_elements.begin(); forced_active!=mpc_forced_active_elements.end(); ++forced_active)
                                    {
                                        if (geometry_neighbours_surface_aligned[k]->Id() == forced_active->GetGeometry().GetGeometryParent(0).Id()){
                                            deactivate = false;
                                            count_forced_active+=1;
                                            break;
                                        }
                                    }
                                    if (deactivate){
                                        for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                                        {   
                                            if (geometry_neighbours_surface_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                                mpc_perspective_forced_active_elements.push_back(*active.base());
                                                break;
                                            }
                                        }
                                    }
                                }
                                    
                            }
                        }

                        // check if the condition within the element can be deactivated
                        if (counter>=max_aligned_surface_neighbors){
                            mpc_inactive_elements.push_back(*i.base());

                            // check if the neighbors have enough active neighbors such that they could be deactivated
                            mpc_not_inactive_elements.clear();
                            mpc_principal_inactive_elements.clear();
                            for(auto perspective_forced_active=mpc_perspective_forced_active_elements.begin(); perspective_forced_active!=mpc_perspective_forced_active_elements.end(); ++perspective_forced_active){
                                auto& geometry_neighbours_surface_aligned_of_neighbors = perspective_forced_active->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_SURFACE_ALIGNED);

                                int counter2 = 0;
                                for (IndexType m = 0; m < geometry_neighbours_surface_aligned_of_neighbors.size(); ++m) {
                                    bool active_element=false;
                                    for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                                    {
                                        if (geometry_neighbours_surface_aligned_of_neighbors[m]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                            active_element=true;
                                            break;
                                        }
                                    }

                                    if (active_element){
                                        bool count = true;
                                        // check if they are already deactivated
                                        for(auto inactive=mpc_inactive_elements.begin(); inactive!=mpc_inactive_elements.end(); ++inactive)
                                        {
                                            if (geometry_neighbours_surface_aligned_of_neighbors[m]->Id() == inactive->GetGeometry().GetGeometryParent(0).Id()){
                                                count=false;
                                                break;
                                            }
                                        }
                                        if (count)
                                            counter2 +=1;
                                    }
                                }
                                if (counter2<max_aligned_surface_neighbors){
                                    mpc_not_inactive_elements.push_back(*perspective_forced_active.base());
                                }
                                else{
                                    mpc_principal_inactive_elements.push_back(*perspective_forced_active.base());
                                }
                            }

                            // force active those elements which could not be deleted
                            for(auto principal_forced_active=mpc_not_inactive_elements.begin(); principal_forced_active!=mpc_not_inactive_elements.end(); ++principal_forced_active){
                                if (count_forced_active<max_aligned_surface_neighbors){
                                    mpc_forced_active_elements.push_back(*principal_forced_active.base());
                                    count_forced_active+=1;
                                    mpc_checked_elements.push_back(*principal_forced_active.base());
                                }
                            }

                            for(auto perspective_forced_active=mpc_principal_inactive_elements.begin(); perspective_forced_active!=mpc_principal_inactive_elements.end(); ++perspective_forced_active){
                                if (count_forced_active<max_aligned_surface_neighbors){
                                    mpc_forced_active_elements.push_back(*perspective_forced_active.base());
                                    count_forced_active+=1;
                                    mpc_checked_elements.push_back(*perspective_forced_active.base());
                                }
                                else{
                                    mpc_inactive_elements.push_back(*perspective_forced_active.base());
                                    mpc_checked_elements.push_back(*perspective_forced_active.base());

                                    // for constraints which are set inactive sufficient forced active elements are needed
                                    auto& inactive_geometry_neighbours_surface_aligned = perspective_forced_active->GetGeometry().GetGeometryParent(0).GetValue(GEOMETRY_NEIGHBOURS_SURFACE_ALIGNED);
                                    mpc_perspective_forced_active_elements2.clear();
                                    int count_forced_active2 = 0;
                                    int counter2 = 0;
                                    for (IndexType n = 0; n < inactive_geometry_neighbours_surface_aligned.size(); ++n) {
                                        bool active_element=false;
                                        for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                                        {
                                            if (inactive_geometry_neighbours_surface_aligned[n]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                                active_element=true;
                                                break;
                                            }
                                        }

                                        if (active_element){
                                            bool count = true;
                                            // check if they are already deactivated
                                            for(auto inactive=mpc_inactive_elements.begin(); inactive!=mpc_inactive_elements.end(); ++inactive)
                                            {
                                                if (inactive_geometry_neighbours_surface_aligned[n]->Id() == inactive->GetGeometry().GetGeometryParent(0).Id()){
                                                    count=false;
                                                    break;
                                                }
                                            }
                                            
                                            if (count){
                                                // count neighbors which are forced active and search neighbors which could be set forced active, too such that the element is sufficiently constrained 
                                                bool is_forced_active = false;
                                                for(auto forced_active=mpc_forced_active_elements.begin(); forced_active!=mpc_forced_active_elements.end(); ++forced_active)
                                                {
                                                    if (inactive_geometry_neighbours_surface_aligned[n]->Id() == forced_active->GetGeometry().GetGeometryParent(0).Id()){
                                                        is_forced_active = true;
                                                        count_forced_active2+=1;
                                                        break;
                                                    }
                                                }
                                                if (is_forced_active == false){
                                                    for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                                                    {   
                                                        if (inactive_geometry_neighbours_surface_aligned[n]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                                                            mpc_perspective_forced_active_elements2.push_back(*active.base());
                                                            break;
                                                        }
                                                    }
                                                }

                                            } 
                                        }
                                    }
                                    
                                    for(auto perspective_forced_active2=mpc_perspective_forced_active_elements2.begin(); perspective_forced_active2!=mpc_perspective_forced_active_elements2.end(); ++perspective_forced_active2){
                                        if (count_forced_active2<max_aligned_surface_neighbors){
                                            mpc_forced_active_elements.push_back(*perspective_forced_active2.base());
                                            count_forced_active2+=1;
                                            mpc_checked_elements.push_back(*perspective_forced_active2.base());
                                        }
                                        else{
                                            KRATOS_WATCH("I AM HERE")
                                        }
                                    }

                                }
                            }

                            // // deactivate all edge aligned active neighbors
                            // for (IndexType k = 0; k < geometry_neighbours_edge_aligned.size(); ++k) 
                            // {
                            //     bool deactivate = true;
                            //     for(auto checked=mpc_checked_elements.begin(); checked!=mpc_checked_elements.end(); ++checked)
                            //     {
                            //         if (geometry_neighbours_edge_aligned[k]->Id() == checked->GetGeometry().GetGeometryParent(0).Id())
                            //             deactivate=false;
                            //     }  
                            //     if (deactivate){
                            //         for(auto active=mpc_active_elements.begin(); active!=mpc_active_elements.end(); ++active)
                            //         {   
                            //             if (geometry_neighbours_edge_aligned[k]->Id() == active->GetGeometry().GetGeometryParent(0).Id()){
                            //                 mpc_inactive_elements.push_back(*active.base());
                            //                 mpc_checked_elements.push_back(*active.base());
                            //             }
                            //         }
                            //     } 

                            // }

                            

                        }
                        else{
                            mpc_forced_active_elements.push_back(*i.base());
                        }

                        
                    }
                } //Choose geometry


            } //loop conditions   
            

            KRATOS_WATCH(mpc_inactive_elements.size())
            KRATOS_WATCH(mpc_edge_elements.size())
            KRATOS_WATCH(mpc_active_elements.size())
            KRATOS_WATCH(mpc_forced_active_elements.size())

            // deactivate conditions which are in the inactive defined elements
            for (int i_cond = 0; i_cond < static_cast<int>(rConditionModelPart.Conditions().size()); ++i_cond) {
                auto condition_itr = rConditionModelPart.Conditions().begin() + i_cond;

                for(auto i=mpc_inactive_elements.begin(); i!=mpc_inactive_elements.end(); ++i)
                {
                    if (i->GetGeometry().GetGeometryParent(0).Id() == condition_itr->GetGeometry().GetGeometryParent(0).Id()){
                        condition_itr->Reset(ACTIVE);  
                    }
                        
                }

                // shape function values should be modified for elements cut parallel to an edge
                for(auto i=mpc_modify_elements.begin(); i!=mpc_modify_elements.end(); ++i)
                {
                    if (i->GetGeometry().GetGeometryParent(0).Id() == condition_itr->GetGeometry().GetGeometryParent(0).Id())
                        condition_itr->Set(MODIFIED);
                }
                
            }
        
        }
    }



} // end namespace MPMSearchElementUtility

} // end namespace Kratos

#endif // KRATOS_MPM_SEARCH_ELEMENT_UTILITY

