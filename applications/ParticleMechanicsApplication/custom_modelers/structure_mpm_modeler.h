//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//


#if !defined(KRATOS_STRUCTURE_MPM_MODELER_H_INCLUDED )
#define  KRATOS_STRUCTURE_MPM_MODELER_H_INCLUDED


// System includes

// External includes

// Project includes
#include "modeler/modeler.h"

#include "geometries/geometry_data.h"
#include "geometries/coupling_geometry.h"
#include "geometries/line_2d_2.h"
#include "geometries/quadrature_point_curve_on_surface_geometry.h"
#include "utilities/quadrature_points_utility.h"
#include "utilities/binbased_fast_point_locator.h"
#include "particle_mechanics_application_variables.h"
#include "custom_utilities/mpm_search_element_utility.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Finds realtion between structure FEM and MPM to make Mortar Mapping.
/* Detail class definition.
*/
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) StructureMpmModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(StructureMpmModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;

    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    StructureMpmModeler()
        : Modeler()
    {
        mpModelOrigin = nullptr;
        mpModelDest = nullptr;
    }

    /// Constructor.
    StructureMpmModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters())
        : Modeler(rModel, ModelerParameters)
    {
        mpModelOrigin = &rModel;
        if (mpModelOrigin->HasModelPart("Background_Grid")) mIsOriginMpm = true;
        else mIsOriginMpm = false;

        mpModelDest = nullptr;
    }

    /// Destructor.
    virtual ~StructureMpmModeler() = default;

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<StructureMpmModeler>(rModel, ModelParameters);
    }

    /// Adds the second model part to the modeler.
    void GenerateNodes(ModelPart& ThisModelPart) override
    {
        mpModelDest = &ThisModelPart.GetModel();
    }

    ///@}
    ///@name Stages
    ///@{

    // Initially setup the relations
    void SetupGeometryModel() override;

    void UpdateGeometryModel();

    // Update relations
    void PrepareGeometryModel() override
    {
        this->UpdateGeometryModel();
    }

    ///@}
    ///@name Model Part generation
    ///@{

    template<
        class TLineGeometriesList,
        class TQuadraturePointGeometriesList = TLineGeometriesList>
    void CreateStructureQuadraturePointGeometries(
        TLineGeometriesList& rInputLineGeometries,
        TQuadraturePointGeometriesList& rOuputQuadraturePointGeometries,
        GeometryData::IntegrationMethod ThisIntegrationMethod
        )
    {
        for (IndexType i = 0; i < rInputLineGeometries.size(); ++i) {
            std::vector<GeometryPointerType> qudrature_point_geometries =
                CreateQuadraturePointsUtility<Node<3>>::Create(rInputLineGeometries[i], ThisIntegrationMethod);
            for (IndexType j = 0; j < qudrature_point_geometries.size(); ++j) {
                qudrature_point_geometries[j]->SetGeometryParent(rInputLineGeometries[i].get());
                rOuputQuadraturePointGeometries.push_back(qudrature_point_geometries[j]);
            }
        }
    }

    void CreateExactlySegmentedStructureQuadraturePointGeometries(
        ModelPart& rInterfaceModelPart,
        std::vector<GeometryPointerType>& rFEMInterfaceGeometries,
        ModelPart& rBackgroundGrid,
        std::vector<GeometryPointerType>& rSegmentedQuadraturePoints);

    template<SizeType TDimension, class TQuadraturePointGeometriesList>
    void CreateMpmQuadraturePointGeometries(
        const TQuadraturePointGeometriesList& rInputQuadraturePointGeometries,
        TQuadraturePointGeometriesList& rOuputQuadraturePointGeometries,
        ModelPart& rBackgroundGridModelPart) {

        if (rOuputQuadraturePointGeometries.size() != rInputQuadraturePointGeometries.size()) {
            rOuputQuadraturePointGeometries.resize(rInputQuadraturePointGeometries.size());
        }

        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
        SearchStructure.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);
        typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

        double fem_quadrature_points_edge_length = 0;
        double mpm_quadrature_points_edge_length = 0;

        const double tolerance = mParameters["minimum_shape_function_value"].GetDouble();
        for (size_t i = 0; i < rInputQuadraturePointGeometries.size(); ++i)
        {
            array_1d<double, 3> coordinates = rInputQuadraturePointGeometries[i]->Center();

            Element::Pointer p_elem;
            Vector N;

            // FindPointOnMesh find the background element in which a given point falls and the relative shape functions
            bool is_found = SearchStructure.FindPointOnMesh(coordinates, N, p_elem, result_begin, 100, 1e-12);

            if (is_found) {
                double integration_weight = rInputQuadraturePointGeometries[i]->IntegrationPoints()[0].Weight();
                //integration_weight *= rInputQuadraturePointGeometries[i]->DeterminantOfJacobian(0);

                array_1d<double, 3> local_coordinates;
                auto& r_geometry = p_elem->GetGeometry();
                r_geometry.PointLocalCoordinates(local_coordinates, coordinates);

                Vector N;
                r_geometry.ShapeFunctionsValues(N, local_coordinates);

                Matrix DN_De;
                r_geometry.ShapeFunctionsLocalGradients(DN_De, local_coordinates);
                Matrix DN_De_non_zero(DN_De.size1(), DN_De.size2());

                typename GeometryType::PointsArrayType points;

                Matrix N_matrix(1, N.size());
                SizeType non_zero_counter = 0;
                for (IndexType i_N = 0; i_N < N.size(); ++i_N) {
                    if (N[i_N] > tolerance)
                    {
                        N_matrix(0, non_zero_counter) = N[i_N];
                        for (IndexType j = 0; j < DN_De.size2(); j++) {
                            DN_De_non_zero(non_zero_counter, j) = DN_De(i_N, j);
                        }
                        points.push_back(r_geometry(i_N));
                        non_zero_counter++;
                    }
                }

                N_matrix.resize(1, non_zero_counter, true);
                DN_De_non_zero.resize(non_zero_counter, DN_De.size2(), true);

                Matrix jacci;
                rInputQuadraturePointGeometries[i]->Jacobian(jacci, 0);
                Vector space_derivatives = column(jacci, 0);

                Matrix inv;
                r_geometry.InverseOfJacobian(inv, local_coordinates);
                KRATOS_DEBUG_ERROR_IF_NOT(inv.size2() == space_derivatives.size()) << "Jacobian and space derivative sizes are mismatched!\n";
                Vector local_tangent = prod(inv, space_derivatives);
                Matrix jacci_the_wacci;
                r_geometry.Jacobian(jacci_the_wacci, local_coordinates);
                double J2 = norm_2(column(jacci_the_wacci, 0) * local_tangent[0] + column(jacci_the_wacci, 1) * local_tangent[1]);
                //integration_weight /= J2;

                IntegrationPoint<3> int_p(local_coordinates, integration_weight);

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    r_geometry.GetDefaultIntegrationMethod(),
                    int_p,
                    N_matrix,
                    DN_De_non_zero);

                rOuputQuadraturePointGeometries[i] = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCurveOnSurface(data_container,
                    points, local_tangent[0], local_tangent[1], p_elem->pGetGeometry().get());

                /*
                #ifdef KRATOS_DEBUG
                std::vector<array_1d<double, 3>> space_derivatives_check(0.0);
                rOuputQuadraturePointGeometries[i]->GlobalSpaceDerivatives(space_derivatives_check, 0, 1);
                Vector tangent_check = space_derivatives_check[1] * local_tangent[0] +
                    space_derivatives_check[2] * local_tangent[1];
                tangent_check.resize(space_derivatives.size(), true);

                if ((norm_2(tangent_check - space_derivatives) > tolerance))
                {

                    KRATOS_WATCH(rInputQuadraturePointGeometries[i]->DeterminantOfJacobian(0));
                    KRATOS_WATCH(rInputQuadraturePointGeometries[i]->DeterminantOfJacobian(0));
                    KRATOS_WATCH(rInputQuadraturePointGeometries[i]->DeterminantOfJacobian(0));
                    KRATOS_WATCH(J2);
                    KRATOS_WATCH(tangent_check);
                    KRATOS_WATCH(space_derivatives);

                    Vector det_jacobian;
                    if (&(rOuputQuadraturePointGeometries[i]->GetGeometryParent(0)) == nullptr)
                    {
                        rOuputQuadraturePointGeometries[i]->DeterminantOfJacobian(det_jacobian);
                    }
                    else
                    {
                        rOuputQuadraturePointGeometries[i]->Calculate(DETERMINANTS_OF_JACOBIAN_PARENT, det_jacobian);
                    }
                    KRATOS_WATCH(det_jacobian);
                    KRATOS_WATCH(rOuputQuadraturePointGeometries[i]->GetGeometryParent(0));
                        rOuputQuadraturePointGeometries[i]->DeterminantOfJacobian(det_jacobian);
                    KRATOS_WATCH(det_jacobian);


                    int terwadf = 1;

                }
                //KRATOS_ERROR_IF(norm_2(tangent_check - space_derivatives) > tolerance)
                //    << "CreateMpmQuadraturePointGeometries | Line and quadrature point tangents not equal."
                //    << "\nFEM boundary line tangent = " << space_derivatives
                //    << "\nMPM quad point on curve tangent = " << tangent_check << "\n";
                #endif
                */

                fem_quadrature_points_edge_length += rInputQuadraturePointGeometries[i]->IntegrationPoints()[0].Weight() *
                    rInputQuadraturePointGeometries[i]->DeterminantOfJacobian(0);
                Vector det_jacobian;
                if (&(rOuputQuadraturePointGeometries[i]->GetGeometryParent(0)) == nullptr) rOuputQuadraturePointGeometries[i]->DeterminantOfJacobian(det_jacobian);
                else  rOuputQuadraturePointGeometries[i]->Calculate(DETERMINANTS_OF_JACOBIAN_PARENT, det_jacobian);
                mpm_quadrature_points_edge_length += rOuputQuadraturePointGeometries[i]->IntegrationPoints()[0].Weight() * det_jacobian[0];
            }
        }

        KRATOS_CHECK_NEAR(fem_quadrature_points_edge_length, mpm_quadrature_points_edge_length, 1e-9);
    }

    template<SizeType TDimension,
        class TConditionsList>
    void UpdateMpmQuadraturePointGeometries(
        TConditionsList& rInputConditions,
        ModelPart& rBackgroundGridModelPart) {

        const IndexType mpm_index = (mIsOriginMpm) ? 0 : 1;
        const IndexType fem_index = 1- mpm_index;

        auto cond_begin = rInputConditions.ptr_begin();
        std::vector< Geometry<Node<3>>*> bin_search_quad_geoms;
        const double tolerance = mParameters["minimum_shape_function_value"].GetDouble();

        bool use_neighbour_search = false;
        if (use_neighbour_search)
        {
            IndexPartition<>(rInputConditions.size()).for_each([&](SizeType i)
                {
                    //for (size_t i = 0; i < rInputConditions.size(); i++)
                    //{

                    Geometry<Node<3>>* p_quad_geom = (&((*(cond_begin + i))->GetGeometry()));
                    array_1d<double, 3> coordinates = p_quad_geom->GetGeometryPart(fem_index).Center();
                    Element::Pointer p_elem;
                    array_1d<double, 3> local_coordinates;

                    // Try neighbour search first
                    bool is_found = false;
                    GeometryType& r_found_geom = MPMSearchElementUtility::FindGridGeom(p_quad_geom->GetGeometryPart(mpm_index).GetGeometryParent(0),
                        rBackgroundGridModelPart, tolerance, coordinates, local_coordinates,
                        rBackgroundGridModelPart.GetProcessInfo(), is_found);

                    if (is_found) {
                        CreateQuadraturePointsUtility<NodeType>::UpdateFromLocalCoordinates(
                            p_quad_geom->pGetGeometryPart(mpm_index),
                            local_coordinates,
                            p_quad_geom->IntegrationPoints()[0].Weight(),
                            r_found_geom);
                    }
                    else {
                        // Add quad geom to do slower bin search later
                        #pragma omp critical
                        bin_search_quad_geoms.push_back(p_quad_geom);
                    }

                    //}
                }
            );
        }
        else
        {
            for (size_t i = 0; i < rInputConditions.size(); ++i)
            {
                Geometry<Node<3>>* p_quad_geom = (&((*(cond_begin + i))->GetGeometry()));
                bin_search_quad_geoms.push_back(p_quad_geom);
            }
        }


        // Do slow search of remaining quad points
        if (bin_search_quad_geoms.size() > 0)
        {
            Vector N;
            array_1d<double, 3> local_coordinates;
            Element::Pointer p_elem;
            BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();
            typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);
            typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

            for (size_t i = 0; i < bin_search_quad_geoms.size(); ++i)
            {
                array_1d<double, 3> coordinates = bin_search_quad_geoms[i]->GetGeometryPart(fem_index).Center();
                bool is_found = SearchStructure.FindPointOnMesh(coordinates, N, p_elem, result_begin, 100, tolerance);

                if (is_found)
                {
                    GeometryType& r_found_geom = p_elem->GetGeometry();
                    r_found_geom.PointLocalCoordinates(local_coordinates, coordinates);

                    for (size_t geom_node = 0; geom_node < r_found_geom.PointsNumber(); ++geom_node)
                    {
                        if (!r_found_geom(geom_node).get()->Is(ACTIVE))
                        {
                            //KRATOS_INFO("modeler update") << "Mapping to inactive mpm grid node!\n";
                            std::cout << "\n\nInactive grid node:"
                                << "\n\tX = " << r_found_geom(geom_node).get()->X()
                                << "\n\tY = " << r_found_geom(geom_node).get()->Y()
                                << "\n";
                        }
                    }

                    CreateQuadraturePointsUtility<NodeType>::UpdateFromLocalCoordinates(
                        bin_search_quad_geoms[i]->pGetGeometryPart(mpm_index),
                        local_coordinates,
                        bin_search_quad_geoms[i]->IntegrationPoints()[0].Weight(),
                        r_found_geom);
                }
                else KRATOS_ERROR << "Coupling quadrature point moved outside MPM background grid!";
            }
        }
    }


    template<SizeType TDimension,
        class TConditionsList>
        void UpdateMpmQuadraturePointGeometriesWithFilter(
            TConditionsList& rInputConditions,
            ModelPart& rBackgroundGridModelPart) {

        std::cout << "\n\n ---------- UpdateMpmQuadraturePointGeometriesWithFilter ------------- \n\n";

        const IndexType mpm_index = (mIsOriginMpm) ? 0 : 1;
        const IndexType fem_index = 1 - mpm_index;

        auto cond_begin = rInputConditions.ptr_begin();
        std::vector< Geometry<Node<3>>*> bin_search_quad_geoms;
        const double tolerance = mParameters["minimum_shape_function_value"].GetDouble();

        // Do slow bin search
        for (size_t i = 0; i < rInputConditions.size(); ++i)
        {
            Geometry<Node<3>>* p_quad_geom = (&((*(cond_begin + i))->GetGeometry()));
            bin_search_quad_geoms.push_back(p_quad_geom);
        }

        // Do slow search of remaining quad points
        if (bin_search_quad_geoms.size() > 0)
        {
            Vector N;
            Matrix DN_De;
            array_1d<double, 3> local_coordinates;
            Element::Pointer p_elem;
            BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundGridModelPart);
            SearchStructure.UpdateSearchDatabase();
            typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);
            typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();

            for (size_t i = 0; i < bin_search_quad_geoms.size(); ++i)
            {
                array_1d<double, 3> coordinates = bin_search_quad_geoms[i]->GetGeometryPart(fem_index).Center();
                bool is_found = SearchStructure.FindPointOnMesh(coordinates, N, p_elem, result_begin, 100, tolerance);
                auto& quad_fem = bin_search_quad_geoms[i]->GetGeometryPart(fem_index);

                if (is_found)
                {
                    auto& r_geometry = p_elem->GetGeometry();
                    r_geometry.PointLocalCoordinates(local_coordinates, coordinates);
                    r_geometry.ShapeFunctionsValues(N, local_coordinates);
                    r_geometry.ShapeFunctionsLocalGradients(DN_De, local_coordinates);

                    double integration_weight = quad_fem.IntegrationPoints()[0].Weight();
                    integration_weight *= quad_fem.DeterminantOfJacobian(0);

                    Matrix DN_De_non_zero(DN_De.size1(), DN_De.size2());

                    typename GeometryType::PointsArrayType points;

                    Matrix N_matrix(1, N.size());
                    SizeType non_zero_counter = 0;
                    for (IndexType i_N = 0; i_N < N.size(); ++i_N) {
                        if (N[i_N] > tolerance && r_geometry(i_N)->Is(ACTIVE))
                        {
                            N_matrix(0, non_zero_counter) = N[i_N];
                            for (IndexType j = 0; j < DN_De.size2(); j++) {
                                DN_De_non_zero(non_zero_counter, j) = DN_De(i_N, j);
                            }
                            points.push_back(r_geometry(i_N));
                            non_zero_counter++;
                        }
                    }

                    N_matrix.resize(1, non_zero_counter, true);
                    DN_De_non_zero.resize(non_zero_counter, DN_De.size2(), true);

                    Matrix jacci;
                    quad_fem.Jacobian(jacci, 0);
                    Vector space_derivatives = column(jacci, 0);

                    Matrix inv;
                    r_geometry.InverseOfJacobian(inv, local_coordinates);
                    KRATOS_DEBUG_ERROR_IF_NOT(inv.size2() == space_derivatives.size()) << "Jacobian and space derivative sizes are mismatched!\n";
                    Vector local_tangent = prod(inv, space_derivatives); // error
                    Matrix jacci_the_wacci;
                    r_geometry.Jacobian(jacci_the_wacci, local_coordinates);
                    double J2 = norm_2(column(jacci_the_wacci, 0) * local_tangent[0] + column(jacci_the_wacci, 1) * local_tangent[1]);
                    integration_weight /= J2;

                    IntegrationPoint<3> int_p(local_coordinates, integration_weight);

                    GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                        r_geometry.GetDefaultIntegrationMethod(),
                        int_p,
                        N_matrix,
                        DN_De_non_zero);

                    GeometryPointerType p_updated_mpm_quad = CreateQuadraturePointsUtility<NodeType>::CreateQuadraturePointCurveOnSurface(data_container,
                        points, local_tangent[0], local_tangent[1], p_elem->pGetGeometry().get());

                    bin_search_quad_geoms[i]->SetGeometryPart(mpm_index, p_updated_mpm_quad);
                }
                else KRATOS_ERROR << "Coupling quadrature point moved outside MPM background grid!";
            }
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "StructureMpmModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream & rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream & rOStream) const override
    {
    }

    ///@}

private:
    Model* mpModelOrigin;
    Model* mpModelDest;
    bool mIsOriginMpm;

    void CopySubModelPart(ModelPart& rDestinationMP, ModelPart& rReferenceMP)
    {
        rDestinationMP.SetNodes(rReferenceMP.pNodes());
        ModelPart& coupling_conditions = rReferenceMP.GetSubModelPart("coupling_conditions");
        rDestinationMP.SetConditions(coupling_conditions.pConditions());
    }

    void CreateInterfaceLineCouplingConditions(ModelPart& rInterfaceModelPart,
        std::vector<GeometryPointerType>& rGeometries,
        ModelPart& rBackgroundGrid);

    void CheckParameters();

    void FixMPMDestInterfaceNodes(ModelPart& rMPMDestInterfaceModelPart)
    {
        const bool is_gauss_seidel = mParameters["is_gauss_seidel"].GetBool();
        if (!mIsOriginMpm && is_gauss_seidel)
        {
            block_for_each(rMPMDestInterfaceModelPart.Nodes(), [&](Node<3>& rNode)
                {
                    rNode.Fix(DISPLACEMENT_X);
                    rNode.Fix(DISPLACEMENT_Y);
                    rNode.Fix(DISPLACEMENT_Z);
                }
            );
        }
    }

    void ReleaseMPMDestInterfaceNodes(ModelPart& rMPMDestInterfaceModelPart)
    {
        const bool is_gauss_seidel = mParameters["is_gauss_seidel"].GetBool();
        if (!mIsOriginMpm && is_gauss_seidel)
        {
            block_for_each(rMPMDestInterfaceModelPart.Nodes(), [&](Node<3>& rNode)
                {
                    rNode.Free(DISPLACEMENT_X);
                    rNode.Free(DISPLACEMENT_Y);
                    rNode.Free(DISPLACEMENT_Z);

                }
            );
        }
    }

    Parameters GetModelerDefaultSettings() const
    {
        return Parameters(R"({
            "echo_level"                    : 0,
            "minimum_shape_function_value"  : 1e-9,
            "gauss_integration_order"  : 5,
            "is_gauss_seidel"               : true
        })");
    }

    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

    ///@}
}; // Class StructureMpmModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    StructureMpmModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const StructureMpmModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_STRUCTURE_MPM_MODELER_H_INCLUDED  defined
