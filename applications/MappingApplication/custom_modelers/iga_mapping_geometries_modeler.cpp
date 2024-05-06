//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//                   Andrea Gorgi
//

// Project includes
#include "iga_mapping_geometries_modeler.h"
#include "custom_utilities/iga_mapping_intersection_utilities.h"

#include "geometries/brep_curve_on_surface.h"
#include "geometries/nurbs_surface_geometry.h"
#include "geometries/nurbs_curve_on_surface_geometry.h"
#include "geometries/nurbs_curve_geometry.h"
#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void IgaMappingGeometriesModeler::SetupGeometryModel()
    {
        CheckParameters();

        ModelPart& coupling_model_part = (mpModels[0]->HasModelPart("coupling"))
            ? mpModels[0]->GetModelPart("coupling")
            : mpModels[0]->CreateModelPart("coupling");

        std::string origin_interface_sub_model_part_name;
        std::string destination_interface_sub_model_part_name;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            origin_interface_sub_model_part_name = mParameters["origin_interface_sub_model_part_name"].GetString();
            destination_interface_sub_model_part_name = mParameters["destination_interface_sub_model_part_name"].GetString();
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }

        // create coupling conditions on interface depending on the dimension
        const IndexType dim = 2;
        if (dim == 2)
        {
            CreateInterfaceLineBrepCurveOnSurface(mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));
            CreateInterfaceLineBrepCurveOnSurface(mpModels.back()->GetModelPart(destination_interface_sub_model_part_name));
        }
        else
        {
            KRATOS_ERROR << "Not implemented yet" << std::endl;
        }

        // Transfer everything into the coupling modelpart
        ModelPart& coupling_interface_origin = (coupling_model_part.HasSubModelPart("interface_origin"))
            ? coupling_model_part.GetSubModelPart("interface_origin")
            : coupling_model_part.CreateSubModelPart("interface_origin");
        CopySubModelPart(coupling_interface_origin,
            mpModels[0]->GetModelPart(origin_interface_sub_model_part_name));

        ModelPart& coupling_interface_destination = (coupling_model_part.HasSubModelPart("interface_destination"))
            ? coupling_model_part.GetSubModelPart("interface_destination")
            : coupling_model_part.CreateSubModelPart("interface_destination");
        CopySubModelPart(coupling_interface_destination,
            mpModels[1]->GetModelPart(destination_interface_sub_model_part_name));

        KRATOS_ERROR_IF(coupling_interface_origin.NumberOfConditions() == 0)
            << "Coupling geometries are currently determined by conditions in the coupling sub model parts,"
            << " but there are currently not conditions in the coupling interface origin sub model part. Please specify some."
            << std::endl;
        SizeType working_dim = coupling_interface_origin.ConditionsBegin()->GetGeometry().WorkingSpaceDimension();
        SizeType local_dim = coupling_interface_origin.ConditionsBegin()->GetGeometry().LocalSpaceDimension();

        if (working_dim == 3 && local_dim == 1)
        {
            IgaMappingIntersectionUtilities::IgaCreateBrepCurveOnSurfaceCouplingGeometries(
                coupling_interface_origin,
                coupling_interface_destination,
                coupling_model_part, 1e-6);
            IgaMappingIntersectionUtilities::IgaCreateQuadraturePointsCoupling1DGeometries2D(
                coupling_model_part, 1e-6);
        }
        else
        {
            KRATOS_ERROR << "Creation of coupling quadrature points not yet supported for requested"
                << " working space dimension = " << working_dim << " and local space dimension = "
                << local_dim << std::endl;
        }
    }

    void IgaMappingGeometriesModeler::CheckParameters()
    {
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_model_part_name"))
            << "Missing \"origin_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_model_part_name"))
            << "Missing \"destination_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;

        KRATOS_ERROR_IF_NOT(mParameters.Has("is_interface_sub_model_parts_specified"))
            << "Missing \"is_interface_sub_model_parts_specified\" in MappingGeometriesModeler Parameters." << std::endl;

        if (mParameters["is_interface_sub_model_parts_specified"].GetBool())
        {
            KRATOS_ERROR_IF_NOT(mParameters.Has("origin_interface_sub_model_part_name"))
                << "Missing \"origin_interface_sub_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;

            KRATOS_ERROR_IF_NOT(mParameters.Has("destination_interface_sub_model_part_name"))
                << "Missing \"destination_interface_sub_model_part_name\" in MappingGeometriesModeler Parameters." << std::endl;
        }
    }

    void IgaMappingGeometriesModeler::CreateInterfaceLineBrepCurveOnSurface(ModelPart& rInterfaceModelPart)
    {
        rInterfaceModelPart.CreateSubModelPart("coupling_conditions");
        ModelPart& coupling_conditions = rInterfaceModelPart.GetSubModelPart("coupling_conditions");
        const ModelPart& root_mp = rInterfaceModelPart.GetRootModelPart();

        GeometryPointerType interface_quadrature_point = rInterfaceModelPart.ConditionsBegin()->pGetGeometry(); // You recover a quedrature point in a curve on surface
        auto* interface_brep_curve_on_surface = &interface_quadrature_point->GetGeometryParent(0); //Get the b_rep_curve_on_surface

        // Get the Nurbs curve on surface 
        GeometryPointerType interface_nurbs_curve_on_surface=(interface_quadrature_point->GetGeometryParent(0)).pGetGeometryPart(std::numeric_limits<IndexType>::max() - 2); 
        
        // Downcast it from the base class (geometry) to the derived class (NurbsCurveOnSurfaceGeometry)
        auto interface_nurbs_curve_on_surface_cast = dynamic_pointer_cast<NurbsCurveOnSurfaceGeometry<3, PointerVector<Point>, PointerVector<NodeType>>>(interface_nurbs_curve_on_surface); 

        // Get a pointer to the underlying surface (get a pointer to the base class)
        GeometryPointerType interface_nurbs_surface = (interface_quadrature_point->GetGeometryParent(0)).pGetGeometryPart(GeometryType::BACKGROUND_GEOMETRY_INDEX); 

        // Downcast the pointer to the derived class (NurbsSurfaceGeometry)
        auto interface_nurbs_surface_cast = dynamic_pointer_cast<NurbsSurfaceGeometry<3, PointerVector<NodeType>>>(interface_nurbs_surface); 
        KRATOS_WATCH(interface_nurbs_surface_cast->KnotsU());

        // Get the Nurbs curve in the parameter space 
        auto p_nurbs_curve=interface_nurbs_curve_on_surface_cast->pGetCurve();

        // intersection between curve and surface parameter space in terms of the curve 1D parameter space.
        std::vector<double> curve_and_parameter_space_intersections;
        interface_brep_curve_on_surface->SpansLocalSpace(curve_and_parameter_space_intersections);
        
        // Determine next condition number
        IndexType condition_id = (root_mp.NumberOfConditions() == 0)
            ? 1 : (root_mp.ConditionsEnd() - 1)->Id() + 1;

        for(size_t nurbs_curve_intersection=0;nurbs_curve_intersection<curve_and_parameter_space_intersections.size()-1;nurbs_curve_intersection++){
            NurbsInterval active_range(curve_and_parameter_space_intersections[nurbs_curve_intersection],curve_and_parameter_space_intersections[nurbs_curve_intersection+1]);

            GeometryType::Pointer pBrep_curve_on_surface =Kratos::make_shared< BrepCurveOnSurface<PointerVector<NodeType>, PointerVector<Point>> >(interface_nurbs_surface_cast,p_nurbs_curve,active_range);

            coupling_conditions.CreateNewCondition("BrepCurveOnSurface", condition_id, pBrep_curve_on_surface, root_mp.ElementsBegin()->pGetProperties());
            condition_id += 1;
        }

        
    }
}
