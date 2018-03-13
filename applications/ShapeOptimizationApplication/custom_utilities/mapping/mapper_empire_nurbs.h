//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Altug Emiroglu, 
//                   Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_MAPPER_EMPIRE_NUBS)
#define KRATOS_MAPPER_EMPIRE_NUBS

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

//Application includes
#include "custom_utilities/mapping/mapper_base.h"
#include "custom_utilities/input_output/external_model_import.h"
#include "custom_utilities/input_output/iges_import.h"

// OpenCASCADE includes
#include <gp_Vec.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopExp_Explorer.hxx>

#include <TopOpeBRepBuild_HBuilder.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>

#include <BRepAlgo_Cut.hxx>
#include <BRepAlgo.hxx>

#include <TDF_Data.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelMap.hxx>
#include <TDF_ChildIterator.hxx>
#include <TDF_MapIteratorOfLabelMap.hxx>

#include <TNaming_NamedShape.hxx>
#include <TNaming_Selector.hxx>
#include <TNaming_Tool.hxx>
#include <TNaming_Builder.hxx>
#include <TNaming.hxx>

#include <TopExp.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <GeomConvert.hxx>
#include <Geom2dConvert.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopoDS_Wire.hxx>
#include <ShapeAnalysis.hxx>
// #include <Standard.hxx>
// #include <Standard_DefineHandle.hxx>

// boost includes
#include <boost/algorithm/string.hpp>

// EMPIRE includes
namespace Empire
{
    #include "MapperLib.h"
}


namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A base class for response functions.

class MapperEmpireNURBS: public Mapper
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MapperEmpireNURBS);

    typedef array_1d<double,3> array_3d;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MapperEmpireNURBS( 
        ModelPart& rDesignSurface,
        Parameters& rParameters 
    ) : mrDesignSurface(rDesignSurface), mrParameters(rParameters)
    {
        KRATOS_TRY;

        std::string str; 

        str = "FE_MESH";
        mFEMeshName.clear();
        mFEMeshName.insert(mFEMeshName.end(), str.begin(), str.end());
        mFEMeshName.push_back('\0');

        str = "IGA_MESH";
        mIGAMeshName.clear();
        mIGAMeshName.insert(mFEMeshName.end(), str.begin(), str.end());
        mIGAMeshName.push_back('\0');

        str = "FE_TO_IGA";
        mMapperFEToIGA.clear();
        mMapperFEToIGA.insert(mFEMeshName.end(), str.begin(), str.end());
        mMapperFEToIGA.push_back('\0');

        str = "IGA_TO_FE";
        mMapperIGAToFE.clear();
        mMapperIGAToFE.insert(mFEMeshName.end(), str.begin(), str.end());
        mMapperIGAToFE.push_back('\0');

        Parameters default_params(R"(
            {
                "design_suface_extrusion": 1.0,
                "fe_mesh_domain_size": 2,
                "iga_model_import":
                {
                    "input_type": "PLEASE_SPECIFY_INPUT_MODEL_FILE_TYPE",
                    "input_filename": "PLEASE_SPECIFY_INPUT_FILENAME",
                    "echo_level": 0
                },
                "nurbs_mapping_fe_to_iga":
                {
                    "enforce_consistency": False,
                    "consistency_tolerance": 1e-6,
                    "surface_projection": 1e-2,
                    "numbe_of_initial_guesses": 20,
                    "max_projection_distance_on_different_patches": 1e-3,
                    "newton_raphson":
                    {
                        "global_iterations": 40,
                        "global_tolerance": 1e-5,
                        "boundary_iterations": 0,
                        "boundary_tolerance": 1e-1,
                        "bisection_iterations": 80,
                        "bisection_tolerance": 1e-5
                    },
                    "gauss_integration_points_triangle": 16,
                    "gauss_integration_points_quad": 25,
                    "weak_dirichlet_conditions_curve"
                    {
                        "curves": False,
                        "surfaces": False,
                        "auto_penalty": False,
                        "alpha_primary": 0.0,
                        "alpha_secondary_bending": 0.0,
                        "alpha_secondary_twisting": 0,0
                    },
                    "patch_coupling":
                    {
                        "weak_patch_continuity": False,
                        "auto_penalty": False,
                        "alpha_primary": 0.0,
                        "alpha_secondary_bending": 0.0,
                        "alpha_secondary_twisting": 0.0
                    },
                    "error_computation":
                    {
                        "error_computation": False,
                        "domain_error": False,
                        "interface_error": False,
                        "curve_error": False
                    }
                },
                "nurbs_mapping_iga_to_fe":
                {
                    "enforce_consistency": False,
                    "consistency_tolerance": 1e-6,
                    "surface_projection": 1e-2,
                    "numbe_of_initial_guesses": 20,
                    "max_projection_distance_on_different_patches": 1e-3,
                    "newton_raphson":
                    {
                        "global_iterations": 40,
                        "global_tolerance": 1e-5,
                        "boundary_iterations": 0,
                        "boundary_tolerance": 1e-1,
                        "bisection_iterations": 80,
                        "bisection_tolerance": 1e-5
                    },
                    "gauss_integration_points_triangle": 16,
                    "gauss_integration_points_quad": 25,
                    "weak_dirichlet_conditions_curve"
                    {
                        "curves": False,
                        "surfaces": False,
                        "auto_penalty": False,
                        "alpha_primary": 0.0,
                        "alpha_secondary_bending": 0.0,
                        "alpha_secondary_twisting": 0,0
                    },
                    "patch_coupling":
                    {
                        "weak_patch_continuity": False,
                        "auto_penalty": False,
                        "alpha_primary": 0.0,
                        "alpha_secondary_bending": 0.0,
                        "alpha_secondary_twisting": 0.0
                    },
                    "error_computation":
                    {
                        "error_computation": False,
                        "domain_error": False,
                        "interface_error": False,
                        "curve_error": False
                    }
                }
            })");
        
        rParameters.ValidateAndAssignDefaults(default_params);

        mDesignSufaceExtrusionLength = rParameters["design_suface_extrusion"].GetDouble();
        mrDesignSurface = rDesignSurface;

        mpExternalModel = new IGESExternalModelImport(rParameters["iga_model_import"]);

        //initiating the sequence to create the mapping between IGA mesh and FE mesh
        mpExternalModel->Execute();

        // makes the finite element mesh
        this->MakeFEMesh();

        // makes IGA mesh
        for (TopoDS_Shape& it : mpExternalModel->GetCompoundNURBSList())
            this->MakeIGAMeshFromNURBSSHape(it);
        
        AssignMappingIds();

        GenerateMapper(mMapperFEToIGA);
        GenerateMapper(mMapperIGAToFE);

        for (TopoDS_Shape& current_NURB : mpExternalModel->GetCompoundNURBSList())
            this->DefinePatchContinuityConditionsOnNurbsShape( mIGAMeshName, current_NURB);

        KRATOS_CATCH("");

    }

    /// Destructor.
    virtual ~MapperEmpireNURBS()
    {
        delete mpExternalModel;
    }

    ///@}
    ///@name Operators
    ///@{
  // ==============================================================================
  void MapToDesignSpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInDesignSpace ) override
  {
      boost::timer mapping_time;
      std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to design space..." << std::endl;

      std::vector<double> variables_in_geometry_space;
      std::vector<double> variables_in_design_space;

      variables_in_geometry_space.resize(mrDesignSurface.Nodes().size()*3);

      for(auto& node_i : mrDesignSurface.Nodes())
      {
          int i = node_i.GetValue(MAPPING_ID);
          array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
          variables_in_geometry_space[3*i] = nodal_variable[0];
          variables_in_geometry_space[3*i+1] = nodal_variable[1];
          variables_in_geometry_space[3*i+2] = nodal_variable[2];
      }

      std::vector<TopoDS_Shape> compound_nurbs_list = mpExternalModel->GetCompoundNURBSList();

      GetAllControlPoints(compound_nurbs_list[0], variables_in_design_space);

      Empire::doConsistentMapping(
          mMapperFEToIGA.data(),
          3, 
          variables_in_geometry_space.size(), 
          variables_in_geometry_space.data(), 
          variables_in_design_space.size(), 
          variables_in_design_space.data()
      );

      Empire::doConsistentMapping(
          mMapperIGAToFE.data(),
          3, 
          variables_in_design_space.size(), 
          variables_in_design_space.data(),
          variables_in_geometry_space.size(), 
          variables_in_geometry_space.data()
      );


      for(auto& node_i : mrDesignSurface.Nodes())
      {
          int i = node_i.GetValue(MAPPING_ID);

          Vector node_vector = ZeroVector(3);
          node_vector(0) = variables_in_geometry_space[3*i];
          node_vector(1) = variables_in_geometry_space[3*i+1];
          node_vector(2) = variables_in_geometry_space[3*i+2];
          node_i.FastGetSolutionStepValue(rNodalVariableInDesignSpace) = node_vector;
      }

      std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
  }

  // --------------------------------------------------------------------------
  void MapToGeometrySpace( const Variable<array_3d> &rNodalVariable, const Variable<array_3d> &rNodalVariableInGeometrySpace ) override
  {
      boost::timer mapping_time;
      std::cout << "\n> Starting to map " << rNodalVariable.Name() << " to geometry space..." << std::endl;

      std::vector<double> variables_in_geometry_space;
      std::vector<double> variables_in_design_space;

      variables_in_geometry_space.resize(mrDesignSurface.Nodes().size()*3);

      for(auto& node_i : mrDesignSurface.Nodes())
      {
          int i = node_i.GetValue(MAPPING_ID);
          array_3d& nodal_variable = node_i.FastGetSolutionStepValue(rNodalVariable);
          variables_in_geometry_space[3*i] = nodal_variable[0];
          variables_in_geometry_space[3*i+1] = nodal_variable[1];
          variables_in_geometry_space[3*i+2] = nodal_variable[2];
      }      

      std::vector<TopoDS_Shape> compound_nurbs_list = mpExternalModel->GetCompoundNURBSList();
      
      GetAllControlPoints(compound_nurbs_list[0], variables_in_design_space);

      Empire::doConsistentMapping(
          mMapperFEToIGA.data(),
          3, 
          variables_in_geometry_space.size(), 
          variables_in_geometry_space.data(),
          variables_in_design_space.size(), 
          variables_in_design_space.data()
      );

      Empire::doConsistentMapping(
          mMapperIGAToFE.data(),
          3, 
          variables_in_design_space.size(), 
          variables_in_design_space.data(),
          variables_in_geometry_space.size(), 
          variables_in_geometry_space.data()
      );

      for(auto& node_i : mrDesignSurface.Nodes())
      {
          int i = node_i.GetValue(MAPPING_ID);

          Vector node_vector = ZeroVector(3);
          node_vector(0) = variables_in_geometry_space[3*i];
          node_vector(1) = variables_in_geometry_space[3*i+1];
          node_vector(2) = variables_in_geometry_space[3*i+2];
          node_i.FastGetSolutionStepValue(rNodalVariableInGeometrySpace) = node_vector;
      }

      std::cout << "> Time needed for mapping: " << mapping_time.elapsed() << " s" << std::endl;
  }

  // ==============================================================================    
    ///@}
    ///@name Operations
    ///@{

    ///@}

protected:
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}

private:
    ///@name Member Variables
    ///@{
    ModelPart& mrDesignSurface;
    Parameters& mrParameters;
    
    ExternalModelImport *mpExternalModel;
    
    double mDesignSufaceExtrusionLength;
    int mEchoLevel;
    std::vector<char> mFEMeshName;
    std::vector<char> mIGAMeshName;
    std::vector<char> mMapperFEToIGA;
    std::vector<char> mMapperIGAToFE;

    

    ///@}
    ///@name Private Operators
    ///@{

    void ReadDesignSurfaceFEMesh(std::map<unsigned int, std::vector<double> >& rNodes, std::map<unsigned int, std::vector<unsigned int> >& rElements)
    {
        rNodes.clear();
        rElements.clear();
        if (mrDesignSurface.GetProcessInfo().GetValue(DOMAIN_SIZE)==2)
        {
            // add the existing nodes of the design surface
            unsigned int max_node_id = 0;
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                
                std::vector<double> current_node_cordinates;
                current_node_cordinates.push_back( node_i->X() );
                current_node_cordinates.push_back( node_i->Y() );
                current_node_cordinates.push_back( node_i->Z() );

                rNodes.insert(std::pair<unsigned int,std::vector<double> >(node_i->Id(), current_node_cordinates ) );
                
                if (max_node_id < node_i->Id())
                    max_node_id = node_i->Id();
            }

            //add imaginary nodes for the extrusion of the elements
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                
                std::vector<double> current_node_cordinates;
                current_node_cordinates.push_back( node_i->X() );
                current_node_cordinates.push_back( node_i->Y() );
                current_node_cordinates.push_back( node_i->Z() + mDesignSufaceExtrusionLength);

                rNodes.insert(std::pair<unsigned int,std::vector<double> >(node_i->Id() + max_node_id, current_node_cordinates ) );
            }

            // add the existing elements of the design surface
            unsigned int element_count = 0;
            for (ModelPart::ElementIterator elem_i = mrDesignSurface.ElementsBegin(); elem_i != mrDesignSurface.ElementsEnd(); ++elem_i)
            {
                element_count++;
                std::vector<unsigned int> element_nodes_1;

                ModelPart::ConditionType::GeometryType &elem_geometry = elem_i->GetGeometry();
                const unsigned int numberOfNodes = elem_geometry.size();
                if (numberOfNodes!=2)
                    KRATOS_THROW_ERROR(std::invalid_argument,"Invalid elements found", "Element ID " + std::to_string(elem_i->Id()) + " should only have 2 nodes." );

                //Add the extended node for 3D extrusion
                element_nodes_1.push_back(elem_geometry[1].Id() + max_node_id);
                //Add existing element nodes
                for (unsigned int i = 0; i < 2; i++)
                    element_nodes_1.push_back(elem_geometry[i].Id());
                
                rElements.insert(std::pair<unsigned int, std::vector<unsigned int> >(element_count, element_nodes_1));

                element_count++;
                std::vector<unsigned int> element_nodes_2;
                
                //Add the extended nodes for 3D extrusion
                element_nodes_2.push_back(elem_geometry[1].Id() + max_node_id);
                element_nodes_2.push_back(elem_geometry[0].Id() + max_node_id);
                element_nodes_2.push_back(elem_geometry[0].Id());
                
                rElements.insert(std::pair<unsigned int, std::vector<unsigned int> >(element_count, element_nodes_2));                
            }                       

        }
        else if (mrDesignSurface.GetProcessInfo().GetValue(DOMAIN_SIZE)==3)
        {
            // add the existing nodes of the design surface
            for (ModelPart::NodeIterator node_i = mrDesignSurface.NodesBegin(); node_i != mrDesignSurface.NodesEnd(); ++node_i)
            {
                
                std::vector<double> current_node_cordinates;
                current_node_cordinates.push_back( node_i->X() );
                current_node_cordinates.push_back( node_i->Y() );
                current_node_cordinates.push_back( node_i->Z() );

                rNodes.insert(std::pair<unsigned int,std::vector<double> >(node_i->Id(), current_node_cordinates ) );
            }

            // add the existing elements of the design surface
            for (ModelPart::ElementIterator elem_i = mrDesignSurface.ElementsBegin(); elem_i != mrDesignSurface.ElementsEnd(); ++elem_i)
            {
                std::vector<unsigned int> element_nodes;

                ModelPart::ConditionType::GeometryType &elem_geometry = elem_i->GetGeometry();
                const unsigned int numberOfNodes = elem_geometry.size();
                if (numberOfNodes!=3)
                    KRATOS_THROW_ERROR(std::invalid_argument,"Invalid elements found", "Element ID " + std::to_string(elem_i->Id()) + " should only have 3 nodes." );

                //Add existing element nodes
                for (unsigned int i = 0; i < 3; i++)
                    element_nodes.push_back(elem_geometry[i].Id());
                
                rElements.insert(std::pair<unsigned int, std::vector<unsigned int> >(elem_i->Id(), element_nodes));
            }             
            
        }
        else {

        }
    }        

    void MakeFEMesh()
    {
        
        std::map<unsigned int, std::vector<double> > node_list;
        std::map<unsigned int, std::vector<unsigned int> > element_list;

        ReadDesignSurfaceFEMesh(node_list, element_list);

        const unsigned int number_of_nodes = node_list.size();
        const unsigned int number_of_elements = element_list.size();

        Empire::initFEMesh(mFEMeshName.data(), number_of_nodes, number_of_elements, false);

        std::vector<int> node_id_list;
        std::vector<double> node_values_list;        
        node_id_list.clear();
        node_values_list.clear();

        for(std::map<unsigned int, std::vector<double> >::iterator it =node_list.begin(); it != node_list.end(); it++)
        {
            int key = (int)it->first;
            std::vector<double> value = it->second;
            
            node_id_list.push_back(key);

            for (unsigned int i=0; i<3; i++)
                node_values_list.push_back(value[i]);
        }
        
        Empire::setNodesToFEMesh(mFEMeshName.data(), node_id_list.data(), node_values_list.data());
        
        std::vector<int> number_of_nodes_per_elem;
        std::vector<int> element_node_ids;

        for(std::map<unsigned int, std::vector<unsigned int> >::iterator it = element_list.begin(); it != element_list.end(); it++)
        {
            std::vector<unsigned int> value = it->second;
            number_of_nodes_per_elem.push_back( value.size() );

            for (unsigned int i=0; i < value.size(); i++)
                element_node_ids.push_back(value[i]);
        }        

        Empire::setElementsToFEMesh(mFEMeshName.data(), number_of_nodes_per_elem.data(), element_node_ids.data());
    }

    void MakeIGAMeshFromNURBSSHape(TopoDS_Shape& rNURBSShape)
    {
        KRATOS_TRY;

        Empire::initIGAMesh(mIGAMeshName.data());
        TopExp_Explorer face_explorer(rNURBSShape, TopAbs_FACE);
        int patch_index = 0;
        int unique_id_counter = 0;

        std::vector<std::vector<int> > patch_unique_id_nets;

        while (face_explorer.More())
        {
            patch_index++;
            const TopoDS_Shape current_face = face_explorer.Current();
            face_explorer.Next();

            if ( mEchoLevel > 0)
                std::cout<<" Gathering information of patch index: "<<patch_index-1<<std::endl;
            
            int p_degree, q_degree, number_u_cps, number_v_cps, number_face_wires;
            std::vector<double> u_knots_vector;
            std::vector<double> v_knots_vector;
            std::vector<double> cp_net;
            std::vector<int> unique_id_net;

            GetFaceData( 
                            current_face,
                            unique_id_counter,
                            number_face_wires,
                            p_degree,
                            q_degree,
                            number_u_cps,
                            number_v_cps,
                            u_knots_vector,
                            v_knots_vector,
                            cp_net,
                            unique_id_net
                        );
            
            patch_unique_id_nets.push_back(unique_id_net);
            unique_id_counter += unique_id_net.size();

            if (mEchoLevel>0)
                std::cout<<" adding patch to EMPIRE"<<std::endl;

            Empire::addPatchToIGAMesh(
                mIGAMeshName.data(),
                p_degree,
                u_knots_vector.size(),
                u_knots_vector.data(),
                q_degree,
                v_knots_vector.size(),
                v_knots_vector.data(),
                number_u_cps,
                number_v_cps,
                cp_net.data(),
                unique_id_net.data()
            );

            TopoDS_Wire outer_wire = ShapeAnalysis::OuterWire(TopoDS::Face(current_face));

            TopExp_Explorer face_wire_explorer(current_face, TopAbs_WIRE);
            while (face_wire_explorer.More())
            {
                const TopoDS_Shape current_face_wire = face_wire_explorer.Current();
                face_wire_explorer.Next();

                int number_trimming_curves;
                bool is_inner_curve;

                GetWireData(
                    outer_wire,
                    current_face_wire,
                    is_inner_curve,
                    number_trimming_curves
                );

                if (mEchoLevel>0)
                    std::cout<<" adding trimming loop to EMPIRE"<<std::endl;
                
                Empire::addTrimmingLoopToPatch(
                    mIGAMeshName.data(),
                    patch_index-1,
                    is_inner_curve,
                    number_trimming_curves

                );

                TopExp_Explorer face_wire_edge_explorer(current_face_wire, TopAbs_EDGE);
                while (face_wire_edge_explorer.More())
                {
                    const TopoDS_Shape current_face_wire_edge = face_wire_edge_explorer.Current();
                    face_wire_edge_explorer.Next();

                    int direction, p_degree, number_u_cps;
                    double active_range_begin, active_range_end;
                    std::vector<double> u_knots_vector;
                    std::vector<double> cp_net;
                    
                    GetEdgeDataWithOrientation (
                        current_face_wire_edge,
                        current_face_wire,
                        current_face,
                        direction,
                        p_degree,
                        u_knots_vector,
                        number_u_cps,
                        cp_net,
                        active_range_begin,
                        active_range_end
                    );

                    if (mEchoLevel>0)
                        std::cout<<" adding trimming curve to trimming loop"<<std::endl;
                    
                    Empire::addTrimmingCurveToTrimmingLoop(
                        mIGAMeshName.data(),
                        patch_index-1,
                        direction,
                        p_degree,
                        u_knots_vector.size(),
                        u_knots_vector.data(),
                        number_u_cps,
                        cp_net.data()
                    );
                }                

            }

            Empire::linearizeTrimmingLoops(mIGAMeshName.data(), patch_index-1);
        }

        KRATOS_CATCH("");
    }

    void GetFaceData(
        const TopoDS_Shape& Face,
        int StartID,
        int& rNoOfFaceWires,
        int& rPDegree,
        int& rQDegree,
        int& rUNoCPs,
        int& rVNoCPs,
        std::vector<double>& rUKnotVector,
        std::vector<double>& rVKnotVector,
        std::vector<double>& rCPNet,
        std::vector<int>& rUniqueIDNet
    )
    {
        KRATOS_TRY;
        
        rNoOfFaceWires = 0;
        TopExp_Explorer face_wire_explorer(Face, TopAbs_WIRE);
        while (face_wire_explorer.More())
        {
            face_wire_explorer.Next();
            rNoOfFaceWires++;
        }

        // Convert a topological TopoDS_FACE into Geom_BSplineSurface
        Handle(Geom_Surface) surface_handle = BRep_Tool::Surface(TopoDS::Face(Face));
        Geom_BSplineSurface& bspline_surface = *GeomConvert::SurfaceToBSplineSurface(surface_handle);
          
        // No of Knots
        int number_u_knots = bspline_surface.NbUKnots();
        int number_v_knots = bspline_surface.NbVKnots();

        rUKnotVector.clear();
        for (int i=1; i<number_u_knots+1; i++)
        {
            for (int j=0; j < bspline_surface.UMultiplicity(i); j++)
            {
                rUKnotVector.push_back(bspline_surface.UKnot(i));
            }
        }

        rVKnotVector.clear();
        for (int i=1; i<number_v_knots+1; i++)
        {
            for (int j=0; j < bspline_surface.VMultiplicity(i); j++)
            {
                rUKnotVector.push_back(bspline_surface.VKnot(i));
            }
        }

        if (bspline_surface.IsUPeriodic())
        {
            if (mEchoLevel>0)
                std::cout<<" found periodic U."<<std::endl;
            rUKnotVector.insert(rUKnotVector.begin(), rUKnotVector[0]);
            rUKnotVector.push_back(rUKnotVector[rUKnotVector.size()-1]);
            bspline_surface.SetUNotPeriodic();
        }
        if (bspline_surface.IsVPeriodic())
        {
            if (mEchoLevel>0)
                std::cout<<" found periodic V."<<std::endl;
            rVKnotVector.insert(rVKnotVector.begin(), rVKnotVector[0]);
            rVKnotVector.push_back(rVKnotVector[rVKnotVector.size()-1]);
            bspline_surface.SetVNotPeriodic();
        }
        
        // No of CPs
        rUNoCPs = bspline_surface.NbUPoles();
        rVNoCPs = bspline_surface.NbVPoles();
        
        // CPs
        TColgp_Array2OfPnt occ_cp_net(1, rUNoCPs, 1, rVNoCPs);
        bspline_surface.Poles(occ_cp_net);

        // CP weights
        TColStd_Array2OfReal occ_cp_weight_net(1, rUNoCPs, 1, rVNoCPs);
        bspline_surface.Weights(occ_cp_weight_net);

        rCPNet.clear();
        rUniqueIDNet.clear();

        for (int iVCP=1; iVCP < rVNoCPs+1; iVCP++)
            for (int iUCP=1; iUCP < rUNoCPs+1; iUCP++)
            {
                const gp_Pnt occ_cp = occ_cp_net.Value(iUCP, iVCP);
                double occ_cp_weight = occ_cp_weight_net.Value(iUCP,iVCP);
                rCPNet.push_back(occ_cp.X());
                rCPNet.push_back(occ_cp.Y());
                rCPNet.push_back(occ_cp.Z());
                rCPNet.push_back(occ_cp_weight);
                rUniqueIDNet.push_back(StartID++);
            }
        
        KRATOS_CATCH("");
    }

    void GetWireData(
        const TopoDS_Wire& rOuterWire, 
        const TopoDS_Shape& rCurrentWire,
        bool& rIsInnerWire,
        int& rNoTrimmingCurves
    )
    {
        KRATOS_TRY;

        if (!rCurrentWire.IsEqual(rOuterWire))
            rIsInnerWire = true;
        else
            rIsInnerWire = false;

        rNoTrimmingCurves = 0;
        TopExp_Explorer edge_explorer(rCurrentWire, TopAbs_EDGE);
        while (edge_explorer.More())
        {
            rNoTrimmingCurves++;
            edge_explorer.Next();
        }

        KRATOS_CATCH("");
    }

    void GetEdgeDataWithOrientation(
        const TopoDS_Shape& rEdge,
        const TopoDS_Shape& rWire,
        const TopoDS_Shape& rFace,
        int& rDirection,
        int& rPDegree,
        std::vector<double>& rUKnotVector,
        int& rUNoCPs,
        std::vector<double>& rCPNet,
        double& rActiveRangeBegin,
        double& rAcriveRangeEnd
    )
    {
        KRATOS_TRY;

        // Convert a topological TopoDS_EDGE into Geom2d_BSplineCurve
        Handle(Geom2d_Curve) curve_2d_handle = BRep_Tool::CurveOnSurface(TopoDS::Edge(rEdge), TopoDS::Face(rFace), rActiveRangeBegin, rAcriveRangeEnd);
        Handle(Geom2d_BSplineCurve) bspline_curve_2d = Geom2dConvert::CurveToBSplineCurve(curve_2d_handle);
        //Geom2d_BSplineCurve bspline_curve_2d(Geom2dConvert::CurveToBSplineCurve(curve_2d_handle));
        
        if (bspline_curve_2d->IsPeriodic())
        {
            std::cout<<" Found periodic curve. setting it to non-periodic curve."<<std::endl;
            bspline_curve_2d->SetNotPeriodic();
        }
        
        // Curve direction 
        rDirection = (rWire.Orientation() + rEdge.Orientation() + 1) % 2;
        // pDegree
        rPDegree = bspline_curve_2d->Degree();
        // No of Knots
        int number_of_knots = bspline_curve_2d->NbKnots();
        // Knot vector
        TColStd_Array1OfInteger occ_u_knot_multi(1, number_of_knots);
        bspline_curve_2d->Multiplicities(occ_u_knot_multi);

        int number_u_knots = 0;
        for (int ctr=1; ctr<number_of_knots+1; ctr++)
            number_u_knots += occ_u_knot_multi.Value(ctr);

        TColStd_Array1OfReal occ_u_knot_sequence(1, number_u_knots);
        bspline_curve_2d->KnotSequence(occ_u_knot_sequence);

        rUKnotVector.clear();
        for (int i=1; i<number_u_knots+1; i++)
            rUKnotVector.push_back(occ_u_knot_sequence.Value(i));

        // No of CPs
        rUNoCPs = bspline_curve_2d->NbPoles();
        // CPs
        TColgp_Array1OfPnt2d occ_cp_net(1, rUNoCPs);
        bspline_curve_2d->Poles(occ_cp_net);
        // CP weights
        TColStd_Array1OfReal occ_cp_weight_net(1, rUNoCPs);
        bspline_curve_2d->Weights(occ_cp_weight_net);
        
        rCPNet.clear();
        for (int i=1; i<rUNoCPs+1; i++)
        {
            const gp_Pnt2d occ_cp = occ_cp_net.Value(i);
            double occ_cp_weight = occ_cp_weight_net.Value(i);
            rCPNet.push_back(occ_cp.X());
            rCPNet.push_back(occ_cp.Y());
            rCPNet.push_back(0.0);
            rCPNet.push_back(occ_cp_weight);
        }
        
        KRATOS_CATCH("");
    }

    void GenerateMapper(std::vector<char> rMapperName)
    {
        std::string input_mapper_name(rMapperName.data());

        char* mapper_name=NULL;
        char* mesh_name_a=NULL;
        char* mesh_name_b=NULL;
        std::string settings;
        
        if (input_mapper_name.compare(mMapperIGAToFE.data())==0) 
        {
            mapper_name = mMapperIGAToFE.data();
            mesh_name_a = mIGAMeshName.data();
            mesh_name_b = mFEMeshName.data();
            settings = "nurbs_mapping_iga_to_fe";
        }
        else if (input_mapper_name.compare(mMapperFEToIGA.data())==0)
        {
            mapper_name = mMapperFEToIGA.data();
            mesh_name_a = mFEMeshName.data();
            mesh_name_b = mIGAMeshName.data();
            settings = "nurbs_mapping_fe_to_iga";
        }


        Empire::initIGAMortarMapper(mapper_name, mesh_name_a, mesh_name_b);
        SetMapperParameters(rMapperName, mrParameters[settings]);

        Empire::initialize(mapper_name);
        Empire::buildCouplingMatrices(mapper_name);  

    }    

    void SetMapperParameters(std::vector<char>& rMapperName, Parameters rParameters)
    {
        char *mapper_name = rMapperName.data();

        bool enforceConsistency = rParameters["enforce_consistency"].GetBool();
        double tolConsistency = rParameters["enforce_consistency"].GetDouble();
        double maxProjectionDistance = rParameters["surface_projection"].GetDouble();
        int noInitialGuess = rParameters["numbe_of_initial_guesses"].GetInt();
        double maxProjectionDistanceOnDifferentPatches = rParameters["max_projection_distance_on_different_patches"].GetDouble();

        int NRmaxNumOfIterations = rParameters["newton_raphson"]["global_iterations"].GetInt();
        double NRtolerance = rParameters["newton_raphson"]["global_tolerance"].GetDouble();
        int BSmaxNumOfIterations = rParameters["newton_raphson"]["boundary_iterations"].GetInt();
        double BStolerance = rParameters["newton_raphson"]["boundary_tolerance"].GetDouble();
        int NRBmaxNumOfIterations = rParameters["newton_raphson"]["bisection_iterations"].GetInt();
        double NRBtolerance = rParameters["newton_raphson"]["bisection_tolerance"].GetDouble();

        int gpTria = rParameters["gauss_integration_points_triangle"].GetInt();
        int gpQuad = rParameters["gauss_integration_points_quad"].GetInt();

        bool isWeakDirichletCurve = rParameters["weak_dirichlet_conditions_curve"]["curves"].GetBool();
        bool isWeakDirichletSurface = rParameters["weak_dirichlet_conditions_curve"]["surfaces"].GetBool();
        bool WDCisAutoPenalty = rParameters["weak_dirichlet_conditions_curve"]["auto_penalty"].GetBool();
        double WDCalphaPrim = rParameters["weak_dirichlet_conditions_curve"]["alpha_primary"].GetDouble();
        double WDCalphaSecBending = rParameters["weak_dirichlet_conditions_curve"]["alpha_secondary_bending"].GetDouble();
        double WDCalphaSecTwisting = rParameters["weak_dirichlet_conditions_curve"]["alpha_secondary_twisting"].GetDouble();

        bool isWeakPatchContCond = rParameters["patch_coupling"]["weak_patch_continuity"].GetBool();
        bool isAutoPenalty = rParameters["patch_coupling"]["auto_penalty"].GetBool();
        double alphaPrim = rParameters["patch_coupling"]["alpha_primary"].GetDouble();
        double alphaSecBending = rParameters["patch_coupling"]["alpha_secondary_bending"].GetDouble();
        double alphaSecTwisting = rParameters["patch_coupling"]["alpha_secondary_twisting"].GetDouble();

        bool isErrorComputation = rParameters["error_computation"]["error_computation"].GetBool();
        bool domainError = rParameters["error_computation"]["domain_error"].GetBool();
        bool interfaceError = rParameters["error_computation"]["interface_error"].GetBool();
        bool curveError = rParameters["error_computation"]["curve_error"].GetBool();

        Empire::setParametersConsistency(mapper_name, enforceConsistency, tolConsistency);
        Empire::setParametersProjection(mapper_name, maxProjectionDistance, noInitialGuess, maxProjectionDistanceOnDifferentPatches);
        Empire::setParametersNewtonRaphson(mapper_name, NRmaxNumOfIterations, NRtolerance);
        Empire::setParametersNewtonRaphsonBoundary(mapper_name, NRBmaxNumOfIterations, NRBtolerance);
        Empire::setParametersBisection(mapper_name, BSmaxNumOfIterations, BStolerance);
        Empire::setParametersIntegration(mapper_name, gpTria, gpQuad);
        Empire::setParametersWeakDirichletConditions(mapper_name, isWeakDirichletCurve, isWeakDirichletSurface, WDCisAutoPenalty, WDCalphaPrim, WDCalphaSecBending, WDCalphaSecTwisting);
        Empire::setParametersWeakPatchContinuityConditions(mapper_name, isWeakPatchContCond, isAutoPenalty, alphaPrim, alphaSecBending, alphaSecTwisting);
        Empire::setParametersErrorComputation(mapper_name, isErrorComputation, domainError, interfaceError, curveError);
    }

    void DefinePatchContinuityConditionsOnNurbsShape(std::vector<char>& rMeshName, TopoDS_Shape& rNURBSShape)
    {
        // Create an empty IndexedDataMapOfShapeListOfShape
        TopTools_IndexedDataMapOfShapeListOfShape map_of_edges_to_faces;

        // Make a map containing the ancestor faces of each edge #4Face #6Edge and fill it into mapfOfShapes
        // TopTools_IndexedMapOfShapeListOfShape uses TopoDS_Shape as a key and TopTools_ListOfShape as a value
        TopExp::MapShapesAndAncestors(rNURBSShape, TopAbs_EDGE, TopAbs_FACE, map_of_edges_to_faces);
    
        // Create a faceExplorer to find out the index of each face
        TopExp_Explorer face_explorer (rNURBSShape, TopAbs_FACE);
    
        std::map<int, std::vector<TopoDS_Shape> > coupled_faces_on_edge;
        std::map<int, std::vector<int> > coupled_faces_on_edge_indices;
        int iConnection = 0;
        
        // Loop over all edges
        for (int iEdge=1; iEdge < map_of_edges_to_faces.Extent()+1; iEdge++)
        {

        
            // If the edge has more than one ancestors (faces) then it is a coupling edge
            if ( map_of_edges_to_faces.FindFromIndex(iEdge).Extent() > 1)
            {
                std::vector<TopoDS_Shape> vec_topods_shape;
                std::vector<int> vec_uint;

                // Loop over the ancestor faces of each edge
                for (int iCoupledFace = 1; iCoupledFace < map_of_edges_to_faces.FindFromIndex(iEdge).Extent()+1; iCoupledFace++)
                {
                    // Create a faceExplorer to find out the index of each face
                    face_explorer.ReInit();
                    int faceCtr = 0;
                    // Loop over all faces of the model
                    while (face_explorer.More())
                    {
                        const TopoDS_Shape current_face = face_explorer.Current();
                        face_explorer.Next();
                        // If the considered face matches with the ancestor face then add it to the coupled face list of the edge
                        if (current_face.IsEqual(map_of_edges_to_faces.FindFromIndex(iEdge).First()))
                        {
                            vec_topods_shape.push_back(map_of_edges_to_faces.FindFromIndex(iEdge).First());
                            vec_uint.push_back(faceCtr);
                            
                            // Remove the first ancestor from the list
                            map_of_edges_to_faces.ChangeFromIndex(iEdge).RemoveFirst();
                            break;
                        }
                        
                        // Increment the face counter
                        faceCtr += 1;
                    }
                }

                coupled_faces_on_edge.insert(std::pair<int,std::vector<TopoDS_Shape> >(iEdge, vec_topods_shape));
                coupled_faces_on_edge_indices.insert(std::pair<int, std::vector<int> >(iEdge, vec_uint));
            }
        }
        
        // Loop over the coupling edges and faces to get their definitions on the corresponding patches
        for(auto const &it : coupled_faces_on_edge) 
        {
            int iEdge = it.first; 
            std::vector<TopoDS_Shape> coupled_faces = it.second;

            // Get the coupling edge
            TopoDS_Shape coupling_edge = map_of_edges_to_faces.FindKey(iEdge);
            
            // Define the first patch found to be the master patch containing the master curve
            int master_face_index = coupled_faces_on_edge_indices[iEdge][0];
            TopoDS_Shape& master_face = coupled_faces[0];
            
            // Find the master wire index
            int master_face_wire_ctr = 0;
            bool found_master_wire = false;
            int master_face_wire_edge_ctr = 0;
            TopExp_Explorer master_face_wire_explorer(master_face, TopAbs_WIRE);
            while (master_face_wire_explorer.More())
            {
                const TopoDS_Shape master_wire = master_face_wire_explorer.Current();
                master_face_wire_explorer.Next();
                TopExp_Explorer master_face_wire_edge_explorer(master_wire, TopAbs_EDGE);
                while (master_face_wire_edge_explorer.More())
                {
                    const TopoDS_Shape master_edge = master_face_wire_edge_explorer.Current();
                    master_face_wire_edge_explorer.Next();
                    if (master_edge.IsSame(coupling_edge))
                    {
                        found_master_wire = true;
                        break;
                    }
                    else
                        master_face_wire_edge_ctr += 1;
                }
                
                if (found_master_wire)
                    break;
                else
                    master_face_wire_ctr += 1;
            }
            for (int iSlaveFace = 0; iSlaveFace < (int)coupled_faces.size(); iSlaveFace++)     
            {
                int slave_face_index = coupled_faces_on_edge_indices[iEdge][iSlaveFace];
                TopoDS_Shape slave_face = coupled_faces[iSlaveFace];
                
                // Find the slave wire index
                int slave_face_wire_ctr = 0;
                bool found_slave_wire = false;
                int slave_face_wire_edge_ctr = 0;
                TopExp_Explorer slave_face_wire_explorer(slave_face, TopAbs_WIRE);
                while (slave_face_wire_explorer.More())
                {
                    const TopoDS_Shape slave_wire = slave_face_wire_explorer.Current();
                    slave_face_wire_explorer.Next();
                    
                    TopExp_Explorer slave_face_wire_edge_explorer(slave_wire, TopAbs_EDGE);
                    while (slave_face_wire_edge_explorer.More())
                    {
                        const TopoDS_Shape slave_edge = slave_face_wire_edge_explorer.Current();
                        slave_face_wire_edge_explorer.Next();
                        if (slave_edge.IsSame(coupling_edge))
                        {
                            found_slave_wire = false;
                            break;
                        }
                        else 
                            slave_face_wire_edge_ctr += 1;
                    }
                
                    if (found_slave_wire)
                        break;
                    else
                        slave_face_wire_ctr += 1;
                }
                // print "Defining patch continuity condition: ", iConnection
                // print " Master Patch index: ", masterFaceIndex, " BL index: ", masterFaceWireCtr, " curve index: ", masterFaceWireEdgeCtr
                // print " Slave Patch index: ", slaveFaceIndex, " BL index: ", slaveFaceWireCtr, " curve index: ", slaveFaceWireEdgeCtr
                
                Empire::addPatchContinuityConditionToIGAMesh(
                    rMeshName.data(), 
                    iConnection,
                    master_face_index, 
                    master_face_wire_ctr, 
                    master_face_wire_edge_ctr,
                    slave_face_index, 
                    slave_face_wire_ctr, 
                    slave_face_wire_edge_ctr
                );
                iConnection += 1;
            }
        }
    }

    void GetAllControlPoints(TopoDS_Shape& rNURBSShape, std::vector<double>& rControlPoints )
    {

        // Face explorer to gather the CPs
        TopExp_Explorer face_explorer(rNURBSShape, TopAbs_FACE);
        while (face_explorer.More())
        {
            const TopoDS_Shape currentFace = face_explorer.Current();
            face_explorer.Next();
            std::vector<gp_Pnt> face_control_points;
            GetFaceControlPoints(currentFace, face_control_points);

            for (gp_Pnt& control_point : face_control_points)
            {
                rControlPoints.push_back(control_point.X());
                rControlPoints.push_back(control_point.Y());
                rControlPoints.push_back(control_point.Z());
            }
        }
    }

    void GetFaceControlPoints(const TopoDS_Shape& face, std::vector<gp_Pnt>& rFaceControlPoints)
    {
        // Convert a topological TopoDS_FACE into Geom_BSplineSurface
        Handle(Geom_BSplineSurface) bspline_surface = GeomConvert::SurfaceToBSplineSurface(BRep_Tool::Surface(TopoDS::Face(face)));

        // Get the number of CPs in each direction
        int UnoCP = bspline_surface->NbUPoles();
        int VnoCP = bspline_surface->NbVPoles();
        
        // Get the CPs of the Geom_BSplineSurface
        TColgp_Array2OfPnt patch_cp_net(1, UnoCP,1, VnoCP);
        bspline_surface->Poles(patch_cp_net);
        
        // Add all the CPs into a list
        rFaceControlPoints.clear();
        for (int j=1; j < VnoCP+1 ; j++)
        {
            for (int i=1; i < UnoCP+1; i++)
            {
                rFaceControlPoints.push_back(patch_cp_net.Value(i,j));
            }
        }   
    }

    // --------------------------------------------------------------------------
    void AssignMappingIds()
    {
        unsigned int i = 0;
        for(auto& node_i : mrDesignSurface.Nodes())
            node_i.SetValue(MAPPING_ID,i++);
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION defined */
