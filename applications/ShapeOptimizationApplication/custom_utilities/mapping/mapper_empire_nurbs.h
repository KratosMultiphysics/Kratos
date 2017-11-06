//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
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


#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <GeomConvert.hxx>
#include <Geom2dConvert.hxx>
#include <Geom_BSplineSurface.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TopoDS_Wire.hxx>
#include <ShapeAnalysis.hxx>

// boost includes
#include <boost/algorithm/string.hpp>

// EMPIRE includes
namespace Empire
{
    #include <MapperLib.hpp>
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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MapperEmpireNURBS( 
        ModelPart& rDesignSurface,
        Parameters& rParameters )
    {
        KRATOS_TRY();

        Parameters default_params(R"(
            {
                "design_suface_extrusion": 1.0,
                "iga_model_import":
                {
                    "input_type": "PLEASE_SPECIFY_INPUT_MODEL_FILE_TYPE",
                    "input_filename": "PLEASE_SPECIFY_INPUT_FILENAME",
                    "echo_level": 0
                }
            })");
        
        rParameters.ValidateAndAssignDefaults(default_params);

        mDesignSufaceExtrusionLength = rParameters["design_suface_extrusion"].GetDouble();
        mrDesignSurface = rDesignSurface;



        //initiating the sequence to create the mapping between IGA mesh and FE mesh
        rExternalModel->Execute();
        this->MakeIGAMeshFromNURBSSHape();

        KRATOS_CATCH("");

    }

    /// Destructor.
    virtual ~MapperEmpireNURBS()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ReadDesignSurfaceFEMesh(std::map<unsigned int, std::vector<double> >& rNodes, std::map<unsigned int, std::vector<unsigned int> >& rElements)
    {
        rNodes.clean();
        rElements.clean();
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
            for (ModelPart::NodeIterator node_i = model_part.NodesBegin(); node_i != model_part.NodesEnd(); ++node_i)
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

        // libEmpireMapper_api = cdll.LoadLibrary(os.environ['EMPIRE_MAPPER_LIBSO_ON_MACHINE'])
        // initFEMesh = Empire::initFEMesh
        // setNodesToFEMesh = Empire::setNodesToFEMesh
        // setElementsToFEMesh = Empire::setElementsToFEMesh

        const std::string& mesh_name = mFEMeshName;

        ReadDesignSurfaceFEMesh(node_list, element_list);

        const unsigned int number_of_nodes = node_list.size();
        const unsigned int number_of_elements = element_list.size();

        // numNodes_c = len(_nodes))
        // numElems_c = len(_elems))

        Empire::initFEMesh(mesh_name.c_str(), number_of_nodes, number_of_elements, false);
        // initFEMesh(meshName_c,numNodes_c,numElems_c,False)


        std::vector<int> node_id_list;
        std::vector<double> node_values_list;        
        node_id_list.clear();
        node_values_list.clear();

        for(std::map<unsigned int, std::vector<double> >::const_iterator it =node_list.begin(); it != node_list.end(); it++)
        {
            int key = (int)it->first;
            std::vector<double>& value = it->second;
            
            node_id_list.push_back(key);

            for (unsigned int i=0; i<3; i++)
                node_values_list.push_back(value[i]);
        }
        

        // nodeIDs = _nodes.keys()
        // _nodes = sum(_nodes.values(),[])
        // nodeIDs_c = (len(nodeIDs))(0)
        // nodes_c = (len(_nodes))(0.0)
        // for nodeIDctr in range(0,len(nodeIDs)):
        //   nodeIDs_c[nodeIDctr] = nodeIDs[nodeIDctr]
        // for nodeCtr in range(0,len(_nodes)):
        //   nodes_c[nodeCtr] = _nodes[nodeCtr]

        Empire::setNodesToFEMesh(mesh_name.c_str(), node_id_list.data(), node_values_list.data());
       
        std::vector<int> number_of_nodes_per_elem;
        std::vector<int> element_node_ids;

        for(std::map<unsigned int, std::vector<double> >::const_iterator it = element_list.begin(); it != element_list.end(); it++)
        {
            std::vector<int>& value = it->second;
            number_of_nodes_per_elem.push_back( value.size() );

            for (unsigned int i=0; i < value.size(); i++)
                element_node_ids.push_back(value[i]);
        }        

        // numNodesPerElem = []
        // elemsToSend = []
        // for elemID in _elems:
        //   numNodesPerElem.append(len(_elems[elemID]['nodes']))
        //   for nodeID in _elems[elemID]['nodes']:
        //     elemsToSend.append(int(nodeID))
            
        // numNodesPerElem_c = (len(numNodesPerElem))(0)
        // for someCtr in range(0,len(numNodesPerElem)):
        //   numNodesPerElem_c[someCtr] = int(numNodesPerElem[someCtr])
        
        // elemsToSend_c = (len(elemsToSend))(0)
        // for someCtr in range(0,len(elemsToSend)):
        //   elemsToSend_c[someCtr]=elemsToSend[someCtr]
        // print "Setting elements to FEMesh"
        Empire::setElementsToFEMesh(mesh_name.c_str(), number_of_nodes_per_elem.data(), element_node_ids.data());
    }

    
    void MakeIGAMeshFromNURBSSHape()
    {
        KRATOS_TRY;

        const std::string& mesh_name = mIGAMeshName;

        Empire::initIGAMesh(mesh_name.c_str());
        TopExp_Explorer face_explorer(mrExternalModel.GetCompoundNURBSList(), TopAbs_FACE);
        IndexType patch_index = 0;
        IndexType unique_id_counter = 0;

        std::vector<std::vector<double>> patch_unique_id_nets;

        while (face_explorer.More())
        {
            patch_index++;
            TopoDS_Shape& current_shape = face_explorer.Current();
            face_explorer.Next();

            if ( mEchoLevel > 0)
                std::cout<<" Gathering information of patch index: "<<patch_index-1<<" Surface: "<<current_shape<<std::endl;
            
            int p_degree, q_degree, number_u_cps, number_v_cps, number_face_wires;
            std::vector<double> u_knots_vector;
            std::vector<double> v_knots_vector;
            std::vector<double> cp_net;
            std::vector<IndexType> unique_id_net;

            GetFaceData( 
                            current_shape,
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
            
            patch_unique_id_nets.append(unique_id_net);
            unique_id_counter += unique_id_net.size();

            if (mEchoLevel>0)
                std::cout<<" adding patch to EMPIRE"<<std::endl;

            Empire::addPatchToIGAMesh(
                mesh_name.c_str(),
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

            TopoDS_Wire outer_wire = ShapeAnalysis::OuterWire(current_shape);

            TopExp_Explorer face_wire_explorer(current_shape, TopAbs_WIRE);
            while (face_wire_explorer.More())
            {
                TopoDS_Shape& current_face_wire = face_wire_explorer.Current();
                face_wire_explorer.Next();

                if (mEchoLevel>0)
                    std::cout<<" Gathering information of trimming: "<<current_face_wire<<std::endl;
                
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
                    mesh_name.c_str(),
                    patch_index-1,
                    is_inner_curve,
                    number_trimming_curves

                );

                TopExp_Explorer face_wire_edge_explorer(current_face_wire, TopAbs_EDGE);
                while (face_wire_edge_explorer.More())
                {
                    TopoDS_Shape& current_face_wire_edge = face_wire_edge_explorer.Current();
                    face_wire_edge_explorer.Next();

                    if (mEchoLevel>0)
                        std::cout<<" Gathering information of trimming curve edge: "<<current_face_wire_edge<<std::endl;

                    int direction, p_degree, number_u_cps;
                    double active_range_begin, active_range_end;
                    std::vector<double> u_knots_vector;
                    std::vector<double> cp_net;
                    
                    GetEdgeDataWithOrientation (
                        current_face_wire_edge,
                        current_face_wire,
                        current_shape,
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
                        mesh_name.c_str(),
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

            Empire::linearizeTrimmingLoops(mesh_name.c_str(), patch_index-1);
        }

        KRATOS_CATCH("");
    }

    void GetFaceData(
        const TopoDS_Shape& Face,
        IndexType StartID,
        IndexType& rNoOfFaceWires,
        IndexType& rPDegree,
        IndexType& rQDegree,
        IndexType& rUNoCPs,
        IndexType& rVNoCPs,
        std::vector<double>& rUKnotVector,
        std::vector<double>& rVKnotVector,
        std::vector<double>& rCPNet,
        std::vector<IndexType>& rUniqueIDNet
    )
    {
        KRATOS_TRY;

        rNoOfFaceWires = 0;
        TopExpExplorer face_wire_explorer(Face, TopAbs_WIRE);
        while (face_wire_explorer.More())
        {
            face_wire_explorer.Next();
            rNoOfFaceWires++;
        }

        // Convert a topological TopoDS_FACE into Geom_BSplineSurface
        Handle<Geom_Surface>& surface_handle = BRep_Tool.Surface(topods_Face(Face));
        Geom_BSplineSurface& bspline_surface = SurfaceToBSplineSurface(surface_handle).GetObject();
          
        // Gather patch info
        
        // Polynomial degrees
        int rPDegree = bspline_surface.UDegree();
        int rQDegree = bspline_surface.VDegree();

        // No of Knots
        int number_u_knots = bspline_surface.NbUKnots();
        int number_v_knots = bspline_surface.NbVKnots();

        rUKnotVector.clear();
        for (IndexType i=1; i<number_u_knots+1; i++)
        {
            for (IndexType j=0; j < bspline_surface.UMultiplicity(i); j++)
            {
                rUKnotVector.push(bspline_surface.UKnot(i));
            }
        }

        rVKnotVector.clear();
        for (IndexType i=1; i<number_v_knots+1; i++)
        {
            for (IndexType j=0; j < bspline_surface.VMultiplicity(i); j++)
            {
                rUKnotVector.push(bspline_surface.VKnot(i));
            }
        }

        if bspline_surface.IsUPeriodic():
        {
            if (mEchoLevel>0)
                std::cout<<" found periodic U."<<std::endl;
            rUKnotVector.insert(0, rUKnotVector[0]);
            rUKnotVector.push(rUKnotVector[rUKnotVector.size()-1]);
            bspline_surface.SetUNotPeriodic();
        }
        if bspline_surface.IsVPeriodic():
        {
            if (mEchoLevel>0)
                std::cout<<" found periodic V."<<std::endl;
            rVKnotVector.insert(0, rVKnotVector[0]);
            rVKnotVector.push(rVKnotVector[rVKnotVector.size()-1]);
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
        OCC_bsplineSurface.Weights(occ_cp_weight_net);

        rCPNet.clear();
        rUniqueIDNet.clear();

        for (IndexType iVCP=1; iVCP < rVNoCPs+1; iVCP++)
            for (IndexType iUCP=1; iUCP < rUNoCPs+1; iUCP++)
            {
                gp_Pnt& occ_cp = occ_cp_net.Value(iUCP, iVCP);
                double occ_cp_weight = occ_cp_weight_net.Value(iUCP,iVCP);
                rCPNet.push(occ_cp.X());
                rCPNet.push(occ_cp.Y());
                rCPNet.push(occ_cp.Z());
                rCPNet.push(occ_cp_weight);
                rUniqueIDNet.push(StartID++);
            }
    }

    void GetWireData(
        const TopoDS_Shape& rOuterWire, 
        const TopoDS_Shape& rCurrentWire,
        bool& rIsInnerWire,
        IndexType& rNoTrimmingCurves
    )
    {
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
    }

    void GetEdgeDataWithOrientation(
        const TopoDS_Shape& rEdge,
        const TopoDS_Shape& rWire,
        const TopoDS_Shape& rFace,
        IndexType& rDirection,
        IndexType& rPDegree,
        std::vector<double>& rUKnotVector,
        IndexType& rUNoCPs,
        std::vector<double>& rCPNet,
        double& rActiveRangeBegin,
        double& rAcriveRangeEnd
    )
    {
        // Convert a topological TopoDS_EDGE into Geom2d_BSplineCurve
        Handle<Geom2d_Curve> curve_2d_handle = BRep_Tool::CurveOnSurface(rEdge, rFace);
        Geom2d_BSplineCurve& bspline_curve_2d = Geom2dConvert::CurveToBsplineCurve(curve_2d_handle).GetObject();
        
        if bspline_curve_2d.IsPeriodic():
        {
            std::cout<<" Found periodic curve. setting it to non-periodic curve."<<std::endl;
            bspline_curve_2d.SetNotPeriodic();
        }
        
        // Curve direction 
        rDirection = (rWire.Orientation() + rEdge.Orientation() + 1)%2;
        // pDegree
        rPDegree = bspline_curve_2d.Degree();
        // No of Knots
        IndexType number_of_knots = bspline_curve_2d.NbKnots();
        // Knot vector
        TColStd_Array1OfInteger occ_u_knot_multi(1, number_of_knots);
        bspline_curve_2d.Multiplicities(occ_u_knot_multi);

        IndexType number_u_knots = 0;
        for (IndexType ctr=1; ctr<number_of_knots+1; ctr++)
            number_u_knots += occ_u_knot_multi.Value(ctr);

        TColStd_Array1OfReal occ_u_knot_sequence(1, number_u_knots);
        bspline_curve_2d.KnotSequence(occ_u_knot_sequence);

        rUKnotVector.clear();
        for (IndexType i=1; i<number_u_knots+1; i++)
            rUKnotVector.push(occ_u_knot_sequence.Value(i));

        // No of CPs
        rUNoCPs = bspline_curve_2d.NbPoles();
        // CPs
        TColgp_Array1OfPnt2d occ_cp_net(1, rUNoCPs);
        bspline_curve_2d.Poles(occ_cp_net);
        // CP weights
        TColStd_Array1OfReal occ_cp_weight_net(1, rUNoCPs);
        bspline_curve_2d.Weights(occ_cp_weight_net)
        
        rCPNet.clear();
        for (IndexType i=1; i<rUNoCPs+1; i++)
        {
            gp_Pnt2d& occ_cp = occ_cp_net.Value(i);
            occ_cp_weight = occ_cp_weight_net.Value(i);
            rCPNet.push(occ_cp.X());
            rCPNet.push(occ_cp.Y());
            rCPNet.push(0.0);
            rCPNet.push(occ_cp_weight);
        }
    }

    void GenerateMapper(const std::string& rrMapperName.c_str())
    {
        std::string mapper_name;
        std::string mesh_name_a;
        std::string mesh_name_b;
     
        if (rrMapperName.c_str().compare(mMapperIGAToFe)==0) 
        {
            mapper_name = rrMapperName.c_str();
            mesh_name_a = mIGAMeshName;
            mesh_name_b = mFEMeshName;
        }
        else if (rrMapperName.c_str().compare(mMapperFEToIGA)==0)
        {
            mapper_name = rrMapperName.c_str();
            mesh_name_a = mFEMeshName;
            mesh_name_b = mIGAMeshName;
        }


        Empire::initIGAMortarMapper(mappe_name.c_str(), mesh_name_a.c_str(), mesh_name_b.c_str());
        Empire::SetMapperParameters(mapper_name);

        Empire::initialize(mapper_name.c_str());
        Empire::buildCouplingMatrices(mapper_name.c_str());  
        
        
        IGAName, FEMName, rMapperName.c_str(), \
        enforceConsistency, tolConsistency, \
        isIGA2FEM, \
        maxProjectionDistance, noInitialGuess, maxProjectionDistanceOnDifferentPatches, \
        NRmaxNumOfIterations, NRtolerance, \
        NRBmaxNumOfIterations, NRBtolerance, \
        BSmaxNumOfIterations, BStolerance, \
        gpTria, gpQuad, \
        isWeakDirichletCurve, isWeakDirichletSurface, WDCisAutoPenalty, WDCalphaPrim, WDCalphaSecBending, WDCalphaSecTwisting, \
        isWeakPatchContCond, isAutoPenalty, alphaPrim, alphaSecBending, alphaSecTwisting, \
        isErrorComputation, domainError, interfaceError, curveError = getIGAMortarMapperParameters()
          
    }

    void SetMapperParameters(const std::string& rMapperName)
    {
        bool enforceConsistency;
        double tolConsistency;
        double maxProjectionDistance;
        int noInitialGuess;
        double maxProjectionDistanceOnDifferentPatches;
        int NRmaxNumOfIterations;
        double NRtolerance;
        int BSmaxNumOfIterations;
        double BStolerance;
        int gpTria;
        int gpQuad;
        bool isWeakDirichletCurve;
        bool isWeakDirichletSurface;
        bool WDCisAutoPenalty;
        double WDCalphaPrim;
        double WDCalphaSecBending;
        double WDCalphaSecTwisting;
        bool isWeakPatchContCond;
        double alphaPrim;
        double alphaSecBending;
        double alphaSecTwisting;
        bool isErrorComputation;
        bool domainError;
        bool interfaceError;
        bool curveError;



        Empire::setParametersConsistency(rMapperName.c_str(), enforceConsistency, tolConsistency);
        Empire::setParametersProjection(rMapperName.c_str(), maxProjectionDistance, noInitialGuess, maxProjectionDistanceOnDifferentPatches);
        Empire::setParametersNewtonRaphson(rMapperName.c_str(), NRmaxNumOfIterations, NRtolerance);
        Empire::setParametersNewtonRaphsonBoundary(rMapperName.c_str(), NRBmaxNumOfIterations, NRBtolerance);
        Empire::setParametersBisection(rMapperName.c_str(), BSmaxNumOfIterations, BStolerance);
        Empire::setParametersIntegration(rMapperName.c_str(), gpTria, gpQuad);
        Empire::setParametersWeakDirichletConditions(rMapperName.c_str(), isWeakDirichletCurve, isWeakDirichletSurface, WDCisAutoPenalty, WDCalphaPrim, WDCalphaSecBending, WDCalphaSecTwisting);
        Empire::setParametersWeakPatchContinuityConditions(rMapperName.c_str(), isWeakPatchContCond, isAutoPenalty, alphaPrim, alphaSecBending, alphaSecTwisting);
        Empire::setParametersErrorComputation(rMapperName.c_str(), isErrorComputation, domainError, interfaceError, curveError);
    }

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
    ExternalModelImport& mrExternalModel;
    ModelPart& mrDesignSurface;
    double mDesignSufaceExtrusionLength;
    int mEchoLevel;


    const std::string mFEMeshName = "FE_MESH";
    const std::string mIGAMeshName = "IGA_MESH";
    const std::string mMapperFEToIGA = "FE_TO_IGA";
    const std::string mMapperIGAToFe = "IGA_TO_FE";
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION defined */
