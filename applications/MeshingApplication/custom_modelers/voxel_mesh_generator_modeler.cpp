#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/model_part.h"

#include "geometries/bounding_box.h"
#include "geometries/point.h"

#include "utilities/parallel_utilities.h"
#include "voxel_mesh_generator_modeler.h"
//#include "additive_manufacturing_application_variables.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "includes/dem_variables.h"


namespace Kratos
{
    typedef std::size_t SizeType;
/**
 * Constructor.
 */
    /// Default constructor.
    VoxelMeshGeneratorModeler::VoxelMeshGeneratorModeler() : Modeler()
    {
    }

    /// Constructor.
    VoxelMeshGeneratorModeler::VoxelMeshGeneratorModeler(
        Model& rModel, Parameters rParameters)
        : Modeler(rModel, rParameters )
        , mpModel(&rModel)
    {
        rParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    }

    ///
    const Parameters VoxelMeshGeneratorModeler::GetDefaultParameters() const
    {
        return Parameters(R"(
        {
            "key_plane_generator": {},
            "entities_generator_list": [],
            "coloring_settings_list": [],
            "output_filename" : "",
            "mdpa_file_name" : "",
            "output_model_part_name" : "",
            "input_model_part_name" : "",
            "voxel size" : 0.001
        }  )");
    }

    ///
    void VoxelMeshGeneratorModeler::SetupGeometryModel() {

        KRATOS_ERROR_IF_NOT( mParameters.Has("input_model_part_name") )
            << "Missing \"input_model_part_name\" in VoxelMeshGeneratorModeler Parameters." << std::endl;

        const std::string input_model_part_name = mParameters["input_model_part_name"].GetString();
        mpInputModelPart = mpModel->HasModelPart(input_model_part_name)
            ? &mpModel->GetModelPart(input_model_part_name)
            : &mpModel->CreateModelPart(input_model_part_name);

        KRATOS_ERROR_IF_NOT(mParameters.Has("mdpa_file_name")) << "mdpa_file_name not defined" << std::endl;

        

        const std::string DataFileName =  mParameters["mdpa_file_name"].GetString();

        KRATOS_INFO_IF("::[VoxelMeshGeneratorModeler]::", mEchoLevel > 0) << "Importing Cad Model from: " << DataFileName << std::endl;

        // Load the mdpa
        if(DataFileName!=""){
            ModelPartIO(DataFileName).ReadModelPart(*mpInputModelPart);
        }

        // Create the target model
        const std::string output_model_part_name = mParameters["output_model_part_name"].GetString();
        ModelPart& main_model_part = mpModel->HasModelPart(output_model_part_name)
             ? mpModel->GetModelPart(output_model_part_name)
             : mpModel->CreateModelPart(output_model_part_name);

        KRATOS_INFO("Modeler") << "Generating Key planes..." << std::endl;
        VoxelMeshKeyPlaneGeneratorBySize(mParameters["key_plane_generator"]["Parameters"]);
        KRATOS_INFO("Modeler") << "Key Planes generated" << std::endl;

        KRATOS_INFO("Modeler") << "Preparing Internal Data Structure" << std::endl;
        PraparingTheInternalDataStructure();
        KRATOS_INFO("Modeler") << "Internal Data Structure prepared" << std::endl;

        if(mParameters.Has("coloring_settings_list")){
            ApplyColoring(mParameters["coloring_settings_list"]);
        }
        // mColors.WriteParaViewVTR("voxel_colors.vtr");
        if(mParameters.Has("entities_generator_list")){
            GenerateEntities(main_model_part, mParameters["entities_generator_list"]);
        } 
        // I would still keep this commented line for future debugging.
        mColors.WriteParaViewVTR("Coloring.vtr");
    }


    void VoxelMeshGeneratorModeler::VoxelMeshKeyPlaneGeneratorBySize(Parameters KeyPlaneGeneratorParameters){
        KRATOS_ERROR_IF_NOT(KeyPlaneGeneratorParameters.Has("voxel_sizes")) << "voxel_sizes should be defined for this generator as an array of 3 sizes for x,y,z directions" << std::endl;
        KRATOS_ERROR_IF_NOT(KeyPlaneGeneratorParameters.Has("min_point")) << "min_point should be defined for this generator as an array of 3 coordinates of the min point" << std::endl;
        KRATOS_ERROR_IF_NOT(KeyPlaneGeneratorParameters.Has("max_point")) << "max_point should be defined for this generator as an array of 3 coordinates of the max point" << std::endl;

        array_1d<double,3> voxel_sizes = KeyPlaneGeneratorParameters["voxel_sizes"].GetVector();
        array_1d<double,3> min_point = KeyPlaneGeneratorParameters["min_point"].GetVector();
        array_1d<double,3> max_point = KeyPlaneGeneratorParameters["max_point"].GetVector();

        for(int i_direction = 0 ; i_direction < 3 ; i_direction++){
            double input_voxel_size = voxel_sizes[i_direction];
            const double min_coordinate = min_point[i_direction];
            const double max_coordinate = max_point[i_direction];
            
            KRATOS_ERROR_IF(input_voxel_size == 0.00) << "voxel_sizes in direction " << i_direction << " cannot be 0.00";

            const double length = (max_coordinate - min_coordinate);

            KRATOS_ERROR_IF_NOT(length>0.0) << "Negative or zero length of voxelization bounding box in " << i_direction << " direction" << std::endl;

            std::size_t number_of_divisions = static_cast<std::size_t>(std::round(length / input_voxel_size));
            number_of_divisions= (number_of_divisions == 0) ? 1 : number_of_divisions;

            const double voxel_size = length / number_of_divisions;
            
            for( std::size_t i = 0 ; i < number_of_divisions + 1 ; i++){
                mKeyPlanes[i_direction].push_back(min_coordinate + i*voxel_size);
            }
        }

    }


    void VoxelMeshGeneratorModeler::PraparingTheInternalDataStructure(){
        mColors.SetCoordinates(mKeyPlanes[0], mKeyPlanes[1], mKeyPlanes[2]);
        mMeshingData.SetNumberOfDivisions(mKeyPlanes[0].size(), mKeyPlanes[1].size(), mKeyPlanes[2].size());

        array_1d<double,3> margins = ZeroVector(3);
        for(int i = 0 ; i < 3 ; i++) {
            margins[i] = (mKeyPlanes[i].back() - mKeyPlanes[i].front()) * 1.0e-2;
        }

        const double margin = *std::max_element(margins.begin(), margins.end());
        
        mColors.ExtendBoundingBox(mpInputModelPart->Nodes(), margin);
    }

    void VoxelMeshGeneratorModeler::ApplyColoring(Parameters ColoringParameters){
        Timer::Start("Voxel Mesh Coloring");

        double outside_color = 1.0;
        if(ColoringParameters.Has("outside_color")){
            outside_color = ColoringParameters["outside_color"].GetDouble();
        }
        mColors.SetAllColors(outside_color);

        for(auto parameters : ColoringParameters){

            Parameters default_parameters(R"(
                {
                    "type" : "cells_with_inside_center",
                    "model_part_name": "Undefined",
                    "color": -1,
                    "cell_color": -1,
                    "input_entities": "elements"
                }  )");

            parameters.ValidateAndAssignDefaults(default_parameters);
            
            std::string model_part_name = parameters["model_part_name"].GetString();

            std::string type = parameters["type"].GetString();

            KRATOS_INFO("Modeler") << "Applying color to " << type << " for " << model_part_name  << " model part" << std::endl;

            if(type == "cells_with_inside_center") {
                ModelPart& skin_part = mpModel->GetModelPart(model_part_name);
                ApplyColorToCellsWithInsideCenter(skin_part, parameters, outside_color);
            }
            else if(type == "cells_in_touch") {
                ModelPart& skin_part = mpModel->GetModelPart(model_part_name);
                ApplyColorToCellsInTouch(skin_part, parameters, outside_color);
            }
            else if(type == "cells_faces") {
                ModelPart& skin_part = mpModel->GetModelPart(model_part_name);
                ApplyColorToCellsFaces(skin_part, parameters, outside_color);
            }
            else if(type == "outer_faces_of_cells_with_color") {
                ApplyColorToOuterFacesOfCellsWithColor(parameters, outside_color);
            }
            else if(type == "cells_with_center_in_sphere_arround_nodes") {
                ApplyColorToCellsWithCenterInSphereArroundNodes(parameters, outside_color);
            }
            else if(type == "connected_cells_in_touch") {
                ModelPart& skin_part = mpModel->GetModelPart(model_part_name);
                ApplyColorToConnectedCellsInTouch(skin_part, parameters);
            }
            else{
                KRATOS_ERROR << "The coloring type " << parameters["type"] << " is not supported. The supported types are " 
                << "\"cells_with_inside_center\"" << ", \"cells_in_touch\"" 
                << ", \"cells_faces\"" << ", \"outer_faces_of_cells_with_color\"" 
                << ", \"cells_with_center_in_sphere_arround_nodes\"" << ", and \"connected_cells_in_touch\"" << std::endl;
            }
        }
        Timer::Stop("Voxel Mesh Coloring");
    }

    void VoxelMeshGeneratorModeler::ApplyColorToCellsWithInsideCenter(ModelPart const& TheSkinModelPart, Parameters parameters, double OutsideColor) {

        double inside_color = parameters["color"].GetDouble();
        std::string input_entities = parameters["input_entities"].GetString();
        
        array_1d< std::size_t, 3 > min_ray_position = zero_vector<std::size_t>(3);
        array_1d< std::size_t, 3 > max_ray_position = zero_vector<std::size_t>(3);
            
        mColors.CalculateMinMaxCenterOfElementPositions(TheSkinModelPart.Nodes(), min_ray_position, max_ray_position);

        mColors.InitializeRays(min_ray_position, max_ray_position, "center_of_elements");

        if(input_entities == "elements"){
            for(auto& element : TheSkinModelPart.Elements())
            {
                Element::GeometryType& r_geometry = element.GetGeometry();
                CheckGeometryType(r_geometry.GetGeometryType());
                mColors.AddGeometry(r_geometry, false);
            }
        }
        else if(input_entities == "conditions"){
            for(auto& condition : TheSkinModelPart.Conditions())
            {
                Condition::GeometryType& r_geometry = condition.GetGeometry();
                CheckGeometryType(r_geometry.GetGeometryType());
                mColors.AddGeometry(r_geometry, false);
            }
        }
        else{
                KRATOS_ERROR << "The input_entities  " << parameters["input_entities"] << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
        }

        mColors.CalculateElementalRayColors(min_ray_position, max_ray_position, inside_color, OutsideColor);
    }

    void VoxelMeshGeneratorModeler::ApplyColorToCellsInTouch(ModelPart const& TheSkinModelPart, Parameters parameters, double OutsideColor) {

        double inside_color = parameters["color"].GetDouble();
        std::string input_entities = parameters["input_entities"].GetString();
            
        if(input_entities == "elements"){
            for(auto& i_geometrical_object : TheSkinModelPart.Elements()) {
                auto& r_geometry = i_geometrical_object.GetGeometry(); 
                CheckGeometryType(r_geometry.GetGeometryType());
                ApplyColorToCellsInTouchWithGeometry(r_geometry, inside_color, OutsideColor);
            }
        }
        else if(input_entities == "conditions"){
            for(auto& i_geometrical_object : TheSkinModelPart.Conditions()) {
                auto& r_geometry = i_geometrical_object.GetGeometry(); 
                CheckGeometryType(r_geometry.GetGeometryType());
                ApplyColorToCellsInTouchWithGeometry(r_geometry, inside_color, OutsideColor);
            }
        }
        else{
                KRATOS_ERROR << "The input_entities  " << parameters["input_entities"] << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
        }
    }

    void VoxelMeshGeneratorModeler::ApplyColorToCellsInTouchWithGeometry(Element::GeometryType& TheGeometry, double InsideColor, double OutsideColor) {
        array_1d<std::size_t, 3> min_position(3,0);
        array_1d<std::size_t, 3> max_position(3,0);
        mColors.CalculateOuterMinMaxNodePositions(TheGeometry, min_position, max_position);
        for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
            for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    Point cell_min_point = mColors.GetPoint(i,j,k);
                    Point cell_max_point = mColors.GetPoint(i+1,j+1,k+1);
                    if(TheGeometry.HasIntersection(cell_min_point,cell_max_point)){
                        mColors.GetElementalColor(i,j,k) = InsideColor;
                    }
                }
            }
        }
    }

    void VoxelMeshGeneratorModeler::ApplyColorToCellsFaces(ModelPart const& TheSkinModelPart, Parameters parameters, double OutsideColor) {

        double interface_color = parameters["color"].GetDouble();
        double cell_color = parameters["cell_color"].GetDouble();
        std::string input_entities = parameters["input_entities"].GetString();

        array_1d< std::size_t, 3 > min_ray_position = zero_vector<std::size_t>(3);
        array_1d< std::size_t, 3 > max_ray_position = zero_vector<std::size_t>(3);
            
        mColors.CalculateMinMaxCenterOfElementPositions(TheSkinModelPart.Nodes(), min_ray_position, max_ray_position);

        mColors.InitializeRays(min_ray_position, max_ray_position, "center_of_elements");

        if(input_entities == "elements"){
            for(auto& element : TheSkinModelPart.Elements())
            {
                Element::GeometryType& r_geometry = element.GetGeometry();
                CheckGeometryType(r_geometry.GetGeometryType());
                mColors.AddGeometry(r_geometry, false);
            }
        }
        else if(input_entities == "conditions"){
            for(auto& condition : TheSkinModelPart.Conditions()) {
                auto& r_geometry = condition.GetGeometry(); 
                CheckGeometryType(r_geometry.GetGeometryType());
                mColors.AddGeometry(r_geometry, false);
            }
        }

        mColors.CalculateElementalFaceColors(min_ray_position, max_ray_position, interface_color, OutsideColor, cell_color);
    }

    void VoxelMeshGeneratorModeler::ApplyColorToOuterFacesOfCellsWithColor(Parameters parameters, double OutsideColor) {

        double interface_color = parameters["color"].GetDouble();
        double cell_color = parameters["cell_color"].GetDouble();

        array_1d<std::size_t, 3> number_of_cells;
        for(int i = 0 ; i < 3 ; i++){
            number_of_cells[i]=mKeyPlanes[i].size() - 1;
        }
    
        array_1d<std::size_t, 3> cell_indices;
        for (cell_indices[2] = 0; cell_indices[2] < number_of_cells[2]; cell_indices[2]++) {
            for (cell_indices[1] = 0; cell_indices[1] < number_of_cells[1]; cell_indices[1]++) {
                for (cell_indices[0] = 0; cell_indices[0] < number_of_cells[0]; cell_indices[0]++) {
                    if(mColors.GetElementalColor(cell_indices[0],cell_indices[1],cell_indices[2]) == cell_color){
                        ApplyColorIfOuterFace(interface_color, cell_color, cell_indices);
                    }
                }
            }
        }
    }

    void VoxelMeshGeneratorModeler::ApplyColorToCellsWithCenterInSphereArroundNodes(Parameters parameters, double OutsideColor) {

        double color = parameters["color"].GetDouble();

        for(auto& i_node : mpInputModelPart->Nodes()){
            array_1d< std::size_t, 3 > min_position;
            array_1d< std::size_t, 3 > max_position;
            double radius = i_node.GetSolutionStepValue(RADIUS);
            double radius2 = radius * radius; // storing radius**2 to avoid sqrt in comparison

            for(int i = 0; i < 3; i++ ) {
                min_position[i] = mColors.CalculateCenterOfElementPosition(i_node[i] - radius, i);
                max_position[i] = mColors.CalculateCenterOfElementPosition(i_node[i] + radius, i) + 1;
            }

            for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
                for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                    for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                        Point cell_center = mColors.GetCenterOfElement(i,j,k);
                        array_1d<double, 3> distance_vector = i_node.Coordinates() -  cell_center.Coordinates();
                        if(inner_prod(distance_vector, distance_vector) < radius2){
                            mColors.GetElementalColor(i,j,k) = color;
                        }                       
                    }
                }
            }
        }
    }

    void VoxelMeshGeneratorModeler::ApplyColorToConnectedCellsInTouch(ModelPart const& TheSkinModelPart, Parameters parameters) {

        double inside_color = parameters["color"].GetDouble();
        double cell_color = parameters["cell_color"].GetDouble();
        std::string input_entities = parameters["input_entities"].GetString();
            
        if(input_entities == "elements"){
            for(auto& i_geometrical_object : TheSkinModelPart.Elements()) {
                auto& r_geometry = i_geometrical_object.GetGeometry(); 
                CheckGeometryType(r_geometry.GetGeometryType());
                ApplyColorToConnectedCellsInTouchWithGeometry(r_geometry, inside_color, cell_color);
            }
        }
        else if(input_entities == "conditions"){
            for(auto& i_geometrical_object : TheSkinModelPart.Conditions()) {
                auto& r_geometry = i_geometrical_object.GetGeometry(); 
                CheckGeometryType(r_geometry.GetGeometryType());
                ApplyColorToConnectedCellsInTouchWithGeometry(r_geometry, inside_color, cell_color);
            }
        }
        else{
                KRATOS_ERROR << "The input_entities  " << parameters["input_entities"] << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
        }
    }

    void VoxelMeshGeneratorModeler::ApplyColorToConnectedCellsInTouchWithGeometry(Element::GeometryType& TheGeometry, double InsideColor, double CellColor) {
        array_1d<std::size_t, 3> min_position(3,0);
        array_1d<std::size_t, 3> max_position(3,0);
        mColors.CalculateOuterMinMaxNodePositions(TheGeometry, min_position, max_position);
        for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
            for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
                for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                    auto cell_color = mColors.GetElementalColor(i,j,k);
                    if(cell_color == CellColor){
                        Point cell_min_point = mColors.GetPoint(i,j,k);
                        Point cell_max_point = mColors.GetPoint(i+1,j+1,k+1);
                        if(TheGeometry.HasIntersection(cell_min_point,cell_max_point)){
                            ColorConnectedCellsToThisCell(i,j,k,InsideColor, CellColor);
                        }
                    }
                }
            }
        }
    }

    void VoxelMeshGeneratorModeler::ColorConnectedCellsToThisCell(std::size_t I, std::size_t J, std::size_t K, double InsideColor, double CellColor) {
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> stack;
        stack.emplace_back(I,J,K);
        std::size_t i_size = mColors.GetNodalCoordinates(0).size();
        std::size_t j_size = mColors.GetNodalCoordinates(1).size();
        std::size_t k_size = mColors.GetNodalCoordinates(2).size();
        while(!stack.empty()){
            auto& indices = stack.back();
            const auto i = std::get<0>(indices);
            const auto j = std::get<1>(indices);
            const auto k = std::get<2>(indices);
            stack.pop_back();
            mColors.GetElementalColor(i,j,k) = InsideColor;
            if(i > 0){
                if(mColors.GetElementalColor(i-1, j, k) == CellColor){
                    stack.emplace_back(i-1, j, k);
                }
            }
            if(j > 0){
                if(mColors.GetElementalColor(i, j-1, k) == CellColor){
                    stack.emplace_back(i, j-1, k);
                }
            }
            if(k > 0){
                if(mColors.GetElementalColor(i, j, k-1) == CellColor){
                    stack.emplace_back(i, j, k-1);
                }
            }
            if(i+1 < i_size){
                if(mColors.GetElementalColor(i+1, j, k) == CellColor){
                    stack.emplace_back(i+1, j, k);
                }
            }
            if(j+1 < j_size){
                if(mColors.GetElementalColor(i, j+1, k) == CellColor){
                    stack.emplace_back(i, j+1, k);
                }
            }
            if(k+1 < k_size){
                if(mColors.GetElementalColor(i, j, k+1) == CellColor){
                    stack.emplace_back(i, j, k+1);
                }
            }
        }
    }

    void VoxelMeshGeneratorModeler::ApplyColorIfOuterFace(double InterfaceColor, double CellColor, array_1d<std::size_t, 3> CellIndices){
        auto& faces_color = mColors.GetElementalFaceColor(CellIndices[0],CellIndices[1],CellIndices[2]);
        for(int i_direction = 0 ; i_direction < 3 ; i_direction++){
            if(CellIndices[i_direction] == 0) {  // It is the first cell and has outer face
                faces_color[i_direction] = InterfaceColor;
            }
            else{
                array_1d<std::size_t, 3> neighbor_cell_indices = CellIndices;
                neighbor_cell_indices[i_direction]--;
                if(mColors.GetElementalColor(neighbor_cell_indices[0],neighbor_cell_indices[1],neighbor_cell_indices[2]) != CellColor){
                    faces_color[i_direction] = InterfaceColor;
                }
            }
            if(CellIndices[i_direction] == (mKeyPlanes[i_direction].size() - 2)) {  // It is the last cell has outer face
                faces_color[i_direction + 3] = InterfaceColor;
            }
            else{
                array_1d<std::size_t, 3> neighbor_cell_indices = CellIndices;
                neighbor_cell_indices[i_direction]++;
                if(mColors.GetElementalColor(neighbor_cell_indices[0],neighbor_cell_indices[1],neighbor_cell_indices[2]) != CellColor){
                    faces_color[i_direction + 3] = InterfaceColor;
                }
            }
        }
    }

    void VoxelMeshGeneratorModeler::CheckGeometryType(const GeometryData::KratosGeometryType &geometry_type){
        KRATOS_ERROR_IF_NOT(geometry_type == GeometryData::KratosGeometryType::Kratos_Triangle3D3) << " Input entities must be of type Triangle3D3" << std::endl;
    }

    void VoxelMeshGeneratorModeler::SetStartIds(ModelPart& rTheVolumeModelPart){
        ModelPart& root_model_part = rTheVolumeModelPart.GetRootModelPart();
        if(root_model_part.NodesArray().empty()){
            mStartNodeId = 1;
        }
        else{
            mStartNodeId = root_model_part.NodesArray().back()->Id() + 1;
        }
        if(root_model_part.ElementsArray().empty()){
            mStartElementId = 1;
        }
        else{
            mStartElementId = root_model_part.ElementsArray().back()->Id() + 1;
        }
        if(root_model_part.ConditionsArray().empty()){
            mStartConditionId = 1;
        }
        else{
            mStartConditionId = root_model_part.ConditionsArray().back()->Id() + 1;
        }
    }

    void VoxelMeshGeneratorModeler::GenerateEntities(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters){
        Timer::Start("Voxel Mesh Entities Generation");

        for(auto parameters : EntityGeneratorParameters){

            Parameters default_parameters(R"(
                {
                    "type" : "elements_with_cell_color",
                    "model_part_name": "PLEASE SPECIFY IT",
                    "color": -1,
                    "properties_id": 1
                }  )");

            parameters.ValidateAndAssignDefaults(default_parameters);
            
            std::string model_part_name = parameters["model_part_name"].GetString();

            ModelPart& volume_part = CreateAndGetModelPart(model_part_name);
            SetStartIds(volume_part);

            KRATOS_ERROR_IF_NOT(parameters.Has("type")) << " The entities generator should have a type field" << std::endl;

            std::string type = parameters["type"].GetString();
            
            KRATOS_INFO("Modeler") << "Generate " << type << " in " << volume_part.Name()  << " model part" << std::endl;

            if(type == "elements_with_cell_color"){
                GenerateElementsWithCellColor(volume_part, parameters);
            }
            else if(type == "conditions_with_face_color"){
                GenerateConditionsWithFaceColor(volume_part, parameters);
            }
            else if(type == "tetrahedral_elements_with_cell_color")
            {
                GenerateTetrahedralElementsWithCellColor(volume_part, parameters);
            }
            else if(type == "triangular_conditions_with_face_color"){
                GenerateTriangularConditionsWithFaceColor(volume_part, parameters);
            }
            else {
                KRATOS_ERROR << "The etity generator type " << parameters["type"] << " is not supported. The supported types are  elements_with_cell_color" << std::endl;
            }
        }
        Timer::Stop("Voxel Mesh Entities Generation");
    }

    ModelPart& VoxelMeshGeneratorModeler::CreateAndGetModelPart(std::string const& FullName){
            std::istringstream iss(FullName);
            std::string name;
            std::getline(iss, name, '.');
            if(!mpModel->HasModelPart(name))
                mpModel->CreateModelPart(name);
            ModelPart* p_current_model_part = &mpModel->GetModelPart(name);
            while (std::getline(iss, name, '.'))
            {
                if(!p_current_model_part->HasSubModelPart(name)){
                    p_current_model_part->CreateSubModelPart(name);
                }
                p_current_model_part = &(p_current_model_part->GetSubModelPart(name));
            }
            
            return *p_current_model_part;
    }

    void VoxelMeshGeneratorModeler::GenerateElementsWithCellColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters){
         double inside_color = EntityGeneratorParameters["color"].GetDouble();
        std::size_t properties_id = EntityGeneratorParameters["properties_id"].GetInt();

        if(!rTheVolumeModelPart.HasProperties(properties_id)){
            rTheVolumeModelPart.CreateNewProperties(properties_id);
        }    
        Properties::Pointer p_properties = rTheVolumeModelPart.pGetProperties(properties_id);

        std::size_t cell_index = 0;
        array_1d<std::size_t, 3> number_of_cells;
        for(int i = 0 ; i < 3 ; i++){
            number_of_cells[i]=mKeyPlanes[i].size() - 1;
        }

        ModelPart::NodesContainerType new_nodes;
        ModelPart::ElementsContainerType new_elements;
        auto& r_prototype_element = KratosComponents<Element>::Get("Element3D8N");
            
    
        Element::NodesArrayType cell_nodes(8);
        for (std::size_t k = 0; k < number_of_cells[2]; k++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                    if(mColors.GetElementalColor(i,j,k) == inside_color){
                        cell_nodes(0) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j  , k);
                        cell_nodes(1) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j  , k);
                        cell_nodes(2) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j+1, k);
                        cell_nodes(3) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j+1, k);
                        cell_nodes(4) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j  , k+1);
                        cell_nodes(5) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j  , k+1);
                        cell_nodes(6) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j+1, k+1);
                        cell_nodes(7) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j+1, k+1);

                        //create the new element
                        Element::Pointer p_element = r_prototype_element.Create(mStartElementId + cell_index, cell_nodes, p_properties);
                        new_elements.push_back(p_element);
                        cell_index++;
                    }
                }
            }
        }
    
        rTheVolumeModelPart.AddNodes(new_nodes.begin(), new_nodes.end());
        rTheVolumeModelPart.AddElements(new_elements.begin(), new_elements.end());
   }

    void VoxelMeshGeneratorModeler::GenerateConditionsWithFaceColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters){
         double inside_color = EntityGeneratorParameters["color"].GetDouble();
        std::size_t properties_id = EntityGeneratorParameters["properties_id"].GetInt();

        const std::size_t x_offset[24]={0,0,0,0, 0,1,1,0, 0,0,1,1, 1,1,1,1, 0,0,1,1, 0,1,1,0};
        const std::size_t y_offset[24]={0,0,1,1, 0,0,0,0, 0,1,1,0, 0,1,1,0, 1,1,1,1, 0,0,1,1};
        const std::size_t z_offset[24]={0,1,1,0, 0,0,1,1, 0,0,0,0, 0,0,1,1, 0,1,1,0, 1,1,1,1};

        if(!rTheVolumeModelPart.HasProperties(properties_id)){
            rTheVolumeModelPart.CreateNewProperties(properties_id);
        }    
        Properties::Pointer p_properties = rTheVolumeModelPart.pGetProperties(properties_id);

        std::size_t condition_index = 0;
        array_1d<std::size_t, 3> number_of_cells;
        for(int i = 0 ; i < 3 ; i++){
            number_of_cells[i]=mKeyPlanes[i].size() - 1;
        }
        ModelPart::NodesContainerType new_nodes;
        ModelPart::ConditionsContainerType new_conditions;
        auto& r_prototype_condition = KratosComponents<Condition>::Get("SurfaceCondition3D4N");
    
        Element::NodesArrayType face_nodes(4);
        for (std::size_t k = 0; k < number_of_cells[2]; k++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                    auto& faces_color = mColors.GetElementalFaceColor(i,j,k);
                    for(std::size_t i_face = 0; i_face < 6; i_face++){
                        if(faces_color[i_face] == inside_color){
                            std::size_t base_index=i_face*4;
                            face_nodes(0) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index]  , j+y_offset[base_index]    , k+z_offset[base_index]  );
                            face_nodes(1) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+1], j+y_offset[base_index+1]  , k+z_offset[base_index+1]);
                            face_nodes(2) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+2], j+y_offset[base_index+2]  , k+z_offset[base_index+2]);
                            face_nodes(3) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+3], j+y_offset[base_index+3]  , k+z_offset[base_index+3]);
                            //create the new condition
                            Condition::Pointer p_condition = r_prototype_condition.Create(mStartConditionId + condition_index, face_nodes, p_properties);
                            new_conditions.push_back(p_condition);
                            condition_index++;
                        }
                    }
                }
            }
        }
        rTheVolumeModelPart.AddNodes(new_nodes.begin(), new_nodes.end());

        rTheVolumeModelPart.AddConditions(new_conditions.begin(), new_conditions.end());
   }

    void VoxelMeshGeneratorModeler::GenerateTriangularConditionsWithFaceColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters){
         double inside_color = EntityGeneratorParameters["color"].GetDouble();
        std::size_t properties_id = EntityGeneratorParameters["properties_id"].GetInt();

        const std::size_t x_offset[24]={0,0,0,0, 0,1,1,0, 0,0,1,1, 1,1,1,1, 0,0,1,1, 0,1,1,0};
        const std::size_t y_offset[24]={0,0,1,1, 0,0,0,0, 0,1,1,0, 0,1,1,0, 1,1,1,1, 0,0,1,1};
        const std::size_t z_offset[24]={0,1,1,0, 0,0,1,1, 0,0,0,0, 0,0,1,1, 0,1,1,0, 1,1,1,1};

        if(!rTheVolumeModelPart.HasProperties(properties_id)){
            rTheVolumeModelPart.CreateNewProperties(properties_id);
        }    
        Properties::Pointer p_properties = rTheVolumeModelPart.pGetProperties(properties_id);

        std::size_t condition_index = 0;
        array_1d<std::size_t, 3> number_of_cells;
        for(int i = 0 ; i < 3 ; i++){
            number_of_cells[i]=mKeyPlanes[i].size() - 1;
        }
        ModelPart::NodesContainerType new_nodes;
        ModelPart::ConditionsContainerType new_conditions;
        auto& r_prototype_condition = KratosComponents<Condition>::Get("SurfaceCondition3D3N");
    
        Element::NodesArrayType face_nodes_1(3);
        Element::NodesArrayType face_nodes_2(3);
        for (std::size_t k = 0; k < number_of_cells[2]; k++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                    auto& faces_color = mColors.GetElementalFaceColor(i,j,k);
                    for(std::size_t i_face = 0; i_face < 6; i_face++){
                        if(faces_color[i_face] == inside_color){
                            std::size_t base_index=i_face*4;
                            face_nodes_1(0) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index]  , j+y_offset[base_index]    , k+z_offset[base_index]  );
                            face_nodes_1(1) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+1], j+y_offset[base_index+1]  , k+z_offset[base_index+1]);
                            face_nodes_1(2) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+2], j+y_offset[base_index+2]  , k+z_offset[base_index+2]);

                            face_nodes_2(0) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index]  , j+y_offset[base_index]    , k+z_offset[base_index]  );
                            face_nodes_2(1) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+2], j+y_offset[base_index+2]  , k+z_offset[base_index+2]);
                            face_nodes_2(2) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+x_offset[base_index+3], j+y_offset[base_index+3]  , k+z_offset[base_index+3]);
                            //create the new conditions
                            Condition::Pointer p_condition_1 = r_prototype_condition.Create(mStartConditionId + condition_index, face_nodes_1, p_properties);
                            new_conditions.push_back(p_condition_1);
                            condition_index++;
                            Condition::Pointer p_condition_2 = r_prototype_condition.Create(mStartConditionId + condition_index, face_nodes_2, p_properties);
                            new_conditions.push_back(p_condition_2);
                            condition_index++;
                        }
                    }
                }
            }
        }
        rTheVolumeModelPart.AddNodes(new_nodes.begin(), new_nodes.end());

        rTheVolumeModelPart.AddConditions(new_conditions.begin(), new_conditions.end());
   }

   void VoxelMeshGeneratorModeler::GenerateTetrahedralElementsWithCellColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters){
        double inside_color = EntityGeneratorParameters["color"].GetDouble();
        std::size_t properties_id = EntityGeneratorParameters["properties_id"].GetInt();

        if(!rTheVolumeModelPart.HasProperties(properties_id)){
            rTheVolumeModelPart.CreateNewProperties(properties_id);
        }
        Properties::Pointer p_properties = rTheVolumeModelPart.pGetProperties(properties_id);

        std::size_t cell_index = 0;
        array_1d<std::size_t, 3> number_of_cells;
        for(int i = 0 ; i < 3 ; i++){
            number_of_cells[i]=mKeyPlanes[i].size() - 1;
        }
    
        ModelPart::NodesContainerType new_nodes;
        ModelPart::ElementsContainerType new_elements;
        // auto& r_prototype_element = KratosComponents<Element>::Get("Element3D8N");
            
        KRATOS_INFO("Modeler") << "        Generating elements ... " << std::endl;
        Element::NodesArrayType cell_nodes(8);
        for (std::size_t k = 0; k < number_of_cells[2]; k++) {
            for (std::size_t j = 0; j < number_of_cells[1]; j++) {
                for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                    if(mColors.GetElementalColor(i,j,k) == inside_color){
                        cell_nodes(0) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j  , k);
                        cell_nodes(1) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j  , k);
                        cell_nodes(2) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j+1, k);
                        cell_nodes(3) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j+1, k);
                        cell_nodes(4) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j  , k+1);
                        cell_nodes(5) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j  , k+1);
                        cell_nodes(6) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i+1, j+1, k+1);
                        cell_nodes(7) = GenerateOrRetriveNode(rTheVolumeModelPart, new_nodes, i  , j+1, k+1);

                        CreateTetrahedraInCell(rTheVolumeModelPart, cell_nodes,  mStartElementId+cell_index*6, p_properties);
                        
                        cell_index++;
                    }
                }
            }
        }
        rTheVolumeModelPart.AddNodes(new_nodes.begin(), new_nodes.end());
   }

   void VoxelMeshGeneratorModeler::CreateTetrahedraInCell(ModelPart& rTheVolumeModelPart,
                                                          Element::NodesArrayType & rCellNodes,
                                                          const std::size_t StartId,
                                                          Properties::Pointer & pProperties)
   {
        constexpr std::size_t nodes_per_tet = 4;
        constexpr std::size_t tets_per_cell = 6;
        using ConnectivityList = std::array<std::array<std::size_t, nodes_per_tet>, tets_per_cell>;
        static constexpr ConnectivityList local_connectivities {{ { 0,3,6,2 },
                                                                  { 3,6,7,0 },
                                                                  { 4,7,6,0 },
                                                                  { 0,4,5,6 },
                                                                  { 0,1,2,6 },
                                                                  { 1,5,6,0 } }};

        for(std::size_t i = 0; i < tets_per_cell; i++)
        {
            const auto & connectivty = local_connectivities[i];
            Element::NodesArrayType tet_nodes(nodes_per_tet);

            for(std::size_t n=0; n<nodes_per_tet; n++)
            {
                tet_nodes(n) = rCellNodes(connectivty[n]);
            }

            rTheVolumeModelPart.CreateNewElement("Element3D4N", StartId + i, tet_nodes, pProperties);
        }
   }

    Node<3>::Pointer VoxelMeshGeneratorModeler::GenerateOrRetriveNode(ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, std::size_t I, std::size_t J, std::size_t K){
        auto& nodal_data = mMeshingData.GetNodalData(I, J, K);
        if(nodal_data.IsCreated()){
            rThisNodes.push_back(nodal_data.pGetNode());
            return nodal_data.pGetNode();
        }

        Point global_coordinates = mColors.GetPoint(I,J,K);
        Node<3>::Pointer temp_node = Kratos::make_intrusive< Node<3> >( mStartNodeId++, global_coordinates[0], global_coordinates[1], global_coordinates[2]);
        // Giving model part's variables list to the node
        temp_node->SetSolutionStepVariablesList(rTheVolumeModelPart.pGetNodalSolutionStepVariablesList());

        //set buffer size
        temp_node->SetBufferSize(rTheVolumeModelPart.GetBufferSize());

        nodal_data.pSetNode(temp_node);
        rThisNodes.push_back(temp_node);
        return temp_node;
    }

} // namespace Kratos
