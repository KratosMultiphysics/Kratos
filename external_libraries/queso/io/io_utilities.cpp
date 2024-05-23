// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

//// STL includes
#include <map>
//// Project includes
#include "io/io_utilities.h"

namespace queso {

struct PointComparison {
    bool operator() (const PointType& rLhs, const PointType& rRhs) const {
        if( std::abs( rLhs[0] - rRhs[0] ) < SNAPTOL
            && std::abs(rLhs[1] - rRhs[1]) < SNAPTOL ){ // If equal
            return rLhs[2] < rRhs[2] - SNAPTOL;
        } else if ( std::abs( rLhs[0] - rRhs[0] ) < SNAPTOL ){
            return rLhs[1] < rRhs[1] - SNAPTOL; }
        else {
            return (rLhs)[0] <= (rRhs)[0] - SNAPTOL;
        }
    };
};

bool IO::WriteMeshToSTL(const TriangleMesh& rTriangleMesh,
                        const char* Filename,
                        const bool Binary){
    std::ofstream file;
    if(Binary)
        file.open(Filename, std::ios::out | std::ios::binary);
    else
        file.open(Filename);

    if(!file.good()){
        std::cerr << "Warning :: IO::WriteMeshToSTL :: Could not open file: " << Filename << ".\n";
        return false;
    }

    const IndexType num_triangles = rTriangleMesh.NumOfTriangles();
    if(Binary) {
        file << "FileType: Binary                                                                ";
        file.write(reinterpret_cast<const char *>(&static_cast<const int&>(num_triangles)), sizeof(static_cast<const int&>(num_triangles)));

        for(IndexType triangle_id = 0; triangle_id < num_triangles; ++triangle_id) {
            const auto& p1 = rTriangleMesh.P1(triangle_id);
            const auto& p2 = rTriangleMesh.P2(triangle_id);
            const auto& p3 = rTriangleMesh.P3(triangle_id);

            const auto& normal = rTriangleMesh.Normal(triangle_id);

            const float coords[12] = { static_cast<float>(normal[0]), static_cast<float>(normal[1]), static_cast<float>(normal[2]),
                                        static_cast<float>(p1[0]), static_cast<float>(p1[1]), static_cast<float>(p1[2]),
                                        static_cast<float>(p2[0]), static_cast<float>(p2[1]), static_cast<float>(p2[2]),
                                        static_cast<float>(p3[0]), static_cast<float>(p3[1]), static_cast<float>(p3[2]) };

            for(int i=0; i<12; ++i){
                file.write(reinterpret_cast<const char *>(&coords[i]), sizeof(coords[i]));
            }
            file << "  ";
        }
    }
    else
    {
        file << "solid\n";
        for(IndexType triangle_id = 0; triangle_id < num_triangles; ++triangle_id) {
            const auto& p1 = rTriangleMesh.P1(triangle_id);
            const auto& p2 = rTriangleMesh.P2(triangle_id);
            const auto& p3 = rTriangleMesh.P3(triangle_id);

            const auto& normal = rTriangleMesh.Normal(triangle_id);

            file << "facet normal " << normal[0] << ' ' << normal[1] << ' ' << normal[2] << "\nouter loop\n";
            file << "vertex " << p1[0] << ' ' << p1[1] << ' ' << p1[2] << ' ' << "\n";
            file << "vertex " << p2[0] << ' ' << p2[1] << ' ' << p2[2] << ' ' << "\n";
            file << "vertex " << p3[0] << ' ' << p3[1] << ' ' << p3[2] << ' ' << "\n";
            file << "endloop\nendfacet\n";
        }
        file << "endsolid"<<std::endl;
    }
    file.close();

    return true;
}

bool IO::ReadMeshFromSTL(TriangleMesh& rTriangleMesh,
                         const char* Filename){

    // Open file
    if( STLIsInASCIIFormat(Filename) ) { // If the first 5 characters are "solid"
        return ReadMeshFromSTL_Ascii(rTriangleMesh, Filename);
    } else {
        return ReadMeshFromSTL_Binary(rTriangleMesh, Filename);
    }
}

bool IO::WriteMeshToVTK(const TriangleMesh& rTriangleMesh,
                        const char* Filename,
                        const bool Binary) {

    const SizeType num_elements = rTriangleMesh.NumOfTriangles();
    const SizeType num_points = rTriangleMesh.NumOfVertices();

    std::ofstream file;
    if(Binary)
        file.open(Filename, std::ios::out | std::ios::binary);
    else
        file.open(Filename);

    file << "# vtk DataFile Version 4.1" << std::endl;
    file << "vtk output" << std::endl;
    if(Binary)
        file << "BINARY"<< std::endl;
    else
        file << "ASCII"<< std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << num_points << " double" << std::endl;

    const auto& r_vertices = rTriangleMesh.GetVertices();
    const auto v_it_begin = r_vertices.begin();
    for(IndexType i = 0; i < num_points; ++i) {
        const auto v_it = v_it_begin+i;
        if( Binary ){
            double rx = (*v_it)[0];
            double ry = (*v_it)[1];
            double rz = (*v_it)[2];

            WriteBinary(file, rx);
            WriteBinary(file, ry);
            WriteBinary(file, rz);
        }
        else {
            file << (*v_it)[0] << ' ' << (*v_it)[1] << ' ' << (*v_it)[2] << std::endl;
        }
    }
    file << std::endl;

    // Write Cells
    file << "Cells " << num_elements << " " << num_elements*4 << std::endl;

    for(IndexType i = 0; i < num_elements; ++i) {
        if( Binary ){
            int k = 3;
            WriteBinary(file, k);
            for( auto id : rTriangleMesh.VertexIds(i) ){
                k = id;
                WriteBinary(file, k);
            }
        }
        else {
            file << 3;
            for( auto id : rTriangleMesh.VertexIds(i) ){
                file << ' ' << id;
            }
            file << std::endl;
        }
    }
    file << std::endl;

    file << "CELL_TYPES " << num_elements << std::endl;
    for( IndexType i = 0; i < num_elements; ++i){
        if( Binary ){
            int k = 5;
            WriteBinary(file, k);
        }
        else {
            file << 5 << std::endl;
        }
    }
    file << std::endl;
    file.close();

    return true;
}

bool IO::WriteDisplacementToVTK(const std::vector<Vector3d>& rDisplacement,
                                const char* Filename,
                                const bool Binary){

    const SizeType num_points = rDisplacement.size();

    std::ofstream file;
    if(Binary)
        file.open(Filename, std::ios::app | std::ios::binary);
    else
        file.open(Filename);

    file << "POINT_DATA " << num_points << std::endl;
    file << "VECTORS Displacement double" << std::endl;
    for(IndexType i = 0; i < num_points; ++i){
        if( Binary ){
            double rw1 = rDisplacement[i][0];
            WriteBinary(file, rw1);
            double rw2 = rDisplacement[i][1];
            WriteBinary(file, rw2);
            double rw3 = rDisplacement[i][2];
            WriteBinary(file, rw3);
        }
        else {
            std::cerr << "Warning :: IO::DisplacementToVTK :: Ascii export not implemented yet. \n";
            return false;
        }
    }
    file << std::endl;
    file.close();

    return true;
}

bool IO::WriteElementsToVTK(const ElementContainer& rElementContainer, //PolygonMesh
                            const char* Filename,
                            const bool Binary){

    const SizeType num_elements = rElementContainer.size();

    std::ofstream file;
    if(Binary)
        file.open(Filename, std::ios::out | std::ios::binary);
    else
        file.open(Filename);

    file << "# vtk DataFile Version 4.1" << std::endl;
    file << "vtk output" << std::endl;
    if(Binary)
        file << "BINARY"<< std::endl;
    else
        file << "ASCII"<< std::endl;

    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << num_elements*8 << " double" << std::endl;

    const auto begin_el_itr = rElementContainer.begin();
    for( IndexType i = 0; i < rElementContainer.size(); ++i){
        const auto& el_itr = *(begin_el_itr + i);
        const auto& lower_point = el_itr->GetBoundsXYZ().first;
        const auto& upper_point = el_itr->GetBoundsXYZ().second;

        if( Binary ){
            double rx0 = lower_point[0];
            SwapEnd(rx0);
            double rx1 = upper_point[0];
            SwapEnd(rx1);
            double ry0 = lower_point[1];
            SwapEnd(ry0);
            double ry1 = upper_point[1];
            SwapEnd(ry1);
            double rz0 = lower_point[2];
            SwapEnd(rz0);
            double rz1 = upper_point[2];
            SwapEnd(rz1);

            file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz0), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry0), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx1), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz1), sizeof(double));

            file.write(reinterpret_cast<char*>(&rx0), sizeof(double));
            file.write(reinterpret_cast<char*>(&ry1), sizeof(double));
            file.write(reinterpret_cast<char*>(&rz1), sizeof(double));
        }
        else {
            file << lower_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
            file << upper_point[0] << ' ' << lower_point[1] << ' ' << lower_point[2] << std::endl;
            file << upper_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
            file << lower_point[0] << ' ' << upper_point[1] << ' ' << lower_point[2] << std::endl;
            file << lower_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
            file << upper_point[0] << ' ' << lower_point[1] << ' ' << upper_point[2] << std::endl;
            file << upper_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
            file << lower_point[0] << ' ' << upper_point[1] << ' ' << upper_point[2] << std::endl;
        }
    }
    file << std::endl;
    // Write Cells
    file << "Cells " << num_elements << " " << num_elements*9 << std::endl;
    for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
        if( Binary ){
            int k = 8;
            WriteBinary(file, k);
        for( int j = 0; j < 8; ++j){
            k = 8*i+j;
            WriteBinary(file, k);
        }
        }
        else {
            file << 8 << ' ' << 8*i     << ' ' << 8*i + 1 << ' ' << 8*i + 2 << ' ' << 8*i + 3
                        << ' ' << 8*i + 4 << ' ' << 8*i + 5 << ' ' << 8*i + 6 << ' ' << 8*i + 7 << std::endl;
        }
    }
    file << std::endl;

    file << "CELL_TYPES " << rElementContainer.size() << std::endl;
    for( int i = 0; i < static_cast<int>(rElementContainer.size()); ++i){
        if( Binary ){
            int k = 12;
            WriteBinary(file, k);
        }
        else {
            file << 12 << std::endl;
        }
    }
    file << std::endl;
    file.close();

    return true;
}


bool IO::WritePointsToVTK(const ElementContainer& rElementContainer,
                          const char* type,
                          const char* Filename,
                          const bool Binary){

    QuESo_ERROR_IF( std::string(type) != "All") << "Only integration points option 'All' is available. \n";

    const auto& p_points = rElementContainer.pGetPoints(type);
    const IndexType num_points = p_points->size();
    const IndexType num_elements = p_points->size();

    std::ofstream file;
    if(Binary)
        file.open(Filename, std::ios::out | std::ios::binary);
    else
        file.open(Filename);

    file << "# vtk DataFile Version 4.1" << std::endl;
    file << "vtk output" << std::endl;
    if(Binary)
        file << "BINARY"<< std::endl;
    else
        file << "ASCII"<< std::endl;


    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << num_points << " double" << std::endl;

    const auto el_it_ptr_begin = rElementContainer.begin();
    for( IndexType i = 0; i < rElementContainer.size(); ++i){
        const auto& el_ptr = (*(el_it_ptr_begin + i));
        const auto& r_points = el_ptr->GetIntegrationPoints();
        for( const auto& r_point : r_points ){
            auto point_global = el_ptr->PointFromParamToGlobal(r_point);

            if( Binary ){
                WriteBinary(file, point_global[0]);
                WriteBinary(file, point_global[1]);
                WriteBinary(file, point_global[2]);
            }
            else {
                file << point_global[0] << ' ' << point_global[1] << ' ' << point_global[2] << std::endl;
            }
        }

    }
    file << std::endl;

    //Write Cells
    file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
    for( IndexType i = 0; i < num_elements; ++i){
        if( Binary ){
            int k = 1;
            WriteBinary(file, k);
            k = i;
            WriteBinary(file, k);
        }
        else {
            file << 1 << ' ' << i << std::endl;
        }
    }
    file << std::endl;

    file << "CELL_TYPES " << num_elements << std::endl;
    for( IndexType i = 0; i < num_elements; ++i){
        if( Binary ){
            int k = 1;
            WriteBinary(file, k);
        }
        else {
            file << 1 << std::endl;
        }
    }
    file << std::endl;

    file << "POINT_DATA " << num_points << std::endl;
    file << "SCALARS Weights double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for( IndexType i = 0; i < rElementContainer.size(); ++i){
        const auto& el_ptr = (*(el_it_ptr_begin + i));
        const auto& points = el_ptr->GetIntegrationPoints();
        for( const auto& point : points ){
            if( Binary ){
                double rw = point.Weight();
                WriteBinary(file, rw);
            }
            else {
                file << point.Weight() << std::endl;
            }
        }
    }
    file << std::endl;
    file.close();

    return true;
}


template<typename Type>
bool IO::WritePointsToVTK(const std::vector<Type>& rPoints,
                          const char* Filename,
                          const bool Binary){

    const auto begin_points_it_ptr = rPoints.begin();
    const IndexType num_points = rPoints.size();
    const IndexType num_elements = rPoints.size();

    std::ofstream file;
    if(Binary)
        file.open(Filename, std::ios::out | std::ios::binary);
    else
        file.open(Filename);

    file << "# vtk DataFile Version 4.1" << std::endl;
    file << "vtk output" << std::endl;
    if(Binary)
        file << "BINARY"<< std::endl;
    else
        file << "ASCII"<< std::endl;


    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    file << "POINTS " << num_points << " double" << std::endl;

    for(IndexType i = 0; i < num_points; ++i){
        auto points_it = (begin_points_it_ptr + i);
        if( Binary ){
            // Make sure to create copy of points. "WriteBinary" will change them.
            auto p_x = (*points_it)[0];
            auto p_y = (*points_it)[1];
            auto p_z = (*points_it)[2];
            WriteBinary(file, p_x);
            WriteBinary(file, p_y);
            WriteBinary(file, p_z);
        }
        else {
            file << (*points_it)[0] << ' ' << (*points_it)[1] << ' ' << (*points_it)[2] << std::endl;
        }
    }
    file << std::endl;

    //Write Cells
    file << "Cells " << num_elements << " " << num_elements*2 << std::endl;
    for( IndexType i = 0; i < num_elements; ++i){
        if( Binary ){
            int k = 1;
            WriteBinary(file, k);
            k = i;
            WriteBinary(file, k);
        }
        else {
            file << 1 << ' ' << i << std::endl;
        }
    }
    file << std::endl;

    file << "CELL_TYPES " << num_elements << std::endl;
    for( IndexType i = 0; i < num_elements; ++i){
        if( Binary ){
            int k = 1;
            WriteBinary(file, k);
        }
        else {
        file << 1 << std::endl;
        }
    }
    file << std::endl;

    file << "POINT_DATA " << num_points << std::endl;
    file << "SCALARS Weights double 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for(IndexType i = 0; i < num_points; ++i){
        auto points_it = (begin_points_it_ptr + i);

        if( Binary ){
            double rw = points_it->Weight();
            WriteBinary(file, rw);
        }
        else {
            file << points_it->Weight() << std::endl;
        }
    }
    file << std::endl;
    file.close();

    return true;
}

// Instantiation
template bool IO::WritePointsToVTK<IntegrationPoint>(const std::vector<IntegrationPoint>& rPoints,
                                                    const char* Filename,
                                                    const bool Binary);

template bool IO::WritePointsToVTK<BoundaryIntegrationPoint>(const std::vector<BoundaryIntegrationPoint>& rPoints,
                                                            const char* Filename,
                                                            const bool Binary);

////// Private member functions //////

bool IO::STLIsInASCIIFormat(const char* Filename) {
    std::ifstream file(Filename, std::ios::in);

    if( !file.good() ){
        QuESo_ERROR << "Couldnt parse file: " << Filename << std::endl;
    }

    char chars [256];
    file.read(chars, 256);
    std::string message(chars, file.gcount());
    transform(message.begin(), message.end(), message.begin(), ::tolower);
    return message.find ("solid")  != std::string::npos &&
           message.find ("normal") != std::string::npos &&
           message.find ("facet")  != std::string::npos &&
           message.find ("\n")     != std::string::npos;
}


bool IO::ReadMeshFromSTL_Ascii(TriangleMesh& rTriangleMesh,
                                const char* Filename){
    // Open file
    std::ifstream file(Filename);
    if( !file.good() ) {
        QuESo_ERROR << "Couldnt parse file: " << Filename << ".\n";
        return false;
    }

    // Initialize map
    int index = 0;
    std::map<Vector3d, int, PointComparison> index_map;

    rTriangleMesh.Clear();
    rTriangleMesh.Reserve(100000);

    std::string message;
    std::getline(file, message); // Ignore "Solid"

    while( !file.eof() || file.fail() ) {
        std::getline(file, message); // Read normal (normals are computed from vertices later.)
        if(  message.find("endsolid") !=  std::string::npos ){
            break;
        }
        std::string token;
        std::getline(file, message); // Ignore outer loop

        Vector3i triangle{};
        std::array<PointType, 3> vertices;
        for( IndexType i = 0; i < 3; i++){
        std::getline(file, message); // Read vertex

        IndexType j = 0;
        Vector3d vertex{};
        std::stringstream ss(message);
        while( std::getline(ss, token, ' ') ){
            token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
            if( token.size() > 0 && token.find("vertex") == std::string::npos ) {
                float value = 0.0;
                try {
                    value = std::stof(token);
                } catch (const std::invalid_argument& e) {
                    std::cout << "Invalid argument: " << e.what() << 'n';
                } catch (const std::out_of_range& e) {
                    double test_value = std::stod(token);
                    if( std::abs(test_value) > 1e-16 ){
                        std::cout << "IO::ReadMeshFromSTL_Ascii :: Out-of-range-value of vertex: " << token << " is set to zero.\n";
                    }
                    value = 0.0;
                }
                vertex[j] = value;
                ++j;
            }
        }

        // Map is used to ensure unique vertices. Note that STL does not reuse vertices.
        auto index_map_iterator = index_map.insert(std::make_pair(vertex, -1)).first;
        if(index_map_iterator->second == -1) {
            triangle[i] = index;
            index_map_iterator->second = index++;
            rTriangleMesh.AddVertex(vertex);
            vertices[i] = vertex;
        }
        else {
            triangle[i] = index_map_iterator->second;
            vertices[i] = index_map_iterator->first;
        }
        }

        // Add Triangle
        rTriangleMesh.AddTriangle(triangle);

        // Note: STL are often given in single precision. Therefore, we have to compute the normals based on
        // the given vertices.
        const PointType normal = TriangleMesh::Normal(vertices[0], vertices[1], vertices[2]);

        rTriangleMesh.AddNormal( normal );

        std::getline(file, message); // Ignore endloop
        std::getline(file, message); // Ignore endfacet
    }
    file.close();
    return rTriangleMesh.Check();
}

bool IO::ReadMeshFromSTL_Binary(TriangleMesh& rTriangleMesh,
                                const char* Filename){

    // Open file
    std::ifstream file(Filename, std::ios::binary);

    if( !file.good() ) {
        QuESo_ERROR << "Couldnt parse file: " << Filename << ".\n";
        return false;
    }

    // Ignore the first 80 chars of the header
    int position = 0;
    char message;
    std::string test_binary_ascii{};
    while(position < 80) {
        file.read(reinterpret_cast<char*>(&message), sizeof(message));
        if( !file.good() )
            break;
        ++position;
    }

    if( test_binary_ascii == "solid" ) { // If the first 5 characters are "solid"
        QuESo_ERROR << "Warning :: File seems to be in Ascii format :: Please use IO::ReadMeshFromSTL_Binary.\n";
        return false;
    }

    if(position != 80) {
        QuESo_ERROR << "File " << Filename << " is empty.\n";
        return false;
    }

    int index = 0;
    std::map<Vector3d, int, PointComparison> index_map;

    // Read number of triangles
    unsigned int num_triangles;
    if(!(file.read(reinterpret_cast<char*>(&num_triangles), sizeof(num_triangles)))) {
        QuESo_ERROR << "Couldnt read number of triangles. \n";
        return false;
    }
    rTriangleMesh.Clear();
    rTriangleMesh.Reserve(num_triangles);

    // Loop over all triangles
    for(IndexType i=0; i<num_triangles; ++i) {
        // Read normals
        float normal_tmp[3]; // Normals are ignored. They are computed from the vertices later.
        if(!(file.read((char*)(&normal_tmp[0]), sizeof(normal_tmp[0]))) ||
                !(file.read((char*)(&normal_tmp[1]), sizeof(normal_tmp[1]))) ||
                !(file.read((char*)(&normal_tmp[2]), sizeof(normal_tmp[2])))) {
            QuESo_ERROR << "Couldnt read normals. \n";
            return false;
        }

        // Read triangles and vertices. Each vertex is read seperately.
        Vector3i triangle{};
        std::array<PointType, 3> vertices;
        for(int j=0; j<3; ++j) {
            float x,y,z;
            if(!(file.read((char*)(&x), sizeof(x))) ||
                    !(file.read((char*)(&y), sizeof(y))) ||
                    !(file.read((char*)(&z), sizeof(z)))) {
                QuESo_ERROR << "Couldnt read coordinates. \n";
                return false;
            }
            Vector3d vertex = {x,y,z};

            // Map is used to ensure unique vertices. Note that STL does not reuse vertices.
            auto index_map_iterator = index_map.insert(std::make_pair(vertex, -1)).first;
            if(index_map_iterator->second == -1) {
                triangle[j] = index;
                index_map_iterator->second = index++;
                rTriangleMesh.AddVertex(vertex);
                vertices[j] = vertex;
            }
            else {
                triangle[j] = index_map_iterator->second;
                vertices[j] = index_map_iterator->first;
            }
        }

        // Add triangle
        rTriangleMesh.AddTriangle(triangle);

        // Note: STL are often given in single precision. Therefore, we have to compute the normals based on
        // the given vertices.
        const PointType normal = TriangleMesh::Normal(vertices[0], vertices[1], vertices[2]);

        rTriangleMesh.AddNormal(normal);

        // Read so-called attribute byte count and ignore it
        char c;
        if(!(file.read((char*)(&c), sizeof(c))) ||
            !(file.read((char*)(&c), sizeof(c)))) {

            QuESo_ERROR << "Couldnt read attribute byte count.\n";
            return false;
        }
    }
    file.close();
    return rTriangleMesh.Check();
}

} // End namespace queso


