//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "cad_io_modeler.h"
#include "input_output/cad_json_input.h"
#include "input_output/cad_json_output.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void CadIoModeler::SetupGeometryModel()
    {
        // Read the refinements.iga.json
        const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");

        /// MODIFIED
        // Check if true boundary is present in the problem
        std::ifstream infile("txt_files/true_points.txt");
        bool file_true_points_exists = infile.good(); 
        if (file_true_points_exists==true){
            // Create the snakes coordiantes
            CreateTheSnakeCoordinates(refinements_parameters);
        }

        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part_name\" in CadIoModeler Parameters." << std::endl;
        const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
        ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
            ? mpModel->GetModelPart(cad_model_part_name)
            : mpModel->CreateModelPart(cad_model_part_name);

        const std::string DataFileName = mParameters.Has("geometry_file_name")
            ? mParameters["geometry_file_name"].GetString()
            : "geometry.cad.json";

        KRATOS_INFO_IF("::[CadIoModeler]::", mEchoLevel > 0) << "Importing Cad Model from: " << DataFileName << std::endl;

        CadJsonInput<Node, Point>(
            DataFileName, mEchoLevel).ReadModelPart(cad_model_part);
    }

    void CadIoModeler::SetupModelPart()
    {
        if (mParameters.Has("output_geometry_file_name")) {
            std::string DataFileName = mParameters["output_geometry_file_name"].GetString();

            const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
            ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
                ? mpModel->GetModelPart(cad_model_part_name)
                : mpModel->CreateModelPart(cad_model_part_name);

            std::string output_file_text;
            CadJsonOutput::GetCadJsonOutput(cad_model_part, output_file_text, mEchoLevel);

            std::ofstream output_file(DataFileName);
            output_file << output_file_text;
            output_file.close();
        }
    }
    ///@}

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters CadIoModeler::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".iga.json") != 0)
            ? rDataFileName + ".iga.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Physics fil: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }









    /// MODIFIED
    void CadIoModeler::CreateTheSnakeCoordinates(const Parameters refinements_parameters) {
        int insert_nb_per_span_u = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        int insert_nb_per_span_v = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_v"].GetInt();

        std::vector<double> knots_u; 
        std::vector<double> knots_v;  
        knots_u.push_back(0.0) ;
        knots_v.push_back(0.0) ;
        double initial = 0.0; 
        double total = 2.0; 
        double knot_step_u = total/(insert_nb_per_span_u+1) ;
        double knot_step_v = total/(insert_nb_per_span_v+1) ;
        for (IndexType j = 1; j < insert_nb_per_span_u + 1; ++j) {
            knots_u.push_back(initial+ knot_step_u * j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_v + 1; ++j) {
            knots_v.push_back(initial+ knot_step_v *j);
        }
        knots_u.push_back(2.0) ;
        knots_v.push_back(2.0) ;
        KRATOS_WATCH(knots_u)

        // Create the matrix of active/inactive knot spans
        std::vector<std::vector<int>> knot_spans_available(insert_nb_per_span_v+1, std::vector<int>(insert_nb_per_span_u+1));

        // Read the point from the mdpa,read 
        std::ifstream file("txt_files/true_points.txt");
        std::vector<double> x_true_boundary;
        std::vector<double> y_true_boundary;
        double x, y;
        while (file >> x >> y) {
            x_true_boundary.push_back(x);
            y_true_boundary.push_back(y);
        }
        file.close();
        const int numPoints = x_true_boundary.size();

        // Intersections
        std::vector<double> x_intersection;
        std::vector<double> y_intersection;

        // Addition for the last segment of the true boundary (for intersection process)
        x_true_boundary.push_back(x_true_boundary[0]);
        y_true_boundary.push_back(y_true_boundary[0]);

        // internal/external surrogate boundary
        std::vector<double> x_advanced_surrogate_boundary;
        std::vector<double> y_advanced_surrogate_boundary;
        x_advanced_surrogate_boundary.push_back(0.0) ;  // Fictitous for the algorithm
        y_advanced_surrogate_boundary.push_back(0.0) ;



        for (int i = 0; i < numPoints; i = i+2) {
            // Initizalize the pair of true_boundary point
            double x_true_boundary1 = x_true_boundary[i];
            double y_true_boundary1 = y_true_boundary[i];
            double x_true_boundary2 = x_true_boundary[i+1];
            double y_true_boundary2 = y_true_boundary[i+1];
            
            //// (OLD) Search in which knot span belongs --> knot_spans_available (OLD)
            // int number_knot_span_u = 0;
            // int number_knot_span_v = 0;
            // for (int j=0; j< knots_u.size() - 1; j++) {
            //     if (x_true_boundary1 > knots_u[j] && x_true_boundary1 <= knots_u[j+1]) {number_knot_span_u = j; } ;
            // } 
            // for (int j=0; j< knots_v.size() - 1; j++) {
            //     if (y_true_boundary1 > knots_v[j] && y_true_boundary1 <= knots_v[j+1]) {number_knot_span_v = j; } ;
            // } 
            // knot_spans_available[number_knot_span_v][number_knot_span_u] = 1;

            // Find the intersections of the skin boundary with the knot values
            int knot_span_u_1st_point = x_true_boundary1 / knot_step_u ;
            int knot_span_u_2nd_point = x_true_boundary2 / knot_step_u ;
            int knot_span_v_1st_point = y_true_boundary1 / knot_step_v ;
            int knot_span_v_2nd_point = y_true_boundary2 / knot_step_v ;

            // candidate x- and y-coordinate of the "advanced_surrogate_boundary"
            double x_candidate_advanced_surrogate_boundary;
            double y_candidate_advanced_surrogate_boundary;

            if (knot_span_u_1st_point != knot_span_u_2nd_point || knot_span_v_1st_point != knot_span_v_2nd_point) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY
                // Compute the line connecting true_points (i, i+1) -> y = mx+q
                double m = (y_true_boundary2-y_true_boundary1) / (x_true_boundary2 - x_true_boundary1) ;
                double q = y_true_boundary1 - m * x_true_boundary1 ;
                // Check if the true boundary crosses an u or a v knot value
                if (knot_span_u_1st_point != knot_span_u_2nd_point) { // u knot value is crossed
                    // KRATOS_WATCH('// u knot value is crossed')
                    double x_intersection_value = std::max(knot_span_u_1st_point, knot_span_u_2nd_point) * knot_step_u ;
                    double y_intersection_value = m * x_intersection_value + q ;
                    x_intersection.push_back(x_intersection_value);
                    y_intersection.push_back(y_intersection_value);
                    // Find the "knot_spans_available" using the intersection
                    knot_spans_available[knot_span_v_1st_point][knot_span_u_1st_point] = 1;
                    knot_spans_available[knot_span_v_1st_point][knot_span_u_2nd_point] = 1;
                    // Save the candidate for the advanced_surrogate_boundary
                    if (std::abs(y_intersection_value-knot_span_v_1st_point*knot_step_v) > std::abs(y_intersection_value-(knot_span_v_1st_point+1)*knot_step_v)) {
                        x_candidate_advanced_surrogate_boundary = x_intersection_value ;
                        y_candidate_advanced_surrogate_boundary = (knot_span_v_1st_point+1)*knot_step_v ; 
                    }
                    else { 
                        x_candidate_advanced_surrogate_boundary = x_intersection_value ;
                        y_candidate_advanced_surrogate_boundary = knot_span_v_1st_point*knot_step_v ; 
                    }
                }
                if (knot_span_v_1st_point != knot_span_v_2nd_point) { // v knot value is crossed
                    // KRATOS_WATCH('// v knot value is crossed')
                    double y_intersection_value = std::max(knot_span_v_1st_point, knot_span_v_2nd_point) * knot_step_v ;
                    double x_intersection_value = 0; // Initialize
                    // m might be inf !! We need the following if statement
                    if (x_true_boundary2 - x_true_boundary1 != 0)  {x_intersection_value = (y_intersection_value - q) / m ;}
                    else {x_intersection_value = x_true_boundary1 ;}
                    y_intersection.push_back(y_intersection_value);
                    x_intersection.push_back(x_intersection_value);
                    // KRATOS_WATCH(x_intersection_value)
                    // KRATOS_WATCH(y_intersection_value)
                    // Find the "knot_spans_available" using the intersection (Snake_coordinate classic -> External Boundary)
                    knot_spans_available[knot_span_v_1st_point][knot_span_u_1st_point] = 1;
                    knot_spans_available[knot_span_v_2nd_point][knot_span_u_1st_point] = 1;
                    // Save the closest knot value
                    if (std::abs(x_intersection_value-knot_span_u_1st_point*knot_step_u) > std::abs(x_intersection_value-(knot_span_u_1st_point+1)*knot_step_u)) {
                        x_candidate_advanced_surrogate_boundary = (knot_span_u_1st_point+1)*knot_step_u ;
                        y_candidate_advanced_surrogate_boundary = y_intersection_value ;
                    }
                    else { 
                        x_candidate_advanced_surrogate_boundary = knot_span_u_1st_point*knot_step_u ;
                        y_candidate_advanced_surrogate_boundary = y_intersection_value ; 
                    }
                }
                // Check if we are jumping some cut knot spans. If yes we have to take a finer grid. (We can improve here...)
                if (std::abs(knot_span_v_1st_point-knot_span_v_2nd_point) > 1 || std::abs(knot_span_u_1st_point-knot_span_u_2nd_point) > 1 ) {
                    KRATOS_WATCH('Your true boundary is too coarse, please take a finer one')
                    exit(0);
                }
                // Add the x- and y-candidate_advanced_surrogate_boundary if the point has not been added yet
                if (x_candidate_advanced_surrogate_boundary != x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] || 
                    y_candidate_advanced_surrogate_boundary != y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1]) {
                    // Check if we are in the case of a "diagonal" surrogate boundary and in the case add an intermidiate advanced_surrogate_boundary point
                    if (x_candidate_advanced_surrogate_boundary != x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] &&
                            y_candidate_advanced_surrogate_boundary != y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1] && 
                            x_advanced_surrogate_boundary.size() > 1) {
                        // KRATOS_WATCH('Diagional surrogate boundary')
                        // Let's identify the two candidate points, only one will be added to the advanced_surrogate_boundary
                        double x1_candidate = x_candidate_advanced_surrogate_boundary;
                        double y1_candidate = y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1];
                        double x2_candidate = x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1];
                        double y2_candidate = y_candidate_advanced_surrogate_boundary;
                        // KRATOS_WATCH(x1_candidate)
                        // KRATOS_WATCH(y1_candidate)
                        // KRATOS_WATCH(x2_candidate)
                        // KRATOS_WATCH(y2_candidate)
                        // Add the one which is closer to both the intersections : [x_intersection.size() - 1] and [x_intersection.size() - 2]
                        if (      std::abs((x1_candidate)*(x1_candidate)-x_intersection[x_intersection.size() - 1]*x_intersection[x_intersection.size() - 1] + (y1_candidate)*(y1_candidate)-y_intersection[y_intersection.size() - 1]*y_intersection[y_intersection.size() - 1] 
                                + (x1_candidate)*(x1_candidate)-x_intersection[x_intersection.size() - 2]*x_intersection[x_intersection.size() - 2] + (y1_candidate)*(y1_candidate)-y_intersection[y_intersection.size() - 2]*y_intersection[y_intersection.size() - 2]  ) 
                            <     std::abs((x2_candidate)*(x2_candidate)-x_intersection[x_intersection.size() - 1]*x_intersection[x_intersection.size() - 1] + (y2_candidate)*(y2_candidate)-y_intersection[y_intersection.size() - 1]*y_intersection[y_intersection.size() - 1] 
                                + (x2_candidate)*(x2_candidate)-x_intersection[x_intersection.size() - 2]*x_intersection[x_intersection.size() - 2] + (y2_candidate)*(y2_candidate)-y_intersection[y_intersection.size() - 2]*y_intersection[y_intersection.size() - 2])    ) {
                            x_advanced_surrogate_boundary.push_back(x1_candidate) ;
                            y_advanced_surrogate_boundary.push_back(y1_candidate) ;
                        }
                        else {
                            x_advanced_surrogate_boundary.push_back(x2_candidate) ;
                            y_advanced_surrogate_boundary.push_back(y2_candidate) ;
                        }
                    }
                    x_advanced_surrogate_boundary.push_back(x_candidate_advanced_surrogate_boundary) ;
                    y_advanced_surrogate_boundary.push_back(y_candidate_advanced_surrogate_boundary) ;
                    // KRATOS_WATCH(x_candidate_advanced_surrogate_boundary)
                    // KRATOS_WATCH(y_candidate_advanced_surrogate_boundary)
                }
            }
        }
        // Case of diagonal surrogate boundary between last and first point of advanced_surrogate_boundary
        if (x_advanced_surrogate_boundary[1] != x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] &&
            y_advanced_surrogate_boundary[1] != y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1] ) {
            // KRATOS_WATCH('Diagional surrogate boundary between last and first point')
            // Let's identify the two candidate points, only one will be added to the advanced_surrogate_boundary
            double x1_candidate = x_advanced_surrogate_boundary[1];
            double y1_candidate = y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1];
            double x2_candidate = x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1];
            double y2_candidate = y_advanced_surrogate_boundary[1];
            // Add the one which is closer to both the intersections : [x_intersection.size() - 1] and [x_intersection.size() - 2]
            if (      std::abs((x1_candidate)*(x1_candidate)-x_intersection[x_intersection.size() - 1]*x_intersection[x_intersection.size() - 1] + (y1_candidate)*(y1_candidate)-y_intersection[y_intersection.size() - 1]*y_intersection[y_intersection.size() - 1] 
                    + (x1_candidate)*(x1_candidate)-x_intersection[0]*x_intersection[0] + (y1_candidate)*(y1_candidate)-y_intersection[0]*y_intersection[0]  ) 
                <     std::abs((x2_candidate)*(x2_candidate)-x_intersection[x_intersection.size() - 1]*x_intersection[x_intersection.size() - 1] + (y2_candidate)*(y2_candidate)-y_intersection[y_intersection.size() - 1]*y_intersection[y_intersection.size() - 1] 
                    + (x2_candidate)*(x2_candidate)-x_intersection[0]*x_intersection[0] + (y2_candidate)*(y2_candidate)-y_intersection[0]*y_intersection[0])    ) {
                x_advanced_surrogate_boundary.push_back(x1_candidate) ;
                y_advanced_surrogate_boundary.push_back(y1_candidate) ;
            }
            else {
                x_advanced_surrogate_boundary.push_back(x2_candidate) ;
                y_advanced_surrogate_boundary.push_back(y2_candidate) ;
            }
        }
        // *Li scrivo al contrario* -> Convention IgaApplication
        // Print on an external file the new control points coordinate for creating the b_rep_trimming curves
        std::ofstream outputFileCoordinate("txt_files/Snake_coordinates2.txt");
        outputFileCoordinate << std::setprecision(16);
        for (int i = x_advanced_surrogate_boundary.size()-1  ; i > 1 ; i-- ) { // *Li scrivo al contrario*
            // Check is there are singularities along the surrogate boundary 
            if (x_advanced_surrogate_boundary[i]!= x_advanced_surrogate_boundary[i+2] || y_advanced_surrogate_boundary[i]!=y_advanced_surrogate_boundary[i+2]) {
                if (x_advanced_surrogate_boundary[i-1]!= x_advanced_surrogate_boundary[i+1] || y_advanced_surrogate_boundary[i-1]!=y_advanced_surrogate_boundary[i+1]) {
                    outputFileCoordinate << x_advanced_surrogate_boundary[i] << " " << y_advanced_surrogate_boundary[i] << std::endl; 
                }
            }
        }
        // *Li scrivo al contrario* -> Convention IgaApplication
        if (x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] != x_advanced_surrogate_boundary[1] || 
                y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1] != y_advanced_surrogate_boundary[1]) {
            outputFileCoordinate << x_advanced_surrogate_boundary[1] << " " << y_advanced_surrogate_boundary[1] << std::endl;
        }

        // // Print on an external file the new control points coordinate for creating the b_rep_trimming curves
        // std::ofstream outputFileCoordinate("txt_files/Snake_coordinates2.txt");
        // outputFileCoordinate << std::setprecision(16);
        // for (int i = 1  ; i < x_advanced_surrogate_boundary.size()-1 ; i++ ) { // *Li scrivo al contrario*
        //     // Check is there are singularities along the surrogate boundary 
        //     if (x_advanced_surrogate_boundary[i]!= x_advanced_surrogate_boundary[i+2] || y_advanced_surrogate_boundary[i]!=y_advanced_surrogate_boundary[i+2]) {
        //         if (x_advanced_surrogate_boundary[i-1]!= x_advanced_surrogate_boundary[i+1] || y_advanced_surrogate_boundary[i-1]!=y_advanced_surrogate_boundary[i+1]) {
        //             outputFileCoordinate << x_advanced_surrogate_boundary[i] << " " << y_advanced_surrogate_boundary[i] << std::endl; 
        //         }
        //     }
        // }
        // // *Li scrivo al contrario* -> Convention IgaApplication
        // if (x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] != x_advanced_surrogate_boundary[1] || 
        //         y_advanced_surrogate_boundary[y_advanced_surrogate_boundary.size()-1] != y_advanced_surrogate_boundary[1]) {
        //     outputFileCoordinate << x_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] << " " << y_advanced_surrogate_boundary[x_advanced_surrogate_boundary.size()-1] << std::endl;
        // }





        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        for (int i = 0; i < (insert_nb_per_span_u+1); i++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                if (knot_spans_available[j][i] == 1 ) {
                    start_i = i;
                    start_j = j;
                    break;
                }
            }
            // KRATOS_WATCH(knot_spans_available[start_j][start_i])
            if (knot_spans_available[start_j][start_i] == 1 ) { break; }
        }
        // KRATOS_WATCH(start_i)
        // KRATOS_WATCH(start_j)
        // KRATOS_WATCH(knot_spans_available)
        // exit(0);

        std::ofstream outputFile("txt_files/Snake_coordinates.txt");
        outputFile << std::setprecision(16);
        outputFile << knots_u[start_i] << " " << knots_v[start_j] << std::endl;
        // KRATOS_WATCH(knots_u[start_i])
        // KRATOS_WATCH(knots_v[start_j])
        // KRATOS_WATCH(start_i)
        // KRATOS_WATCH(start_j)
        // exit(0);
        
        // Follow the clockwise loop
        int end = 0;
        // We are going orizontally
        int direction = 0 ;      // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int i = start_i;         // Indici del percorso lungo i knot values
        int j = start_j;         // Indici del percorso lungo i knot values
        int I = start_i;         // Indici del percorso lungo i knot spans
        int J = start_j;         // Indici del percorso lungo i knot spans
        int steps = 0;
        while (end == 0 && steps < 1000) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // KRATOS_WATCH(i)
            // KRATOS_WATCH(j)
            // KRATOS_WATCH(I)
            // KRATOS_WATCH(J)
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knot_spans_available[J][I] == 1) {
                // KRATOS_WATCH('trovato, sinistra') 
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; risetto e muovo
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}
                // KRATOS_WATCH(I)
                // KRATOS_WATCH(J)
                if (knot_spans_available[J][I] == 1) {
                    // KRATOS_WATCH('trovato, dritto')

                    // Stiamo andando a Dritti! -> Non scrivo nulla e muovo (i,j)
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    outputFile << knots_u[i] << " " << knots_v[j] << std::endl;
                    // KRATOS_WATCH(knots_u[i])
                    // KRATOS_WATCH(knots_v[j])
                }
                else {
                    // Dobbiamo cercare a destra; resetto e muovo
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}

                    if (knot_spans_available[J][I] == 1) {
                        // KRATOS_WATCH('trovato, destra')

                        // Stiamo andando a DX! -> Prima passo dritto, poi stampo, poi passo a destra (i,j), poi scrivo di nuovo
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        outputFile << knots_u[i] << " " << knots_v[j] << std::endl;
                        // KRATOS_WATCH(knots_u[i])
                        // KRATOS_WATCH(knots_v[j])
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        outputFile << knots_u[i] << " " << knots_v[j] << std::endl;
                        // KRATOS_WATCH(knots_u[i])
                        // KRATOS_WATCH(knots_v[j])
                        direction = direction + 1;
                        if (direction == 4) {direction = 0;}
                    }
                    else { // Super special case of "isolated" knot span to be circumnavigated
                        // KRATOS_WATCH('Super special case of "isolated" knot span')
                        is_special_case = 0;
                        // Resetto e muovo (I,J) "indietro"
                        if (direction == 0) {I-- ; J-- ;} 
                        if (direction == 1) {J++ ; I-- ;}
                        if (direction == 2) {I++ ; J++ ;}
                        if (direction == 3) {J-- ; I++ ;}
                        // Check
                        if (knot_spans_available[J][I] == 1) {
                            // First passo dritto, poi print, then move to right, then print again, then move to the right and print again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }
                            outputFile << knots_u[i] << " " << knots_v[j] << std::endl;
                            // KRATOS_WATCH(knots_u[i])
                            // KRATOS_WATCH(knots_v[j])
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }
                            outputFile << knots_u[i] << " " << knots_v[j] << std::endl;
                            // KRATOS_WATCH(knots_u[i])
                            // KRATOS_WATCH(knots_v[j])
                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }
                            outputFile << knots_u[i] << " " << knots_v[j] << std::endl;
                            // KRATOS_WATCH(knots_u[i])
                            // KRATOS_WATCH(knots_v[j])
                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else{ KRATOS_WATCH('errore nello Snakes Coordinates'); exit(0);}  
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i && J == start_j && is_special_case == 1 ) {
                // End of the while loop
                end = 1;
                KRATOS_WATCH(steps)
                }
                

        }
        // outputFile << knots_u[start_i] << " " << knots_v[start_j] << std::endl;
        outputFile.close();
        KRATOS_WATCH('Snake process has finished')
        // exit(0);
    }
}
