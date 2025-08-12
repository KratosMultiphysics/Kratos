//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti

// Project includes
#include "apply_strong_BCS_extended_gradient_method_process.h"

namespace Kratos
{

    ApplyStrongBCSExtendedGradientMethodProcess::ApplyStrongBCSExtendedGradientMethodProcess(
        Model& rModel, Parameters ThisParameters) : mpModel(&rModel), mParameters(ThisParameters)
    {
        ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::Execute(){
        // Read the background mesh model part name
        std::string background_mesh_model_part_name = mParameters["background_mesh_model_part_name"].GetString();
        mpBackgroundMeshModelPart = &mpModel->GetModelPart(background_mesh_model_part_name);
        // Read the skin model part name
        std::string skin_model_part_name = mParameters["skin_model_part_name"].GetString();
        mpSkinModelPart = &mpModel->GetModelPart(skin_model_part_name);
        // Read the unknown variable for the problem
        mUnknownVariableName = mParameters["variable_name"].GetString();
        // Retrieve the Kratos variable from its string name
        mpUnknownVariable = &Kratos::KratosComponents<Kratos::Variable<double>>::Get(mUnknownVariableName);
        // Read the value to be applied as BC
        std::string value = mParameters["value"].GetString();
        // Read the interpolation scheme 
        mInterpolationScheme = mParameters["interpolation_scheme"].GetString();
        // Create the function to be evaluated
        mpEvalFunction = Kratos::make_unique<GenericFunctionUtility>(value);
        // Define the polinomial order of the MLS interpolation
        mMLSPolinomialOrder = mParameters["MLS_polinomial_order"].GetInt();
        // Get the number of nodes of the background mesh 
        mNumberOfNodesBackgroundMesh = mpBackgroundMeshModelPart->GetRootModelPart().NumberOfNodes();
        // Get the number of interpolation points
        mInterpolationPointsNumber = mParameters["number_of_interpolation_points"].GetInt();
        // Get the initial gradient inside the trimmed elements
        mInitialGradientInsideTrimmedElements = mParameters["initial_gradient_inside_trimmed_elements"].GetDouble();

        std::cout << "Applying the Dirichlet BC for: " + mUnknownVariableName << std::endl;

        // Build the LHS in the first iteration
        if (mIterations == 0){
            // Initialize the mLHS and compute it
            mLHS = ZeroMatrix(mNumberOfNodesBackgroundMesh, mNumberOfNodesBackgroundMesh);
            ComputeLHS(mLHS);
            // Build the RHS Boundary Contribution
            mRHSBoundaryContribution = ZeroVector(mNumberOfNodesBackgroundMesh);
            ComputeRHSBoundaryContribution(mRHSBoundaryContribution);  
        }

        // Build the RHS trimmed elements contribution in each iteration 
        mRHSTrimmedElementsContribution = ZeroVector(mNumberOfNodesBackgroundMesh);
        ComputeRHSTrimmedElementsGradientContribution(mRHSTrimmedElementsContribution);

        // Build the RHS
        Vector mRHS = mRHSBoundaryContribution + mRHSTrimmedElementsContribution;

        // Solve the linear system
        mPhiDir = ZeroVector(mNumberOfNodesBackgroundMesh);
        LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
        // We solve the reference system
        psolver->Solve(mLHS, mPhiDir, mRHS);

        // Apply the BCs
        ApplyStrongBoundaryConditions();   

        mIterations += 1;
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::ComputeLHS(CompressedMatrix& LHS){
        // Initialize the two matrices defining the LHS
        CompressedMatrix LHSBoundaryContribution = ZeroMatrix(mNumberOfNodesBackgroundMesh, mNumberOfNodesBackgroundMesh);
        CompressedMatrix LHSTrimmedElementsGradientContribution = ZeroMatrix(mNumberOfNodesBackgroundMesh, mNumberOfNodesBackgroundMesh);

        // Compute the contributions to the LHS
        ComputeLHSBoundaryContribution(LHSBoundaryContribution);
        ComputeLHSTrimmedElementsGradientContribution(LHSTrimmedElementsGradientContribution);

        LHS = LHSBoundaryContribution + LHSTrimmedElementsGradientContribution;
        VerifyAndModifyDiagonalLHS(LHS);
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::ComputeLHSBoundaryContribution(CompressedMatrix& LHSBoundaryContribution){
        IndexType number_active_control_points_per_gp = mpSkinModelPart->ConditionsBegin()->GetGeometry().PointsNumber();

        // Loop over the conditions in intersected_elements_sub_model_part. Each condition represents a gauss point over a brep curve on surface
        for (auto gauss_point_it = mpSkinModelPart->ConditionsBegin(); gauss_point_it != mpSkinModelPart->ConditionsEnd(); gauss_point_it++){
            const auto p_gauss_point_geometry = gauss_point_it->pGetGeometry();
            
            // Initialize the elemental matrix
            CompressedMatrix LHSBoundaryContribution_elemental = ZeroMatrix(number_active_control_points_per_gp, number_active_control_points_per_gp);

            // Get the integration weight  
            const double weight = p_gauss_point_geometry->IntegrationPoints()[0].Weight();  
            
            // Shape functions evaluated at the GP
            const CompressedMatrix& N = p_gauss_point_geometry->ShapeFunctionsValues();

            // Determinant of jacobian
            Vector det_jacobian;
            p_gauss_point_geometry->DeterminantOfJacobian(det_jacobian);
            
            LHSBoundaryContribution_elemental = prod(trans(N), N) * weight * det_jacobian[0];


            // Assemble the elemental contribution in the global matrix
            std::vector<IndexType> equation_id_vector;
            for (IndexType i = 0; i < number_active_control_points_per_gp; i++){
                equation_id_vector.push_back((gauss_point_it->GetGeometry()[i].Id()));
            }

            for (IndexType j = 0; j < equation_id_vector.size(); j++){
                for (IndexType k = 0; k < equation_id_vector.size(); k++){
                    LHSBoundaryContribution(equation_id_vector[j], equation_id_vector[k]) += LHSBoundaryContribution_elemental(j, k);
                }
            }
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::ComputeLHSTrimmedElementsGradientContribution(CompressedMatrix& LHSTrimmedElementsGradientContribution){
        if (mpSkinModelPart->NumberOfElements() == 0) return;

        IndexType number_active_control_points_per_gp = mpSkinModelPart->ElementsBegin()->GetGeometry().PointsNumber();

        // Loop over the elements in intersected_elements_sub_model_part. Each element represents a gauss point over a NURBS surface
        for (auto gauss_point_it = mpSkinModelPart->ElementsBegin(); gauss_point_it != mpSkinModelPart->ElementsEnd(); gauss_point_it++){
            const auto p_gauss_point_geometry = gauss_point_it->pGetGeometry();
            
            // Initialize the elemental matrix
            Matrix LHSTrimmedElementsContribution_elemental = ZeroMatrix(number_active_control_points_per_gp, number_active_control_points_per_gp);

            // Get the integration weight  
            const double weight = p_gauss_point_geometry->IntegrationPoints()[0].Weight();  
            
            // Shape functions evaluated at the GP
            const GeometryType::ShapeFunctionsGradientsType& dN_dxi = p_gauss_point_geometry->ShapeFunctionsLocalGradients();
            const unsigned int dim = dN_dxi[0].size2();
            Matrix dN_dx(number_active_control_points_per_gp, dim); 
            noalias(dN_dx) = dN_dxi[0];

            // Determinant of jacobian
            Vector det_jacobian;
            p_gauss_point_geometry->DeterminantOfJacobian(det_jacobian);

            double beta = 1.0;
            bool use_full_gradient = false;
            if (use_full_gradient){
                LHSTrimmedElementsContribution_elemental = beta * prod(dN_dx, trans(dN_dx)) * weight * det_jacobian[0];
            } else {
                // Calculate the normal vector of the boundary inside the element
                Vector unit_normal_vector = gauss_point_it -> GetValue(UNIT_NORMAL);
                Vector unit_normal_vector_2d = ZeroVector(2); 
                for (IndexType i = 0; i < 2; ++i) unit_normal_vector_2d[i] = unit_normal_vector[i];

                Matrix normal_prod =  outer_prod(unit_normal_vector_2d, unit_normal_vector_2d);
                Matrix first_prod = prod(dN_dx, normal_prod);
                LHSTrimmedElementsContribution_elemental = beta * prod(first_prod, trans(dN_dx)) * weight * det_jacobian[0];
            }

            // Assemble the elemental contribution in the global matrix
            std::vector<IndexType> equation_id_vector;
            for (IndexType i = 0; i < number_active_control_points_per_gp; i++){
                equation_id_vector.push_back(gauss_point_it->GetGeometry()[i].Id());
            }
            
            for (IndexType j = 0; j < equation_id_vector.size(); j++){
                for (IndexType k = 0; k < equation_id_vector.size(); k++){
                    LHSTrimmedElementsGradientContribution(equation_id_vector[j], equation_id_vector[k]) += LHSTrimmedElementsContribution_elemental(j, k);
                }
            }
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::VerifyAndModifyDiagonalLHS(CompressedMatrix& LHS) {
        for (IndexType i = 0; i < LHS.size1(); ++i) { // Iterate over rows
            // Check if the entire row is zero
            bool is_zero_row = true;
            for (std::size_t j = 0; j < LHS.size2(); ++j) {
                if (LHS(i, j) != 0.0) {
                    is_zero_row = false;
                    break;
                }
            }

            // If the row is full of zeros, modify the diagonal element
            if (is_zero_row && i < LHS.size2()) { // Ensure diagonal exists
                LHS(i, i) = 1.0;
            }
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::ComputeRHSBoundaryContribution(Vector& RHSBoundaryContribution){
        IndexType number_active_control_points_per_gp = mpSkinModelPart->ConditionsBegin()->GetGeometry().PointsNumber();

        // Loop over the conditions in intersected_elements_sub_model_part. Each condition represents a gauss point over a brep curve on surface
        for (auto gauss_point_it = mpSkinModelPart->ConditionsBegin(); gauss_point_it != mpSkinModelPart->ConditionsEnd(); gauss_point_it++){
            const auto p_gauss_point_geometry = gauss_point_it->pGetGeometry();
            const auto gauss_point_position = p_gauss_point_geometry->Center();
            
            // Initialize the elemental vector
            Vector RHSBoundaryContribution_elemental = ZeroVector(number_active_control_points_per_gp);

            // Get the integration weight  
            const double weight = p_gauss_point_geometry->IntegrationPoints()[0].Weight();  
            
            // Shape functions evaluated at the GP
            const CompressedMatrix& N = p_gauss_point_geometry->ShapeFunctionsValues();
            auto N_row = row(N, 0);

            // Determinant of jacobian
            Vector det_jacobian;
            p_gauss_point_geometry->DeterminantOfJacobian(det_jacobian);

            // Evaluation of the boundary condition at the gauss point
            double current_time = mpSkinModelPart->GetProcessInfo()[TIME];
            const double value = mpEvalFunction->CallFunction(gauss_point_position.X(), gauss_point_position.Y(), gauss_point_position.Z(), current_time);

            double beta = 1.0;
            RHSBoundaryContribution_elemental = beta * N_row * value * weight * det_jacobian[0];

            // Assemble the elemental contribution in the global matrix
            std::vector<IndexType> equation_id_vector;
            for (IndexType i = 0; i < number_active_control_points_per_gp; i++){
                equation_id_vector.push_back(gauss_point_it->GetGeometry()[i].Id());
            }
            
            for (IndexType j = 0; j < equation_id_vector.size(); j++){
                RHSBoundaryContribution(equation_id_vector[j]) += RHSBoundaryContribution_elemental(j);
            }
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::ComputeRHSTrimmedElementsGradientContribution(Vector& RHSTrimmedElementsGradientContribution){
        if (mpSkinModelPart->NumberOfElements() == 0) return;

        IndexType number_active_control_points_per_gp = mpSkinModelPart->ElementsBegin()->GetGeometry().PointsNumber();

        // Define the interpolation points and the fields for the rbf interpolation 
        DefineInterpolationPointsAndGradients(mInterpolationPointsCoordinates, mSolution, mSolutionGradientX, mSolutionGradientY, mInterpolationPointsPointers);

        // Loop over the elements in intersected_elements_sub_model_part. Each element represents a gauss point over a NURBS surface
        for (auto gauss_point_it = mpSkinModelPart->ElementsBegin(); gauss_point_it != mpSkinModelPart->ElementsEnd(); gauss_point_it++){
            const auto p_gauss_point_geometry = gauss_point_it->pGetGeometry();

            // Gauss point position
            double gauss_point_x = p_gauss_point_geometry->Center()[0];
            double gauss_point_y = p_gauss_point_geometry->Center()[1];

            // Initialize the elemental vector
            Vector RHSTrimmedElementsContribution_elemental = ZeroVector(number_active_control_points_per_gp);

            // Get the integration weight  
            const double weight = p_gauss_point_geometry->IntegrationPoints()[0].Weight();  
            
            // Shape functions evaluated at the GP
            const GeometryType::ShapeFunctionsGradientsType& dN_dxi = p_gauss_point_geometry->ShapeFunctionsLocalGradients();
            const unsigned int dim = dN_dxi[0].size2();
            Matrix dN_dx(number_active_control_points_per_gp, dim); 
            noalias(dN_dx) = dN_dxi[0];

            // Determinant of jacobian
            Vector det_jacobian;
            p_gauss_point_geometry->DeterminantOfJacobian(det_jacobian);

            // Calculate the interpolated solution gradient at the gp position
            Vector interpolated_solution_gradient = ZeroVector(3);
            Vector exact_solution_gradient = ZeroVector(3);
            double interpolated_solution = 0.0;
            
            if (mIterations == 0){
                interpolated_solution_gradient[0] = mInitialGradientInsideTrimmedElements;
                interpolated_solution_gradient[1] = mInitialGradientInsideTrimmedElements;
            } else {
                InterpolateSolutionGradients(interpolated_solution_gradient, interpolated_solution, p_gauss_point_geometry->Center(), mInterpolationPointsNumber);
            }
            
            exact_solution_gradient[0]=(M_PI*std::cos(M_PI * p_gauss_point_geometry->Center()[0])*std::cos(M_PI*p_gauss_point_geometry->Center()[1]));
            exact_solution_gradient[1] = (-M_PI*std::sin(M_PI * p_gauss_point_geometry->Center()[0])*std::sin(M_PI*p_gauss_point_geometry->Center()[1]));
            // if(mUnknownVariableName == "DISPLACEMENT_X" or mUnknownVariableName == "DISPLACEMENT_Y"){
            //     exact_solution_gradient[0]=0.0;
            //     exact_solution_gradient[1]=0.0;
            // } else if(mUnknownVariableName == "DISPLACEMENT_Z"){
            //     ComputeExactGradient(gauss_point_x, gauss_point_y, exact_solution_gradient);
            // }

            
            // std::cout << "----------------" << std::endl;
            // KRATOS_WATCH(p_gauss_point_geometry->Center())
            // KRATOS_WATCH(exact_solution_gradient[0])
            // KRATOS_WATCH(exact_solution_gradient[1])
            // KRATOS_WATCH(interpolated_solution_gradient[0])
            // KRATOS_WATCH(interpolated_solution_gradient[1])


            double beta = 1.0;
            bool use_full_gradient = false;
            if (use_full_gradient){
                // Full gradient regularization
                Vector gradient_dot_dN = prod(dN_dx, interpolated_solution_gradient);
                RHSTrimmedElementsContribution_elemental = beta * gradient_dot_dN * weight * det_jacobian[0];
            } else{
                // Calculate the normal vector of the boundary inside the element
                Vector unit_normal_vector = gauss_point_it -> GetValue(UNIT_NORMAL);
                Vector unit_normal_vector_2d = ZeroVector(2); 
                for (IndexType i = 0; i < 2; ++i) unit_normal_vector_2d[i] = unit_normal_vector[i];

                // Normal component regularization (original behavior)
                double gp_normal_gradient = inner_prod(interpolated_solution_gradient, unit_normal_vector);
                RHSTrimmedElementsContribution_elemental = beta * prod(dN_dx, unit_normal_vector_2d) * gp_normal_gradient * weight * det_jacobian[0];
            }

            // Assemble the elemental contribution in the global matrix
            std::vector<IndexType> equation_id_vector;
            for (IndexType i = 0; i < number_active_control_points_per_gp; i++){
                equation_id_vector.push_back(gauss_point_it->GetGeometry()[i].Id());
            }
            
            for (IndexType j = 0; j < equation_id_vector.size(); j++){
                RHSTrimmedElementsGradientContribution(equation_id_vector[j]) += RHSTrimmedElementsContribution_elemental(j);
            }
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::InterpolateSolutionGradients(Vector& SolutionGradients, double& Solution, CoordinatesArrayType quadrature_point_position, IndexType NumberOfClosestPoints){
        KRATOS_ERROR_IF(mInterpolationPointsCoordinates.size1() != mSolutionGradientX.size())
        << "The number of points must match the size of solution gradients." << std::endl;

        FindTheNClosestInterpolationPoints(NumberOfClosestPoints, quadrature_point_position);

        double h_MLS = 0.0;

        if (mIterations == 0){
            h_MLS = 0.1;
        }
        else if (mIterations != 0){
            h_MLS = CalculateKernelParameterMLS(quadrature_point_position);
        }
        
        if (mInterpolationScheme == "RBF"){
            // Using RBF
            // Compute the shape parameter (h)
            double h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(mClosestInterpolationPointsCoordinates);

            // Compute RBF shape functions
            Vector rN_RBF;
            RBFShapeFunctionsUtility::CalculateShapeFunctions(mClosestInterpolationPointsCoordinates, quadrature_point_position, h, rN_RBF, nullptr);

            // Interpolate the field value
            double interpolated_solution = RBFShapeFunctionsUtility::CalculateShapeFunctionsAndInterpolation(
                mClosestInterpolationPointsCoordinates, quadrature_point_position, h, rN_RBF, mSolution
            );
            double interpolated_solution_x_gradient = RBFShapeFunctionsUtility::CalculateShapeFunctionsAndInterpolation(
                mClosestInterpolationPointsCoordinates, quadrature_point_position, h, rN_RBF, mClosestPointsSolutionGradientX
            );
            double interpolated_solution_y_gradient = RBFShapeFunctionsUtility::CalculateShapeFunctionsAndInterpolation(
                mClosestInterpolationPointsCoordinates, quadrature_point_position, h, rN_RBF, mClosestPointsSolutionGradientY
            );

            SolutionGradients[0] = interpolated_solution_x_gradient;
            SolutionGradients[1] = interpolated_solution_y_gradient;
        }
        else if (mInterpolationScheme == "MLS")
        {
            // Using MLS (Moving Least Squares)
            Vector rN_MLS;
            if (mMLSPolinomialOrder == 1){
                MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(mClosestInterpolationPointsCoordinates, quadrature_point_position, h_MLS, rN_MLS);
            }
            else if (mMLSPolinomialOrder == 2){
                MLSShapeFunctionsUtility::CalculateShapeFunctions<3,2>(mClosestInterpolationPointsCoordinates, quadrature_point_position, h_MLS, rN_MLS);
            }

            // Interpolate the field
            double interpolatedValue_gradx = 0.0;
            double interpolatedValue_grady = 0.0;
            double interpolatedValue_sol = 0.0;

            for (IndexType i = 0; i < mClosestPointsSolutionGradientX.size(); ++i) {
                interpolatedValue_gradx += rN_MLS[i] * mClosestPointsSolutionGradientX[i];
                interpolatedValue_grady += rN_MLS[i] * mClosestPointsSolutionGradientY[i];
                interpolatedValue_sol += rN_MLS[i] * mClosestPointsSolution[i];
            }

            SolutionGradients[0] = interpolatedValue_gradx;
            SolutionGradients[1] = interpolatedValue_grady;
        }
        else
        {
            std::cout << "Please specify an interpolation scheme for the gradient" << std::endl;
        }
        
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::DefineInterpolationPointsAndGradients(Matrix& InterpolationPoints, Vector& Solution, Vector& SolutionGradientX, Vector& SolutionGradientY, std::vector<GeometryPointerType>& InterpolationPointsPointers){
        // Define the sub model part used for the interpolation
        ModelPart& interpolation_sub_model_part = mpModel->GetModelPart("IgaBackgroundMeshModelPart.active_elements");
        IndexType number_active_control_points_per_gp = interpolation_sub_model_part.ElementsBegin()->GetGeometry().PointsNumber();

        // First, define the number of interpolation points to resize the matrix and vectors
        IndexType interpolation_points_counter = 0;
        for (auto gauss_point_it = interpolation_sub_model_part.ElementsBegin(); gauss_point_it != interpolation_sub_model_part.ElementsEnd(); gauss_point_it++){
            const auto p_gauss_point_geometry = gauss_point_it->pGetGeometry();
            const auto gp_center = p_gauss_point_geometry->Center();

            // For the interpolation, we do not use the gauss points which are already fixed by the algorithm 
            bool skip_element = false;
            for (IndexType i = 0; i < number_active_control_points_per_gp; i++){
                if (gauss_point_it->GetGeometry()[i].GetDof(*mpUnknownVariable).IsFixed() == 1){
                    skip_element = true;
                    break;
                }
            }

            if (skip_element == true) continue;
            interpolation_points_counter += 1;
        }
       
       // Resize the matrices and vectors for: interpolation points, solution, gradient in x and y directions
       InterpolationPoints.clear();
       InterpolationPoints.resize(interpolation_points_counter, 3);
       Solution.clear();
       Solution.resize(interpolation_points_counter);
       SolutionGradientX.clear();
       SolutionGradientX.resize(interpolation_points_counter);
       SolutionGradientY.clear();
       SolutionGradientY.resize(interpolation_points_counter);

       // Fill the matrix and vectors
       IndexType i = 0;
       for (auto gauss_point_it = interpolation_sub_model_part.ElementsBegin(); gauss_point_it != interpolation_sub_model_part.ElementsEnd(); gauss_point_it++){
            const auto p_gauss_point_geometry = gauss_point_it->pGetGeometry();
            const auto gp_center = p_gauss_point_geometry->Center();

            // For the interpolation, we do not use the gauss points which are already fixed by the algorithm 
            bool skip_element = false;
            Vector nodal_solution = ZeroVector(number_active_control_points_per_gp);
            for (IndexType k = 0; k < number_active_control_points_per_gp; k++){
                if (gauss_point_it->GetGeometry()[k].GetDof(*mpUnknownVariable).IsFixed() == 1){
                    skip_element = true;
                    break;
                }
                nodal_solution[k] = gauss_point_it->GetGeometry()[k].GetSolutionStepValue(*mpUnknownVariable);
            }

            if (skip_element == true) continue;

            const GeometryType::ShapeFunctionsGradientsType& dN_dxi = p_gauss_point_geometry->ShapeFunctionsLocalGradients();
            const unsigned int dim = dN_dxi[0].size2();
            CompressedMatrix dN_dx(number_active_control_points_per_gp, dim); 
            noalias(dN_dx) = dN_dxi[0];

            Vector solution_gradients = prod(trans(dN_dx), nodal_solution);

            if (mUnknownVariableName == "DISPLACEMENT_Z"){
                Vector exact_solution_gradient = ZeroVector(3);
                ComputeExactGradient(gp_center.X(), gp_center.Y(), exact_solution_gradient);
                // KRATOS_WATCH(gauss_point_it->Id())
                // KRATOS_WATCH(gp_center)
                // KRATOS_WATCH(solution_gradients)
                // KRATOS_WATCH(exact_solution_gradient)
                SolutionGradientX[i] = solution_gradients[0];
                SolutionGradientY[i] = solution_gradients[1];
            } else{
                // Fill the vectors and matrices with the required information 
                // KRATOS_WATCH(gauss_point_it->Id())
                // KRATOS_WATCH(gp_center)
                SolutionGradientX[i] = solution_gradients[0];
                SolutionGradientY[i] = solution_gradients[1];
            }

            InterpolationPoints(i, 0) = gp_center.X();
            InterpolationPoints(i, 1) = gp_center.Y();
            InterpolationPoints(i, 2) = gp_center.Z();

            InterpolationPointsPointers.push_back(p_gauss_point_geometry);

            i += 1;
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::ApplyStrongBoundaryConditions(){
        IndexType number_active_control_points_per_gp = mpSkinModelPart->ConditionsBegin()->GetGeometry().PointsNumber();

        for (auto gauss_point_it = mpSkinModelPart->ConditionsBegin(); gauss_point_it != mpSkinModelPart->ConditionsEnd(); gauss_point_it++){
            for (IndexType i = 0; i < number_active_control_points_per_gp; i++){
                IndexType node_id = gauss_point_it->GetGeometry()[i].Id();
                double value_to_fix = mPhiDir[node_id];
                if (mLHS(node_id,node_id) == 1.0) continue;
                gauss_point_it->GetGeometry()[i].Fix(*mpUnknownVariable);
                gauss_point_it->GetGeometry()[i].FastGetSolutionStepValue(*mpUnknownVariable) = value_to_fix;
            }
        }
    }

    void ApplyStrongBCSExtendedGradientMethodProcess::FindTheNClosestInterpolationPoints(IndexType NumberOfClosestPoints, CoordinatesArrayType gauss_point_position){
        // The closest points could be, as maximum, the number of interpolation points
        if (NumberOfClosestPoints > mInterpolationPointsCoordinates.size1()){
            NumberOfClosestPoints = mInterpolationPointsCoordinates.size1();
        }

        // Clear and resize the matrix and vector 
        mClosestInterpolationPointsCoordinates.clear(); 
        mClosestInterpolationPointsCoordinates.resize(NumberOfClosestPoints, 3); 

        mClosestPointsSolution.clear(); 
        mClosestPointsSolution.resize(NumberOfClosestPoints);

        mClosestPointsSolutionGradientX.clear(); 
        mClosestPointsSolutionGradientX.resize(NumberOfClosestPoints);

        mClosestPointsSolutionGradientY.clear(); 
        mClosestPointsSolutionGradientY.resize(NumberOfClosestPoints); 

        mClosestInterpolationPointsPointers.clear(); 

        std::vector<std::pair<double, IndexType>> distance_to_point_and_index ;
        distance_to_point_and_index.clear();

        // Fill the distance to point and index 
        for (IndexType i = 0; i < mInterpolationPointsCoordinates.size1(); i++){
            double gp_x_coordinate = mInterpolationPointsCoordinates(i, 0);
            double gp_y_coordinate = mInterpolationPointsCoordinates(i, 1);
            double gp_z_coordinate = mInterpolationPointsCoordinates(i, 2);
            double distance = std::sqrt(std::pow(gp_x_coordinate - gauss_point_position[0], 2) +
                     std::pow(gp_y_coordinate - gauss_point_position[1], 2) +
                     std::pow(gp_z_coordinate - gauss_point_position[2], 2));
            distance_to_point_and_index.emplace_back(distance, i);
        }

        // Sort the vector considering the distance (from smallest to biggest)
        std::sort(distance_to_point_and_index.begin(), distance_to_point_and_index.end(),
              [](const std::pair<double, IndexType>& a, const std::pair<double, IndexType>& b) {
                  return a.first < b.first; 
              });


        // Fill the closest points matrices and vectors
        for (IndexType i = 0; i < NumberOfClosestPoints; i++) {
            IndexType index = distance_to_point_and_index[i].second;

            // Fill the closest points coordinates matrix
            mClosestInterpolationPointsCoordinates(i, 0) = mInterpolationPointsCoordinates(index, 0);
            mClosestInterpolationPointsCoordinates(i, 1) = mInterpolationPointsCoordinates(index, 1);
            mClosestInterpolationPointsCoordinates(i, 2) = mInterpolationPointsCoordinates(index, 2);

            // Fill the solution gradient vectors
            mClosestPointsSolution[i] = mSolution[index];
            mClosestPointsSolutionGradientX[i] = mSolutionGradientX[index];
            mClosestPointsSolutionGradientY[i] = mSolutionGradientY[index];

            mClosestInterpolationPointsPointers.push_back(mInterpolationPointsPointers[index]);
        }   
    }

     double ApplyStrongBCSExtendedGradientMethodProcess::CalculateKernelParameterMLS(CoordinatesArrayType quadrature_point_position){
        double min_distance = 1.0e10;

        for (IndexType i = 0; i < mClosestInterpolationPointsCoordinates.size1(); i++){
            array_1d<double, 3> interpolation_points_coordinates;

            interpolation_points_coordinates[0] = mClosestInterpolationPointsCoordinates(i, 0);
            interpolation_points_coordinates[1] = mClosestInterpolationPointsCoordinates(i, 1);
            interpolation_points_coordinates[2] = mClosestInterpolationPointsCoordinates(i, 2);

            double distance = norm_2(quadrature_point_position - interpolation_points_coordinates);
            if (distance < min_distance){
                min_distance = distance;
            }
        }
        
        return 1.3 * min_distance;
     }

    void ApplyStrongBCSExtendedGradientMethodProcess::ComputeExactGradient(double x, double y, Vector& exact_solution_gradient){
        double sum_x = 0.0;
        double a = 1.0;
        double b = 1.0;
        double q0 = -1.0;
        double E = 3.0e9;
        double D = (E * std::pow(0.01, 3)) / (12.0 * (1.0 - std::pow(0.3, 2)));

        for (int m = 1; m <= 100; m += 2) {  // Only odd values of m
            for (int n = 1; n <= 100; n += 2) {  // Only odd values of n
                double sin_ny = std::sin(n * M_PI * (y-0.5) / b);
                double cos_mx = std::cos(m * M_PI * (x-0.5) / a);
                double term = (m * M_PI / a) * cos_mx * sin_ny;

                double denominator = m * n * std::pow((std::pow(m / a, 2) + std::pow(n / b, 2)), 2);
                sum_x += term / denominator;
            }
        }

        exact_solution_gradient[0] = (16 * q0 / (std::pow(M_PI, 6) * D)) * sum_x;

        double sum_y = 0.0;

        for (int m = 1; m <= 100; m += 2) {  // Only odd values of m
            for (int n = 1; n <= 100; n += 2) {  // Only odd values of n
                double sin_mx = std::sin(m * M_PI * (x-0.5) / a);
                double cos_ny = std::cos(n * M_PI * (y-0.5) / b);
                double term = (n * M_PI / b) * sin_mx * cos_ny;

                double denominator = m * n * std::pow((std::pow(m / a, 2) + std::pow(n / b, 2)), 2);
                sum_y += term / denominator;
            }
        }

        exact_solution_gradient[1] = (16 * q0 / (std::pow(M_PI, 6) * D)) * sum_y;
    }


} // End namespace Kratos
