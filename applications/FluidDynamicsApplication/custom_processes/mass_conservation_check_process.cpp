//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski
//
//

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "mass_conservation_check_process.h"


namespace Kratos
{

/* Public functions *******************************************************/

MassConservationCheckProcess::MassConservationCheckProcess(
        ModelPart& rModelPart,
        const int MassComputationFreq,
        const bool CompareToInitial,
        const bool WriteToLogFile)
    : Process(), mrModelPart(rModelPart) {

    mMassComputationFreq = MassComputationFreq;
    mCompareToInitial = CompareToInitial;
    mWriteToLogFile = WriteToLogFile;
}


MassConservationCheckProcess::MassConservationCheckProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(), mrModelPart(rModelPart) {

    Parameters default_parameters( R"(
    {
        "model_part_name"                        : "default_model_part_name",
        "mass_computation_frequency"             : 20,
        "compare_to_initial_values"              : true,
        "write_to_log_file"                      : true
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMassComputationFreq = rParameters["mass_computation_frequency"].GetInt();
    mCompareToInitial = rParameters["compare_to_initial_values"].GetBool();
    mWriteToLogFile = rParameters["write_to_log_file"].GetBool();
}


bool MassConservationCheckProcess::GetUpdateStatus(){

    const int time_step = mrModelPart.GetProcessInfo()[STEP];
    double interfaceArea;

    if ( mCompareToInitial && time_step == 1){
        // getting initial values and storing as reference
        this->ComputeVolumesOfFluids( mInitialPositiveVolume, mInitialNegativeVolume, interfaceArea );
    }

    if ( time_step % this->mMassComputationFreq == 0){

        // writing an output at a given frequncy
        this->ComputeVolumesOfFluids( mCurrentPositiveVolume, mCurrentNegativeVolume, mCurrentInterfaceArea );
        mIsUpdated = true;
        return true;

    } else {

        mCurrentPositiveVolume = -1.0;
        mCurrentNegativeVolume = -1.0;
        mCurrentInterfaceArea = -1.0;
        return false;
    }
}


/* Private functions ****************************************************/

void MassConservationCheckProcess::ComputeVolumesOfFluids( double& positiveVolume, double& negativeVolume, double& interfaceArea ){

    // useless containers (needed as dummy arguments)
    Matrix rShapeFunctions;
    GeometryType::ShapeFunctionsGradientsType rShapeDerivatives;

    // abbreviation of variable names
    double& posVol = positiveVolume;
    double& negVol = negativeVolume;
    double& interA = interfaceArea;

    // initalisation
    posVol = 0.0;
    negVol = 0.0;
    interA = 0.0;

    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        // iteration over all elements
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;

        auto p_geom = it_elem->pGetGeometry();
        unsigned int ptCountPos = 0;
        unsigned int ptCountNeg = 0;

        // instead of using data.isCut()
        for (unsigned int pt = 0; pt < p_geom->Points().size(); pt++){
            if ( p_geom->GetPoint(pt).FastGetSolutionStepValue(DISTANCE) > 0.0 ){
                ptCountPos++;
            } else {
                ptCountNeg++;
            }
        }

        if ( ptCountPos == p_geom->PointsNumber() ){
            // all nodes are positive
            posVol += p_geom->DomainSize();;
        }
        else if ( ptCountNeg == p_geom->PointsNumber() ){
            // all nodes are negative
            negVol += p_geom->DomainSize();
        }
        else {
            // element is cut by the surface (splitting)
            ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
            Vector w_gauss_pos_side(3, 0.0);
            Vector w_gauss_neg_side(3, 0.0);
            Vector w_gauss_interface(3, 0.0);

            Vector Distance( p_geom->PointsNumber(), 0.0 );
            for (unsigned int i = 0; i < p_geom->PointsNumber(); i++){
                Distance[i] = p_geom->GetPoint(i).FastGetSolutionStepValue(DISTANCE);
            }

            if ( p_geom->PointsNumber() == 3 ){ p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, Distance); }
            else if ( p_geom->PointsNumber() == 4 ){ p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, Distance); }
            else { KRATOS_ERROR << "The process can not be applied on this kind of element" << std::endl; }

            // Call the positive side modified shape functions calculator (Gauss weights woulb be enough)
            // Object p_modified_sh_func has full knowledge of slit geometry
            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                    rShapeFunctions,                    // N
                    rShapeDerivatives,                  // DN
                    w_gauss_pos_side,                   // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration (1 point per triangle)

            for ( unsigned int i = 0; i < w_gauss_pos_side.size(); i++){
                posVol += w_gauss_pos_side[i];
            }

            // Call the negative side modified shape functions calculator
            // Object p_modified_sh_func has full knowledge of slit geometry
            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                    rShapeFunctions,                    // N
                    rShapeDerivatives,                  // DN
                    w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration

            for ( unsigned int i = 0; i < w_gauss_neg_side.size(); i++){
                negVol += w_gauss_neg_side[i];
            }

            p_modified_sh_func->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
                    rShapeFunctions,                    // N
                    rShapeDerivatives,                  // DN
                    w_gauss_interface,                  // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration

            for ( unsigned int i = 0; i < w_gauss_interface.size(); i++){
                interA += std::abs( w_gauss_interface[i] );
            }
        }
    }
}


double MassConservationCheckProcess::ComputeOutletVolumeFlow(){

    // Convention: "mass" is considered as "water", meaning the volumes with a negative distance is considered
    double netInflow = 0.0;
    const double epsilon = 1.0e-7;

    for (int i_cond = 0; i_cond < static_cast<int>(mrModelPart.NumberOfConditions()); ++i_cond){
        // iteration over all elements

        auto condition = mrModelPart.ConditionsBegin() + i_cond;

        if ( condition->Is(OUTLET) ){

            std::cout << "--------------------outlet found --------------------------------------------"

            const int dim = condition->GetGeometry().Dimension();
            const Geometry<Node<3> >& p_geom = condition->GetGeometry();
            Vector Distance( p_geom.PointsNumber(), 0.0 );

            unsigned int negCount = 0;
            unsigned int posCount = 0;
            for (unsigned int i = 0; i < p_geom.PointsNumber(); i++){
                Distance[i] = condition->GetGeometry()[i].GetSolutionStepValue( DISTANCE );
                if ( condition->GetGeometry()[i].GetSolutionStepValue( DISTANCE ) > 0.0 ){
                    posCount++;
                } else {
                    negCount++;
                }
            }

            // leave the current iteration of the condition is completely on positive side
            if ( posCount == p_geom.PointsNumber() ){ continue; }

            if (dim == 2){      // 2D case: condition is a line

                array_1d<double, 3> normal;
                this->CalculateUnitNormal2D( normal, p_geom);
                if( norm_2( normal ) < epsilon ){ continue; }
                else { normal /= norm_2( normal ); }

                if ( negCount == p_geom.PointsNumber() ){       // the condition is completely on the negative side

                    // Gauss point information
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                    const unsigned int NumGauss = IntegrationPoints.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    p_geom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
                    const Matrix Ncontainer = p_geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                    for (unsigned int i_gauss = 0; i_gauss < NumGauss; i_gauss++){
                        const Vector N = row(Ncontainer, i_gauss);
                        double const wGauss = GaussPtsJDet[i_gauss] * IntegrationPoints[i_gauss].Weight();
                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        for (unsigned int n_node = 0; n_node < p_geom.PointsNumber(); n_node++){
                            interpolatedVelocity += N[n_node] * condition->GetGeometry()[n_node].GetSolutionStepValue(VELOCITY);
                        }
                        netInflow += - wGauss * inner_prod( normal, interpolatedVelocity );
                    }

                } else if ( negCount < p_geom.PointsNumber() && negCount < p_geom.PointsNumber() ){      // the condition is cut

                    // Compute the relative coordinate of the intersection point over the edge
                    const double aux_node_rel_location = std::abs ( Distance[0] /( Distance[1]-Distance[0] ));
                    // Shape function values at the position where the surface cuts the element
                    Vector Ncut = ZeroVector(2);
                    Ncut[0] = 1.0 - aux_node_rel_location;
                    Ncut[1] = aux_node_rel_location;

                    // Creation of an auxiliary geometry which covers the negative volume domain only
                    // (imitation of the splitting mechanism for triangles)
                    PointerVectorSet<IndexedPoint, IndexedObject> aux_point_container;
                    aux_point_container.reserve(2);
                    array_1d<double, 3> aux_point1_coords, aux_point2_coords, aux_velocity1, aux_velocity2;

                    IndexedPoint::Pointer paux_point1 = nullptr;
                    IndexedPoint::Pointer paux_point2 = nullptr;
                    // Modify the node with the positive distance
                    for (unsigned int i_node = 0; i_node < p_geom.PointsNumber(); i_node++){
                        if ( p_geom[i_node].GetSolutionStepValue( DISTANCE ) > 0.0 ){
                            // interpolating position for the new node
                            aux_point1_coords[0] = Ncut[0] * p_geom[0].X() + Ncut[1] * p_geom[1].X();
                            aux_point1_coords[1] = Ncut[0] * p_geom[0].Y() + Ncut[1] * p_geom[1].Y();
                            aux_point1_coords[2] = Ncut[0] * p_geom[0].Z() + Ncut[1] * p_geom[1].Z();

                            paux_point1 = Kratos::make_shared<IndexedPoint>(aux_point1_coords, mrModelPart.NumberOfNodes()+1);
                            aux_velocity1 = Ncut[0] * p_geom[0].GetSolutionStepValue( VELOCITY ) + Ncut[1] * p_geom[1].GetSolutionStepValue( VELOCITY );

                        } else {
                            aux_point2_coords[0] = p_geom[i_node].X();
                            aux_point2_coords[1] = p_geom[i_node].Y();
                            aux_point2_coords[2] = p_geom[i_node].Z();

                            paux_point2 = Kratos::make_shared<IndexedPoint>(aux_point2_coords, mrModelPart.NumberOfNodes()+2);
                            aux_velocity2 = p_geom[i_node].GetSolutionStepValue( VELOCITY );
                        }
                    }

                    Geometry <IndexedPoint>::Pointer p_aux_line = Kratos::make_shared< Line3D2 < IndexedPoint > >( paux_point1, paux_point2 );
                    // Gauss point information for auxiliary geometry
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_aux_line->IntegrationPoints( GeometryData::GI_GAUSS_2 );
                    const unsigned int NumGauss = IntegrationPoints.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    p_aux_line->DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
                    const Matrix Ncontainer = p_aux_line->ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                    for (unsigned int i_gauss = 0; i_gauss < NumGauss; i_gauss++){
                        const Vector N = row(Ncontainer, i_gauss);
                        double const wGauss = GaussPtsJDet[i_gauss] * IntegrationPoints[i_gauss].Weight();

                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        interpolatedVelocity = N[0] * aux_velocity1 + N[1] * aux_velocity2;
                        netInflow += - wGauss * inner_prod( normal, interpolatedVelocity );
                    }
                }
            }

            else if (dim == 3){
                // 3D case: condition is a triangle (implemented routines can be used)

                array_1d<double, 3> normal;
                this->CalculateUnitNormal3D( normal, p_geom);
                if( norm_2( normal ) < epsilon ){ continue; }
                else { normal /= norm_2( normal ); }

                if ( negCount == p_geom.PointsNumber() ){       // the condition is completely on the negative side

                    // Gauss point information
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                    const unsigned int NumGauss = IntegrationPoints.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    p_geom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
                    const Matrix Ncontainer = p_geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                    for (unsigned int i_gauss = 0; i_gauss < NumGauss; i_gauss++){
                        const Vector N = row(Ncontainer, i_gauss);
                        double const wGauss = GaussPtsJDet[i_gauss] * IntegrationPoints[i_gauss].Weight();
                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        for (unsigned int n_node = 0; n_node < p_geom.PointsNumber(); n_node++){
                            interpolatedVelocity += N[n_node] * condition->GetGeometry()[n_node].GetSolutionStepValue(VELOCITY);
                        }
                        netInflow += - wGauss * inner_prod( normal, interpolatedVelocity );
                    }

                } else if ( negCount < p_geom.PointsNumber() && negCount < p_geom.PointsNumber() ){      // the condition is cut

                    Matrix rShapeFunctions;
                    GeometryType::ShapeFunctionsGradientsType rShapeDerivatives;
                    Vector w_gauss_neg_side;
                    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
                    auto geomPointer = condition->pGetGeometry();
                    p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(geomPointer, Distance);

                    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                        rShapeFunctions,                    // N
                        rShapeDerivatives,                  // DN
                        w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                        GeometryData::GI_GAUSS_2);          // first order Gauss integration

                    // interating velocity over the negative area of the condition
                    for ( unsigned int i_gauss = 0; i_gauss < w_gauss_neg_side.size(); i_gauss++){
                        const array_1d<double,3>& N = row(rShapeFunctions, i_gauss);

                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        for (unsigned int n_node = 0; n_node < p_geom.PointsNumber(); n_node++){
                            interpolatedVelocity += N[n_node] * condition->GetGeometry()[n_node].GetSolutionStepValue(VELOCITY);
                        }
                        netInflow += - w_gauss_neg_side[i_gauss] * inner_prod( normal, interpolatedVelocity );
                    }
                }
            }
        }
    }
    return netInflow;
}


double MassConservationCheckProcess::ComputeInletVolumeFlow(){

    // Convention: "mass" is considered as "water", meaning the volumes with a negative distance is considered
    double netInflow = 0.0;
    const double epsilon = 1.0e-7;

    for (int i_cond = 0; i_cond < static_cast<int>(mrModelPart.NumberOfConditions()); ++i_cond){
        // iteration over all elements

        auto condition = mrModelPart.ConditionsBegin() + i_cond;

        if ( condition->Is(INLET) ){

            const int dim = condition->GetGeometry().Dimension();
            const Geometry<Node<3> >& p_geom = condition->GetGeometry();
            Vector Distance( p_geom.PointsNumber(), 0.0 );

            unsigned int negCount = 0;
            unsigned int posCount = 0;
            for (unsigned int i = 0; i < p_geom.PointsNumber(); i++){
                Distance[i] = condition->GetGeometry()[i].GetSolutionStepValue( DISTANCE );
                if ( condition->GetGeometry()[i].GetSolutionStepValue( DISTANCE ) > 0.0 ){
                    posCount++;
                } else {
                    negCount++;
                }
            }

            // leave the current iteration of the condition is completely on positive side
            if ( posCount == p_geom.PointsNumber() ){ continue; }

            if (dim == 2){      // 2D case: condition is a line

                array_1d<double, 3> normal;
                this->CalculateUnitNormal2D( normal, p_geom);
                if( norm_2( normal ) < epsilon ){ continue; }
                else { normal /= norm_2( normal ); }

                if ( negCount == p_geom.PointsNumber() ){       // the condition is completely on the negative side

                    // Gauss point information
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                    const unsigned int NumGauss = IntegrationPoints.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    p_geom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
                    const Matrix Ncontainer = p_geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                    for (unsigned int i_gauss = 0; i_gauss < NumGauss; i_gauss++){
                        const Vector N = row(Ncontainer, i_gauss);
                        double const wGauss = GaussPtsJDet[i_gauss] * IntegrationPoints[i_gauss].Weight();
                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        for (unsigned int n_node = 0; n_node < p_geom.PointsNumber(); n_node++){
                            interpolatedVelocity += N[n_node] * condition->GetGeometry()[n_node].GetSolutionStepValue(VELOCITY);
                        }
                        netInflow += - wGauss * inner_prod( normal, interpolatedVelocity );
                    }

                } else if ( negCount < p_geom.PointsNumber() && negCount < p_geom.PointsNumber() ){      // the condition is cut

                    // Compute the relative coordinate of the intersection point over the edge
                    const double aux_node_rel_location = std::abs ( Distance[0] /( Distance[1]-Distance[0] ));
                    // Shape function values at the position where the surface cuts the element
                    Vector Ncut = ZeroVector(2);
                    Ncut[0] = 1.0 - aux_node_rel_location;
                    Ncut[1] = aux_node_rel_location;

                    // Creation of an auxiliary geometry which covers the negative volume domain only
                    // (imitation of the splitting mechanism for triangles)
                    PointerVectorSet<IndexedPoint, IndexedObject> aux_point_container;
                    aux_point_container.reserve(2);
                    array_1d<double, 3> aux_point1_coords, aux_point2_coords, aux_velocity1, aux_velocity2;

                    IndexedPoint::Pointer paux_point1 = nullptr;
                    IndexedPoint::Pointer paux_point2 = nullptr;
                    // Modify the node with the positive distance
                    for (unsigned int i_node = 0; i_node < p_geom.PointsNumber(); i_node++){
                        if ( p_geom[i_node].GetSolutionStepValue( DISTANCE ) > 0.0 ){
                            // interpolating position for the new node
                            aux_point1_coords[0] = Ncut[0] * p_geom[0].X() + Ncut[1] * p_geom[1].X();
                            aux_point1_coords[1] = Ncut[0] * p_geom[0].Y() + Ncut[1] * p_geom[1].Y();
                            aux_point1_coords[2] = Ncut[0] * p_geom[0].Z() + Ncut[1] * p_geom[1].Z();

                            paux_point1 = Kratos::make_shared<IndexedPoint>(aux_point1_coords, mrModelPart.NumberOfNodes()+1);
                            aux_velocity1 = Ncut[0] * p_geom[0].GetSolutionStepValue( VELOCITY ) + Ncut[1] * p_geom[1].GetSolutionStepValue( VELOCITY );

                        } else {
                            aux_point2_coords[0] = p_geom[i_node].X();
                            aux_point2_coords[1] = p_geom[i_node].Y();
                            aux_point2_coords[2] = p_geom[i_node].Z();

                            paux_point2 = Kratos::make_shared<IndexedPoint>(aux_point2_coords, mrModelPart.NumberOfNodes()+2);
                            aux_velocity2 = p_geom[i_node].GetSolutionStepValue( VELOCITY );
                        }
                    }

                    Geometry <IndexedPoint>::Pointer p_aux_line = Kratos::make_shared< Line3D2 < IndexedPoint > >( paux_point1, paux_point2 );
                    // Gauss point information for auxiliary geometry
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_aux_line->IntegrationPoints( GeometryData::GI_GAUSS_2 );
                    const unsigned int NumGauss = IntegrationPoints.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    p_aux_line->DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
                    const Matrix Ncontainer = p_aux_line->ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                    for (unsigned int i_gauss = 0; i_gauss < NumGauss; i_gauss++){
                        const Vector N = row(Ncontainer, i_gauss);
                        double const wGauss = GaussPtsJDet[i_gauss] * IntegrationPoints[i_gauss].Weight();

                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        interpolatedVelocity = N[0] * aux_velocity1 + N[1] * aux_velocity2;
                        netInflow += - wGauss * inner_prod( normal, interpolatedVelocity );
                    }
                }
            }

            else if (dim == 3){
                // 3D case: condition is a triangle (implemented routines can be used)

                array_1d<double, 3> normal;
                this->CalculateUnitNormal3D( normal, p_geom);
                if( norm_2( normal ) < epsilon ){ continue; }
                else { normal /= norm_2( normal ); }

                if ( negCount == p_geom.PointsNumber() ){       // the condition is completely on the negative side

                    // Gauss point information
                    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geom.IntegrationPoints(GeometryData::GI_GAUSS_2);
                    const unsigned int NumGauss = IntegrationPoints.size();
                    Vector GaussPtsJDet = ZeroVector(NumGauss);
                    p_geom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::GI_GAUSS_2);
                    const Matrix Ncontainer = p_geom.ShapeFunctionsValues( GeometryData::GI_GAUSS_2 );

                    for (unsigned int i_gauss = 0; i_gauss < NumGauss; i_gauss++){
                        const Vector N = row(Ncontainer, i_gauss);
                        double const wGauss = GaussPtsJDet[i_gauss] * IntegrationPoints[i_gauss].Weight();
                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        for (unsigned int n_node = 0; n_node < p_geom.PointsNumber(); n_node++){
                            interpolatedVelocity += N[n_node] * condition->GetGeometry()[n_node].GetSolutionStepValue(VELOCITY);
                        }
                        netInflow += - wGauss * inner_prod( normal, interpolatedVelocity );
                    }

                } else if ( negCount < p_geom.PointsNumber() && negCount < p_geom.PointsNumber() ){      // the condition is cut

                    Matrix rShapeFunctions;
                    GeometryType::ShapeFunctionsGradientsType rShapeDerivatives;
                    Vector w_gauss_neg_side;
                    ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
                    auto geomPointer = condition->pGetGeometry();
                    p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(geomPointer, Distance);

                    p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                        rShapeFunctions,                    // N
                        rShapeDerivatives,                  // DN
                        w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                        GeometryData::GI_GAUSS_2);          // first order Gauss integration

                    // interating velocity over the negative area of the condition
                    for ( unsigned int i_gauss = 0; i_gauss < w_gauss_neg_side.size(); i_gauss++){
                        const array_1d<double,3>& N = row(rShapeFunctions, i_gauss);

                        array_1d<double,3> interpolatedVelocity = ZeroVector(3);
                        for (unsigned int n_node = 0; n_node < p_geom.PointsNumber(); n_node++){
                            interpolatedVelocity += N[n_node] * condition->GetGeometry()[n_node].GetSolutionStepValue(VELOCITY);
                        }
                        netInflow += - w_gauss_neg_side[i_gauss] * inner_prod( normal, interpolatedVelocity );
                    }
                }
            }
        }
    }
    return netInflow;
}


void MassConservationCheckProcess::ShiftDistanceField( double deltaDist ){

    ModelPart::NodesContainerType rNodes = mrModelPart.Nodes();
    #pragma omp parallel for
    for(int count = 0; count < static_cast<int>(rNodes.size()); count++)
    {
        ModelPart::NodesContainerType::iterator i_node = rNodes.begin() + count;
        i_node->FastGetSolutionStepValue( DISTANCE ) = i_node->FastGetSolutionStepValue( DISTANCE ) + deltaDist;
    }
}


/// Computes the 2D condition normal
/**
* @param An reference to condition normal vector
*/
void MassConservationCheckProcess::CalculateUnitNormal2D(array_1d<double,3>& An, const Geometry<Node<3> >& pGeometry)
{
    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.0;

}


/// Computes the 3D condition normal
/**
* @param An reference to condition normal vector
*/
void MassConservationCheckProcess::CalculateUnitNormal3D(array_1d<double,3>& An, const Geometry<Node<3> >& pGeometry)
{
    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}




}  // namespace Kratos.
