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

    if ( mCompareToInitial && time_step == 1){
        // getting initial values and storing as reference
        this->ComputeVolumesOfFluids( mInitialPositiveVolume, mInitialNegativeVolume );
    }

    if ( time_step % this->mMassComputationFreq == 0){

        // writing an output at a given frequncy
        this->ComputeVolumesOfFluids( mCurrentPositiveVolume, mCurrentNegativeVolume );
        mIsUpdated = true;
        return true;

    } else {

        mCurrentPositiveVolume = -1.0;
        mCurrentNegativeVolume = -1.0;
        return false;
    }
};


/* Private functions ****************************************************/

void MassConservationCheckProcess::ComputeVolumesOfFluids( double& posVolume, double& negVolume ){

    // initalisation (necessary because no reductio for type reference)
    double posVol = 0.0;
    double negVol = 0.0;

    #pragma omp parallel for reduction(+: posVol, negVol)
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        // iteration over all elements
        auto it_elem = mrModelPart.ElementsBegin() + i_elem;

        // useless containers
        Matrix rShapeFunctionsPos, rShapeFunctionsNeg;
        GeometryType::ShapeFunctionsGradientsType rShapeDerivativesPos, rShapeDerivativesNeg;

        auto p_geom = it_elem->pGetGeometry();
        unsigned int ptCountPos = 0;
        unsigned int ptCountNeg = 0;

        // instead of using data.isCut()
        for (unsigned int pt = 0; pt < p_geom->Points().size(); pt++){
            if ( p_geom->GetPoint(pt).FastGetSolutionStepValue(DISTANCE) > 0.0 ){ ptCountPos++; }
            else { ptCountNeg++; }
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
            // element is cut by the surface
            ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
            Vector w_gauss_pos_side(3, 0.0);
            Vector w_gauss_neg_side(3, 0.0);

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
                    rShapeFunctionsPos,                 // N
                    rShapeDerivativesPos,               // DN
                    w_gauss_pos_side,                   // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration (1 point per triangle)

            for ( unsigned int i = 0; i < w_gauss_pos_side.size(); i++){
                posVol += w_gauss_pos_side[i];
            }

            // Call the negative side modified shape functions calculator
            // Object p_modified_sh_func has full knowledge of slit geometry
            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                    rShapeFunctionsNeg,                 // N
                    rShapeDerivativesNeg,               // DN
                    w_gauss_neg_side,                   // includes the weights of the GAUSS points (!!!)
                    GeometryData::GI_GAUSS_1);          // first order Gauss integration

            for ( unsigned int i = 0; i < w_gauss_neg_side.size(); i++){
                negVol += w_gauss_neg_side[i];
            }
        }
    }

    posVolume = posVol;
    negVolume = negVol;
}

};  // namespace Kratos.