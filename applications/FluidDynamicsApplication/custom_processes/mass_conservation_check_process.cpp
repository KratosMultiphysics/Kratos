//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    under supervision of Ruben Zorrilla
//
//

// System includes
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_h_process.h"

#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "mass_conservation_check_process.h"
#include "fluid_dynamics_application_variables.h"

#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{

/* Public functions *******************************************************/

MassConservationCheckProcess::MassConservationCheckProcess(
        ModelPart& rModelPart, 
        const bool IsSwitchedOn,
        const int MassComputationFreq,
        const bool CompareToInitial, 
        const bool WriteToLogFile)
    : Process(), mrModelPart(rModelPart) {

    mIsSwitchedOn = IsSwitchedOn;
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
        "is_switched_on"                         : true,
        "mass_computation_frequency"             : 20,
        "compare_to_initial_values"              : true,
        "write_to_log_file"                      : true
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mIsSwitchedOn = rParameters["is_switched_on"].GetBool();
    mMassComputationFreq = rParameters["mass_computation_frequency"].GetInt();
    mCompareToInitial = rParameters["compare_to_initial_values"].GetBool();
    mWriteToLogFile = rParameters["write_to_log_file"].GetBool();
}

// void DistanceModificationProcess::ExecuteInitialize() {

//     KRATOS_TRY;

//     // Continuous distance field required variables check
//     if (mContinuousDistance){
//         const auto& r_node = *mrModelPart.NodesBegin();
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_H, r_node);
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);
//     }

//     KRATOS_CATCH("");
// }

void MassConservationCheckProcess::ExecuteFinalizeSolutionStep()
{
  
    if ( this->mIsSwitchedOn ){
            double negVol = 0.0;
            double posVol = 0.0;

        if ( mCompareToInitial && mCounter == 1){
            // getting initial values and storing as reference
            this->ComputeVolumesOfFluids( posVol, negVol );
            mInitialPositiveVolume = posVol;
            mInitialNegativeVolume = negVol;
        }

        if ( mCounter % this->mMassComputationFreq == 0){
            // writing an output at a given frequncy
            this->ComputeVolumesOfFluids( posVol, negVol );
            std::cout << " --- Volume Checking Process --- " << std::endl;
            std::cout << "   positive Volume = " << posVol << std::endl;
            if ( mCompareToInitial ){ 
                std::cout << "   ( " << posVol / mInitialPositiveVolume * 100.0 << "% of initial )" << std::endl; 
            }
            std::cout << "   negative Volume = " << negVol << std::endl;
            if ( mCompareToInitial ){ 
                std::cout << "   ( " << negVol / mInitialNegativeVolume * 100.0 << "% of initial )" << std::endl; 
            }
            std::cout << " --- --- --- --- --- --- --- --- " << std::endl;

            if (mWriteToLogFile){
                std::ofstream myfile;
                myfile.open ("volumeLog.txt", std::ios::out | std::ios::app);
                if ( myfile.is_open() ){

                    myfile << mCounter << "    ";
                    myfile << posVol << "    ";
                    if ( mCompareToInitial ){ 
                        myfile << posVol / mInitialPositiveVolume * 100.0 << "    "; 
                    }
                    myfile << negVol << "    ";
                    if ( mCompareToInitial ){ 
                        myfile << negVol / mInitialNegativeVolume * 100.0 << "\n"; 
                    } 
                    myfile.close();
                }
            }
        }

        mCounter++;
    } else {
        return;
    }
};


/* Private functions ****************************************************/

void MassConservationCheckProcess::ComputeVolumesOfFluids( double &positiveVolume, double &negativeVolume )
{
    // useless containers
    Matrix rShapeFunctionsPos, rShapeFunctionsNeg;
    GeometryType::ShapeFunctionsGradientsType rShapeDerivativesPos, rShapeDerivativesNeg;
    Vector null3(3, 0.0);
    Vector null4(4, 0.0);

    // abbreviation of variable names
    double Vol;
    double& posVol = positiveVolume;
    double& negVol = negativeVolume;

    // initalisation
    posVol = 0.0;
    negVol = 0.0;

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

        if ( ptCountPos == p_geom->Points().size() ){
            // all nodes are positive
            if( p_geom->Points().size() == 3 ){
                Vol = p_geom->Area();
            } else {
                Vol = p_geom->Volume();
            }
            posVol += Vol;
        } 
        else if ( ptCountNeg == p_geom->Points().size() ){
            // all nodes are negative
            if( p_geom->Points().size() == 3 ){
                Vol = p_geom->Area();
            } else {
                Vol = p_geom->Volume();
            }
            negVol += Vol;
        }
        else {
            // element is cut by the surface
            ModifiedShapeFunctions::Pointer p_modified_sh_func = nullptr;
            Vector w_gauss_pos_side(3, 0.0);
            Vector w_gauss_neg_side(3, 0.0);

            // the element is cut by the interfae between the fluids
            if( p_geom->Points().size() == 3 ){
                Vector& Distance = null3;
                for (unsigned int i = 0; i < 3; i++){
                    Distance[i] = p_geom->GetPoint(i).FastGetSolutionStepValue(DISTANCE);
                }
                // the current element is a triangle (2D)
                // constructor executes splitting of geometry based on distance fields
                p_modified_sh_func = Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geom, Distance);
            }
            else if ( p_geom->Points().size() == 4 ){
                Vector& Distance = null4;
                for (unsigned int i = 0; i < 4; i++){
                    Distance[i] = p_geom->GetPoint(i).FastGetSolutionStepValue(DISTANCE);
                }
                // the current element is a tetrahedron (3D)
                // constructor executes splitting of geometry based on distance fields
                p_modified_sh_func = Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(p_geom, Distance);
            } 
            else {
                KRATOS_ERROR << "The process can not be applied on this kind of element (MODIFY MESSAGE!)" << std::endl;
            }

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
}


};  // namespace Kratos.
