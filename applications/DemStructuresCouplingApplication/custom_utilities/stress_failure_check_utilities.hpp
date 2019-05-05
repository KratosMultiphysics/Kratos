//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#ifndef KRATOS_STRESS_FAILURE_CHECK_UTILITIES
#define KRATOS_STRESS_FAILURE_CHECK_UTILITIES

// System includes
#include <fstream>
#include <iostream>
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{

/// Note: For the moment this only works for cylindrical probes with its axis oriented in Z direction.

class StressFailureCheckUtilities
{

public:

KRATOS_CLASS_POINTER_DEFINITION(StressFailureCheckUtilities);

/// Default constructor.

StressFailureCheckUtilities(ModelPart& rModelPart,
                            Parameters& rParameters
                            ) :
                            mrModelPart(rModelPart)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
        {
            "cylinder_center": [0.0,0.0],
            "min_radius": 0.00381,
            "max_radius": 0.00481
        }  )" );

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mCylinderCenter[0] = rParameters["cylinder_center"][0].GetDouble();
    mCylinderCenter[1] = rParameters["cylinder_center"][1].GetDouble();
    mMinRadius = rParameters["min_radius"].GetDouble();
    mMaxRadius = rParameters["max_radius"].GetDouble();

    KRATOS_CATCH("");
}

/// Destructor.

virtual ~StressFailureCheckUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

// Before FEM and DEM solution
void ExecuteInitialize()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

// Before FEM and DEM solution
void ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

// After FEM and DEM solution
void ExecuteFinalizeSolutionStep()
{
    int NElems = static_cast<int>(mrModelPart.Elements().size());
    ModelPart::ElementsContainerType::iterator el_begin = mrModelPart.ElementsBegin();

    std::vector<double> Sigma1(OpenMPUtils::GetNumThreads(), 0.0);
    std::vector<double> Sigma2(OpenMPUtils::GetNumThreads(), 0.0);
    std::vector<double> Sigma3(OpenMPUtils::GetNumThreads(), 0.0);

    #pragma omp parallel for
    for(int i = 0; i < NElems; i++)
    {
        ModelPart::ElementsContainerType::iterator itElem = el_begin + i;

        double X = itElem->GetGeometry().GetPoint(0).X();
        double Y = itElem->GetGeometry().GetPoint(0).Y();


    }


/////////////////////

        // std::vector<double> thread_maxima(OpenMPUtils::GetNumThreads(), 0.0);
        // const int number_of_particles = (int) mListOfSphericContinuumParticles.size();

        // #pragma omp parallel for
        // for (int i = 0; i < number_of_particles; i++) {
        //     double max_sphere = mListOfSphericContinuumParticles[i]->CalculateMaxSearchDistance(has_mpi, r_process_info);
        //     if (max_sphere > thread_maxima[OpenMPUtils::ThisThread()]) thread_maxima[OpenMPUtils::ThisThread()] = max_sphere;
        // }

        // double maximum_across_threads = 0.0;
        // for (int i = 0; i < OpenMPUtils::GetNumThreads(); i++) {
        //     if (thread_maxima[i] > maximum_across_threads) maximum_across_threads = thread_maxima[i];
        // }

/////////////////////

            // BoundedMatrix<double, 3, 3> average_stress_tensor = ZeroMatrix(3,3);
            // for (int i = 0; i < 3; i++) {
            //     for (int j = 0; j < 3; j++) {
            //         average_stress_tensor(i,j) = 0.5 * ((*(element1->mSymmStressTensor))(i,j) + (*(element2->mSymmStressTensor))(i,j));
            //     }
            // }

            // Vector principal_stresses(3);
            // noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(average_stress_tensor);

            // Properties& element1_props = element1->GetProperties();
            // Properties& element2_props = element2->GetProperties();

            // const double mohr_coulomb_c = 1e6 * 0.5*(element1_props[INTERNAL_COHESION] + element2_props[INTERNAL_COHESION]);
            // const double mohr_coulomb_phi = 0.5 * (element1_props[INTERNAL_FRICTION_ANGLE] + element2_props[INTERNAL_FRICTION_ANGLE]);
            // const double mohr_coulomb_phi_in_radians = mohr_coulomb_phi * Globals::Pi / 180.0;
            // const double sinphi = std::sin(mohr_coulomb_phi_in_radians);
            // const double cosphi = std::cos(mohr_coulomb_phi_in_radians);

            // const double max_stress = *std::max_element(principal_stresses.begin(), principal_stresses.end());
            // const double min_stress = *std::min_element(principal_stresses.begin(), principal_stresses.end());

/////////////////////

        // std::fstream PropDataFile;
        // PropDataFile.open ("PropagationData.tcl", std::fstream::out | std::fstream::trunc);
        // PropDataFile.precision(12);
        // PropDataFile << "set PropagationData [list]" << std::endl;
        // PropDataFile << "return $PropagationData" << std::endl;
        // PropDataFile.close();

}



//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const
{
    return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const
{
}


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{

    ModelPart& mrModelPart;
    array_1d<double,2> mCylinderCenter;
    double mMinRadius;
    double mMaxRadius;

///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

private:

///@name Static Member r_variables
///@{


///@}
///@name Member r_variables
///@{
///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{


///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessible methods
///@{

/// Assignment operator.
StressFailureCheckUtilities & operator=(StressFailureCheckUtilities const& rOther);


///@}

}; // Class StressFailureCheckUtilities

}  // namespace Python.

#endif // KRATOS_STRESS_FAILURE_CHECK_UTILITIES
