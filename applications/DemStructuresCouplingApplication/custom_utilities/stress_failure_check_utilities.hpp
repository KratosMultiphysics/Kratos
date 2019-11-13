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
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_elements/spheric_continuum_particle.h"

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
            "cylinder_center": [0.0,0.0,0.0],
            "min_radius": 0.00381,
            "max_radius": 0.00481
        }  )" );

    // Now validate agains defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mCylinderCenter[0] = rParameters["cylinder_center"][0].GetDouble();
    mCylinderCenter[1] = rParameters["cylinder_center"][1].GetDouble();
    mCylinderCenter[2] = rParameters["cylinder_center"][2].GetDouble();
    mMinRadius = rParameters["min_radius"].GetDouble();
    mMaxRadius = rParameters["max_radius"].GetDouble();

    KRATOS_CATCH("");
}

/// Destructor.

virtual ~StressFailureCheckUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

void ExecuteFinalizeSolutionStep()
{
    std::vector<double> ThreadSigma1(OpenMPUtils::GetNumThreads(), 0.0);
    // std::vector<double> ThreadSigma2(OpenMPUtils::GetNumThreads(), 0.0);
    std::vector<double> ThreadSigma3(OpenMPUtils::GetNumThreads(), 0.0);
    std::vector<int> ThreadNParticles(OpenMPUtils::GetNumThreads(), 0);

    ModelPart::ElementsContainerType& rElements = mrModelPart.GetCommunicator().LocalMesh().Elements();

    #pragma omp parallel
    {
        int k = OpenMPUtils::ThisThread();

        #pragma omp for
        for (int i = 0; i < (int)rElements.size(); i++) {

            ModelPart::ElementsContainerType::ptr_iterator ptr_itElem = rElements.ptr_begin() + i;

            const array_1d<double,3> DemPosition = (*ptr_itElem)->GetGeometry()[0].Coordinates();
            const double Distance2 = std::pow(DemPosition[0] - mCylinderCenter[0], 2) + std::pow(DemPosition[1] - mCylinderCenter[1], 2);

            Element* p_element = ptr_itElem->get();
            SphericContinuumParticle* pDemElem = dynamic_cast<SphericContinuumParticle*>(p_element);

            if ((pDemElem->IsNot(DEMFlags::STICKY)) && (Distance2 >= std::pow(mMinRadius,2)) && (Distance2 <= std::pow(mMaxRadius,2))) {

                BoundedMatrix<double, 3, 3> stress_tensor = (*(pDemElem->mSymmStressTensor));
                Vector principal_stresses(3);
                noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod(stress_tensor);
                const double max_stress = *std::max_element(principal_stresses.begin(), principal_stresses.end());
                const double min_stress = *std::min_element(principal_stresses.begin(), principal_stresses.end());
                ThreadSigma1[k] += max_stress;
                ThreadSigma3[k] += min_stress;
                ThreadNParticles[k] += 1;
            }
        }
    }

    double Sigma1Average = 0.0;
    double Sigma3Average = 0.0;
    int NParticles = 0;
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++) {
        Sigma1Average += ThreadSigma1[k];
        Sigma3Average += ThreadSigma3[k];
        NParticles += ThreadNParticles[k];
    }

    if (NParticles > 0) {
        Sigma1Average = Sigma1Average / NParticles;
        Sigma3Average = Sigma3Average / NParticles;
    }

    double CurrentTime = mrModelPart.GetProcessInfo()[TIME];
    mrModelPart.GetProcessInfo()[SIGMA_3_AVERAGE] = Sigma3Average;

    // Note: These two files are written in the "problemname_Graphs" folder. They are erased at the beginning of the simulation.
    std::fstream Sigma1File;
    std::fstream Sigma3File;

    Sigma1File.open("sigma1average_t.txt", std::fstream::out | std::fstream::app);
    Sigma1File.precision(12);
    Sigma1File << CurrentTime << " " << Sigma1Average << std::endl;
    Sigma1File.close();

    Sigma3File.open("sigma3average_t.txt", std::fstream::out | std::fstream::app);
    Sigma3File.precision(12);
    Sigma3File << CurrentTime << " " << Sigma3Average << std::endl;
    Sigma3File.close();
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
    array_1d<double,3> mCylinderCenter;
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
