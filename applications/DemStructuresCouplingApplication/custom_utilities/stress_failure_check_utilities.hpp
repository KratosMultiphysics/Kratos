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

    mSigma1File.open("sigma1average_t.txt", std::fstream::out | std::fstream::trunc);
    mSigma1File << "0.0 0.0" << std::endl;
    mSigma1File.close();

    mSigma3File.open("sigma3average_t.txt", std::fstream::out | std::fstream::trunc);
    mSigma3File << "0.0 0.0" << std::endl;
    mSigma3File.close();

    KRATOS_CATCH("");
}

/// Destructor.

virtual ~StressFailureCheckUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

void ExecuteFinalizeSolutionStep()
{
    // int NElems = static_cast<int>(mrModelPart.Elements().size());
    // ModelPart::ElementsContainerType::iterator el_begin = mrModelPart.ElementsBegin();
    std::vector<double> ThreadSigma1(OpenMPUtils::GetNumThreads(), 0.0);
    // std::vector<double> ThreadSigma2(OpenMPUtils::GetNumThreads(), 0.0);
    std::vector<double> ThreadSigma3(OpenMPUtils::GetNumThreads(), 0.0);
    std::vector<int> ThreadNParticles(OpenMPUtils::GetNumThreads(), 0);

    #pragma omp parallel
    {
        int k = OpenMPUtils::ThisThread();

        #pragma omp for
        for(int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = mrModelPart.ElementsBegin() + i;

            const array_1d<double,3> DemPosition = itElem->GetGeometry().GetPoint(0).Coordinates();
            const double XDistance2 = std::pow(DemPosition[0]-mCylinderCenter[0],2);
            const double YDistance2 = std::pow(DemPosition[1]-mCylinderCenter[1],2);

            if( XDistance2 + YDistance2 > mMinRadius && XDistance2 + YDistance2 < mMaxRadius )
            {
                // BoundedMatrix<double, 3, 3> stress_tensor = (*(itElem->mSymmStressTensor));
                Vector principal_stresses(3);
                noalias(principal_stresses) = AuxiliaryFunctions::EigenValuesDirectMethod((*(itElem->mSymmStressTensor)));
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
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++)
    {
        Sigma1Average += ThreadSigma1[k];
        Sigma3Average += ThreadSigma3[k];
        NParticles += ThreadNParticles[k];
    }

    if(NParticles > 0)
    {
        Sigma1Average = Sigma1Average/NParticles;
        Sigma3Average = Sigma3Average/NParticles;
    }

    double CurrentTime = mrModelPart.GetProcessInfo()[TIME];

    mSigma1File.open("sigma1average_t.txt", std::fstream::out | std::fstream::app);
    mSigma1File << CurrentTime << " " << Sigma1Average << std::endl;
    mSigma1File.close();

    mSigma3File.open("sigma3average_t.txt", std::fstream::out | std::fstream::app);
    mSigma3File << CurrentTime << " " << Sigma3Average << std::endl;
    mSigma3File.close();
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
    std::fstream mSigma1File;
    std::fstream mSigma3File;

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
