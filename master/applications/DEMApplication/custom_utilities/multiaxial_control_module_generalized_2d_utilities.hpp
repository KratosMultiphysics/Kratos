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
//                   Guillermo Casas
//


#ifndef KRATOS_MULTIAXIAL_CONTROL_MODULE_GENERALIZED_2D_UTILITIES
#define KRATOS_MULTIAXIAL_CONTROL_MODULE_GENERALIZED_2D_UTILITIES

// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <map>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"

#include "includes/table.h"
#include "includes/kratos_parameters.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_elements/spheric_continuum_particle.h"

#include "DEM_application_variables.h"


namespace Kratos
{
class KRATOS_API(DEM_APPLICATION) MultiaxialControlModuleGeneralized2DUtilities
{
public:

KRATOS_CLASS_POINTER_DEFINITION(MultiaxialControlModuleGeneralized2DUtilities);

/// Definition of the index type
typedef std::size_t IndexType;

/// Defining a table with double argument and result type as table type.
typedef Table<double,double> TableType;

/// Default constructor.

MultiaxialControlModuleGeneralized2DUtilities(ModelPart& rDemModelPart, ModelPart& rFemModelPart, Parameters& rParameters);

/// Destructor.

virtual ~MultiaxialControlModuleGeneralized2DUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

// Before DEM solution
virtual void ExecuteInitialize();

// Before DEM solution
virtual void ExecuteInitializeSolutionStep();

// After DEM solution
virtual void ExecuteFinalizeSolutionStep();

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> string.

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

    ModelPart& mrDemModelPart;
    ModelPart& mrFemModelPart;
    double mCMDeltaTime;
    double mCMTime;
    double mKDeltaTime;
    double mKTime;
    unsigned int mStep;
    unsigned int mCMStep;
    unsigned int mKStep;
    unsigned int mActuatorCounter;
    double mPerturbationTolerance;
    unsigned int mPerturbationPeriod;
    double mMaxReactionCorrectionFraction;
    double mStiffnessAlpha;
    double mVelocityAlpha;
    double mReactionAlpha;
    double mMaximumInputStressVariationRate;
    std::vector<std::string> mVectorOfActuatorNames;
    std::map<std::string, std::vector<ModelPart*>> mListsOfFEMSubModelPartsForEachActuator; /// FEM SubModelParts associated to each boundary of every actuator
    std::map<std::string, std::vector<ModelPart*>> mListsOfDEMSubModelPartsForEachActuator; /// DEM SubModelParts associated to each boundary of every actuator
    std::map<std::string, std::vector<array_1d<double,3>>> mFEMOuterNormals; /// OuterNormal associated to each FEM boundary of every actuator
    // std::map<std::string, std::vector<array_1d<double,3>>> mDEMOuterNormals; /// OuterNormal associated to each DEM boundary of every actuator. TODO: not used for now...
    std::map<std::string, unsigned int> mTargetStressTableIds; /// TargetStressTableIds associated to every actuator
    Vector mVelocity;
    Vector mAcceleration; // TODO: possible CM enhancement
    Vector mReactionStress;
    Vector mReactionStressOld;
    Vector mDisplacement;
    Vector mDisplacementOld;
    Vector mElasticReactionStress;
    Matrix mStiffness;
    Matrix mDeltaDisplacement;
    Matrix mDeltaReactionStress;
    double mNormOfInitiallyEstimatedStiffness;
    double mMaxNodalVelocityForMultiDofs = 0.0;
    double mCompressionLengthForMultiDofs = 0.0;

///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{

virtual Vector MeasureReactionStress(const Variable<array_1d<double,3>>& rVariable);
Vector GetPerturbations(const Vector& rTargetStress, const double& rTime);
double GetConditionNumber(const Matrix& rInputMatrix, const Matrix& rInvertedMatrix);
void CalculateVelocity(const Vector& r_next_target_stress, const double& r_current_time);
void CalculateAcceleration(const Vector& r_next_target_stress, const double& r_current_time);
void CalculateStiffness();
void CalculateMaximumInputStressVariationRate();
void AddTableToSubModelPart(const unsigned int TableId, const Parameters TableParameters, ModelPart* pSubModelPart);
Parameters GetDefaultParametersForRadialActuator();
Parameters GetDefaultParametersForRadialMultiDofsActuator();
Parameters GetDefaultParametersForZActuator();
Parameters GetDefaultParametersForXOrYActuator();
void SetProvidedInitialVelocityToTheControlledBoundaries();

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
MultiaxialControlModuleGeneralized2DUtilities & operator=(MultiaxialControlModuleGeneralized2DUtilities const& rOther);


///@}

}; // Class MultiaxialControlModuleGeneralized2DUtilities

}  // namespace Python.

#endif // KRATOS_MULTIAXIAL_CONTROL_MODULE_GENERALIZED_2D_UTILITIES
