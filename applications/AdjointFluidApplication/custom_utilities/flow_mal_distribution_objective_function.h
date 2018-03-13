//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION)
#define KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION

// System includes
#include <vector>
#include <string>
#include <fstream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_utilities/objective_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// An objective function for flow mal distribution.
template <unsigned int TDim>
class FlowMalDistributionObjectiveFunction : public ObjectiveFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(FlowMalDistributionObjectiveFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FlowMalDistributionObjectiveFunction(Parameters& rParameters)
    {
        KRATOS_TRY

        Parameters DefaultParamsSurfaceNormal(R"(
        {
            "objective_type": "flow_mal_distribution",
            "surface_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "flow_mal_distribution_direction": "surface_normal",
            "output_file": "flow_mal_distribution"
        })");

        Parameters DefaultParamsCustomDirection(R"(
        {
            "objective_type": "flow_mal_distribution",
            "surface_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "flow_mal_distribution_direction": [1.0, 0.0, 0.0],
            "output_file": "flow_mal_distribution"
        })");

        mSurfaceModelPartName = rParameters["surface_model_part_name"].GetString();
        mOutputFilename = rParameters["output_file"].GetString();

        if (rParameters["flow_mal_distribution_direction"].IsArray() == false)
        {
            rParameters.ValidateAndAssignDefaults(DefaultParamsSurfaceNormal);
            if (rParameters["flow_mal_distribution_direction"].GetString().compare("surface_normal")==0)
                mMalDistributionDirection = MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL;
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                               "flow_mal_distribution_direction only accepts as text inputs \"surface_normal\" or \"absolute_value\":",
                               rParameters.PrettyPrintJsonString())
        }
        else 
        {
            rParameters.ValidateAndAssignDefaults(DefaultParamsCustomDirection);
            if (rParameters["flow_mal_distribution_direction"].size() != 3)
            {
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "flow_mal_distribution_direction vector does not have size 3:",
                                   rParameters.PrettyPrintJsonString())
            }
            else
            {
                mMalDistributionDirection = MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION;
                for (unsigned int d = 0; d < TDim; ++d)
                    mDirection[d] = rParameters["flow_mal_distribution_direction"][d].GetDouble();

                if (std::abs(norm_2(mDirection) - 1.0) > 1e-3)
                {
                    const double magnitude = norm_2(mDirection);
                    if (magnitude == 0.0)
                        KRATOS_THROW_ERROR(std::runtime_error,
                                        "flow_mal_distribution_direction is not properly defined.",
                                        "")

                    std::cout
                        << "WARNING: non unit vector detected in \"flow_mal_distribution_direction\": "
                        << rParameters.PrettyPrintJsonString() << std::endl;
                    std::cout << "normalizing \"flow_mal_distribution_direction\"..." << std::endl;

                    for (unsigned int d = 0; d < TDim; d++)
                        mDirection[d] /= magnitude;
                }
            }
        }        

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~FlowMalDistributionObjectiveFunction()
    {
        if (mOutputFileOpenend)
            mOutputFileStream.close();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY

        if (rModelPart.HasSubModelPart(mSurfaceModelPartName) == false)
        {
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters \"surface_model_part_name\": ",
                mSurfaceModelPartName)
        }
        // initialize the variables to zero.
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

            for (auto it = NodesBegin; it != NodesEnd; ++it)
                it->Set(STRUCTURE, false);
        }

        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);
        // mark objective surface
        for (auto it = rSurfaceModelPart.NodesBegin();
             it != rSurfaceModelPart.NodesEnd();
             ++it)
            it->Set(STRUCTURE, true);
        
        if (mOutputFilename.compare("") != 0) 
        {
            mOutputFileStream.open(mOutputFilename + ".data");
            mOutputFileStream<<"#time       FLOW_MAL_DISTRIBUTION"<<std::endl;
            mOutputFileOpenend = true;
        }
        KRATOS_CATCH("")
    }

    virtual void InitializeSolutionStep(ModelPart& rModelPart)
    {
        KRATOS_TRY

        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        NormalCalculationUtils normal_calculation_obj = NormalCalculationUtils();
        const unsigned int domain_size = rModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
        normal_calculation_obj.CalculateOnSimplex(rSurfaceModelPart, domain_size);
        
        mObjectiveValue = Calculate(rModelPart);

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointVelocityContribution(const Element& rElem,
                                                      const Matrix& rAdjointMatrix,
                                                      Vector& rRHSContribution,
                                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        const unsigned int NumNodes = rElem.GetGeometry().PointsNumber();

        double inv_coeff = 1.0/((mN-1)*mObjectiveValue);

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElem.GetGeometry()[iNode].Is(STRUCTURE))
                {
                    const array_1d<double,3>& normal = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(NORMAL, 0);
                    double magnitude = norm_2(normal);
                    array_1d<double,3> normalized_normal = normal/magnitude;

                    const array_1d<double,3>& velocity = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, 0);
                    double directed_velocity = velocity[0]*normalized_normal[0] + 
                                               velocity[1]*normalized_normal[1] + 
                                               velocity[2]*normalized_normal[2];
                    
                    directed_velocity -= mAverageVelocity;

                    double coeff = inv_coeff * directed_velocity /
                                            rElem.GetGeometry()[iNode].GetValue(NODAL_AREA);
                    
                    for (unsigned int d = 0; d < TDim; d++) 
                        rRHSContribution[LocalIndex++] = coeff * normalized_normal[d];
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = 0.0;
                }

                //Adding pressure DOF                        
                rRHSContribution[LocalIndex++] = 0.0;
            }
            // std::cout<<rElem.GetId()<<rRHSContribution<<std::endl;            
        } 
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElem.GetGeometry()[iNode].Is(STRUCTURE))
                {
                    const array_1d<double,3>& velocity = rElem.GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, 0);
                    double directed_velocity = velocity[0]*mDirection[0] + 
                                               velocity[1]*mDirection[1] + 
                                               velocity[2]*mDirection[2];
                    
                    directed_velocity -= mAverageVelocity;

                    double coeff = inv_coeff * directed_velocity /
                                            rElem.GetGeometry()[iNode].GetValue(NODAL_AREA);

                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = coeff * mDirection[d];
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = 0.0;
                }

                //Adding pressure DOF                        
                rRHSContribution[LocalIndex++] = 0.0;
            }
        }

        
        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
                                                          const Matrix& rAdjointMatrix,
                                                          Vector& rRHSContribution,
                                                          ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    virtual double Calculate(ModelPart& rModelPart)
    {
        double result = 0.0;

        CalculateAverageObjectiveParameters(rModelPart);

        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            #pragma omp parallel reduction(+:result)
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rSurfaceModelPart.Nodes(), NodesBegin, NodesEnd);

                for (auto it = NodesBegin; it != NodesEnd; ++it)
                {
                    const array_1d<double,3>& normal = it->FastGetSolutionStepValue(NORMAL, 0);
                    double magnitude = norm_2(normal);
                    array_1d<double,3> normalized_normal = normal/magnitude;

                    const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);

                    result += pow(
                                    (
                                        normalized_normal[0]*velocity[0] + 
                                        normalized_normal[1]*velocity[1] + 
                                        normalized_normal[2]*velocity[2]
                                    ) - mAverageVelocity,
                                    2
                                );
                }
            }

        } 
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            #pragma omp parallel reduction(+:result)
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rSurfaceModelPart.Nodes(), NodesBegin, NodesEnd);

                for (auto it = NodesBegin; it != NodesEnd; ++it)
                {
                    const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);

                    result += pow(
                                    (
                                        mDirection[0]*velocity[0] + 
                                        mDirection[1]*velocity[1] + 
                                        mDirection[2]*velocity[2]
                                    ) - mAverageVelocity,
                                    2
                                );
                }
            }

        }

        return sqrt(result/(mN-1));
    }

    virtual void FinalizeSolutionStep(ModelPart& rModelPart)
    {
        mObjectiveValue = Calculate(rModelPart);

        if (mOutputFileOpenend) 
        {
            ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
            mOutputFileStream.precision(5);
            mOutputFileStream<<std::scientific<<rProcessInfo[TIME]<<" ";
            mOutputFileStream.precision(15);
            mOutputFileStream<<std::scientific<<mObjectiveValue<<std::endl;        
        }
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mSurfaceModelPartName;
    std::string mOutputFilename;
    std::ofstream mOutputFileStream;
    bool mOutputFileOpenend = false;
    array_1d<double, TDim> mDirection;
    int mMalDistributionDirection;

    double mObjectiveValue;
    double mAverageVelocity;
    int mN;

    const int MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL = 1;
    const int MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION = 2;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{   
    void CalculateAverageObjectiveParameters(ModelPart& rModelPart)
    {
        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        mAverageVelocity = 0.0;
        mN = 0;

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            for (auto it = rSurfaceModelPart.NodesBegin();
                it != rSurfaceModelPart.NodesEnd();
                ++it)
            {
                const array_1d<double,3>& normal = it->FastGetSolutionStepValue(NORMAL);
                double magnitude = norm_2(normal);
                array_1d<double,3> normalized_normal = normal/magnitude;
                const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);

                mAverageVelocity += (
                                        normalized_normal[0]*velocity[0] + 
                                        normalized_normal[1]*velocity[1] + 
                                        normalized_normal[2]*velocity[2]
                                    );
                mN++;
            }
        }
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            for (auto it = rSurfaceModelPart.NodesBegin();
                it != rSurfaceModelPart.NodesEnd();
                ++it)
            {
                const array_1d<double,3>& velocity = it->FastGetSolutionStepValue(VELOCITY, 0);


                mAverageVelocity += (
                                        mDirection[0]*velocity[0] + 
                                        mDirection[1]*velocity[1] + 
                                        mDirection[2]*velocity[2]
                                    );
                mN++;
            }
        }

        //Getting the average velocity by dividing flow rate by total objective surface area
        mAverageVelocity /= mN;
    }

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION defined */
