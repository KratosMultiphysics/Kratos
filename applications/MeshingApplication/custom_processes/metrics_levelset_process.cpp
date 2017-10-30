// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// Project includes
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "custom_processes/metrics_levelset_process.h"

namespace Kratos
{
template<unsigned int TDim>  
ComputeLevelSetSolMetricProcess<TDim>::ComputeLevelSetSolMetricProcess(
        ModelPart& rThisModelPart,
        const Variable<array_1d<double,3>> rVariableGradient,
        Parameters ThisParameters
        ):mThisModelPart(rThisModelPart),
          mVariableGradient(rVariableGradient)
{   
    Parameters DefaultParameters = Parameters(R"(
    {
        "minimal_size"                         : 0.1, 
        "enforce_current"                      : true, 
        "anisotropy_remeshing"                 : true, 
        "anisotropy_parameters": 
        {
            "hmin_over_hmax_anisotropic_ratio"      : 1.0, 
            "boundary_layer_max_distance"           : 1.0, 
            "interpolation"                         : "Linear"
        }
    })" );
    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
    
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mEnforceCurrent = ThisParameters["enforce_current"].GetBool();
    
    // In case we have isotropic remeshing (default values)
    if (ThisParameters["anisotropy_remeshing"].GetBool() == false)
    {
        mAnisRatio = DefaultParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = DefaultParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(DefaultParameters["anisotropy_parameters"]["interpolation"].GetString());
    }
    else
    {
        mAnisRatio = ThisParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = ThisParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(ThisParameters["anisotropy_parameters"]["interpolation"].GetString());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>  
void ComputeLevelSetSolMetricProcess<TDim>::Execute()
{
    // Iterate in the nodes
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    int num_nodes = nodes_array.end() - nodes_array.begin();
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; i++) 
    {
        auto it_node = nodes_array.begin() + i;
        
        if ( it_node->SolutionStepsDataHas( mVariableGradient ) == false )
        {
            KRATOS_ERROR << "Missing gradient variable on node " << it_node->Id() << std::endl;
        }
        
        const double distance = it_node->FastGetSolutionStepValue(DISTANCE);
        array_1d<double, 3> gradient_value = it_node->FastGetSolutionStepValue(mVariableGradient, 0);
        
        const double ratio = CalculateAnisotropicRatio(distance, mAnisRatio, mBoundLayer, mInterpolation);
        
        // For postprocess pourposes
        it_node->SetValue(ANISOTROPIC_RATIO, ratio); 
        
        double element_size = mMinSize;
        const double& NodalH = it_node->FastGetSolutionStepValue(NODAL_H, 0);
        if (((element_size > NodalH) && (mEnforceCurrent == true)) || (std::abs(distance) > mBoundLayer))
        {
            element_size = NodalH;
        }
        
        const double tolerance = 1.0e-12;
        const double norm_gradient_value = norm_2(gradient_value);
        if (norm_gradient_value > tolerance)
        {
            gradient_value /= norm_gradient_value;
        }
        
        // We compute the metric
    #ifdef KRATOS_DEBUG 
        if( it_node->Has(MMG_METRIC) == false) 
        {
            KRATOS_ERROR <<  " MMG_METRIC not defined for node " << it_node->Id();
        }
    #endif     
        Vector& metric = it_node->GetValue(MMG_METRIC);
        
    #ifdef KRATOS_DEBUG 
        if(metric.size() != TDim * 3 - 3) 
        {
            KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << it_node->Id() << " size is " << metric.size() << " expected size was " << TDim * 3 - 3;
        }
    #endif
        
        const double norm_metric = norm_2(metric);
        if (norm_metric > 0.0) // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
        {
            const Vector& old_metric = it_node->GetValue(MMG_METRIC);
            const Vector& new_metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
            
            metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
        }
        else
        {
            metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
Vector ComputeLevelSetSolMetricProcess<2>::ComputeLevelSetMetricTensor(
    const array_1d<double, 3>& GradientValue,
    const double& Ratio,
    const double& ElementSize
)
{
    Vector metric;
    metric.resize(3, false);
    
    const double Coeff0 = 1.0/(ElementSize * ElementSize);
    const double Coeff1 = Coeff0/(Ratio * Ratio);
    
    const double v0v0 = GradientValue[0]*GradientValue[0];
    const double v0v1 = GradientValue[0]*GradientValue[1];
    const double v1v1 = GradientValue[1]*GradientValue[1];
    
    metric[0] = Coeff0*(1.0 - v0v0) + Coeff1*v0v0;
    metric[1] = Coeff0*(    - v0v1) + Coeff1*v0v1;  
    metric[2] = Coeff0*(1.0 - v1v1) + Coeff1*v1v1;
    
    return metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<>  
Vector ComputeLevelSetSolMetricProcess<3>::ComputeLevelSetMetricTensor(
    const array_1d<double, 3>& GradientValue,
    const double& Ratio,
    const double& ElementSize
)
{
    Vector metric;
    metric.resize(6, false);
    
    const double Coeff0 = 1.0/(ElementSize * ElementSize);
    const double Coeff1 = Coeff0/(Ratio * Ratio);
    
    const double v0v0 = GradientValue[0]*GradientValue[0];
    const double v0v1 = GradientValue[0]*GradientValue[1];
    const double v0v2 = GradientValue[0]*GradientValue[2];
    const double v1v1 = GradientValue[1]*GradientValue[1];
    const double v1v2 = GradientValue[1]*GradientValue[2];
    const double v2v2 = GradientValue[2]*GradientValue[2];
    
    metric[0] = Coeff0*(1.0 - v0v0) + Coeff1*v0v0;
    metric[1] = Coeff0*(    - v0v1) + Coeff1*v0v1; 
    metric[2] = Coeff0*(    - v0v2) + Coeff1*v0v2; 
    metric[3] = Coeff0*(1.0 - v1v1) + Coeff1*v1v1; 
    metric[4] = Coeff0*(    - v1v2) + Coeff1*v1v2; 
    metric[5] = Coeff0*(1.0 - v2v2) + Coeff1*v2v2;

    return metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>  
Interpolation ComputeLevelSetSolMetricProcess<TDim>::ConvertInter(const std::string& str)
{
    if(str == "Constant") 
    {
        return Constant;
    }
    else if(str == "Linear") 
    {
        return Linear;
    }
    else if(str == "Exponential") 
    {
        return Exponential;
    }
    else
    {
        return Linear;
    }
}
    
/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>  
double ComputeLevelSetSolMetricProcess<TDim>::CalculateAnisotropicRatio(
    const double& Distance,
    const double& rAnisRatio,
    const double& rBoundLayer,
    const Interpolation& rInterpolation
    )
{
    const double tolerance = 1.0e-12;
    double ratio = 1.0; // NOTE: Isotropic mesh
    if (rAnisRatio < 1.0)
    {                           
        if (std::abs(Distance) <= rBoundLayer)
        {
            if (rInterpolation == Constant)
            {
                ratio = rAnisRatio;
            }
            else if (rInterpolation == Linear)
            {
                ratio = rAnisRatio + (std::abs(Distance)/rBoundLayer) * (1.0 - rAnisRatio);
            }
            else if (rInterpolation == Exponential)
            {
                ratio = - std::log(std::abs(Distance)/rBoundLayer) * rAnisRatio + tolerance;
                if (ratio > 1.0)
                {
                    ratio = 1.0;
                }
            }
        }
    }
    
    return ratio;
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeLevelSetSolMetricProcess<2>;
template class ComputeLevelSetSolMetricProcess<3>;

};// namespace Kratos.
