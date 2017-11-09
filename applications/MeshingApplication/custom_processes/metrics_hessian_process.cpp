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
#include "processes/compute_nodal_gradient_process.h" 
#include "custom_processes/metrics_hessian_process.h" 

namespace Kratos
{
template<unsigned int TDim, class TVarType>  
ComputeHessianSolMetricProcess<TDim, TVarType>::ComputeHessianSolMetricProcess(
        ModelPart& rThisModelPart,
        TVarType& rVariable,
        Parameters ThisParameters
        ):mThisModelPart(rThisModelPart),
          mVariable(rVariable)
{               
    Parameters DefaultParameters = Parameters(R"(
    {
        "minimal_size"                        : 0.1,
        "maximal_size"                        : 10.0, 
        "enforce_current"                     : true, 
        "hessian_strategy_parameters": 
        { 
            "interpolation_error"                  : 1.0e-6, 
            "mesh_dependent_constant"              : 0.28125
        }, 
        "anisotropy_remeshing"                : true, 
        "anisotropy_parameters":
        {
            "hmin_over_hmax_anisotropic_ratio"     : 1.0, 
            "boundary_layer_max_distance"          : 1.0, 
            "interpolation"                        : "Linear"
        }
    })" );
    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mEnforceCurrent = ThisParameters["enforce_current"].GetBool();
    
    // In case we have isotropic remeshing (default values)
    if (ThisParameters["anisotropy_remeshing"].GetBool() == false)
    {
        mInterpError = DefaultParameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble();
        mMeshConstant = DefaultParameters["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble();
        mAnisRatio = DefaultParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = DefaultParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(DefaultParameters["anisotropy_parameters"]["interpolation"].GetString());
    }
    else
    {
        mInterpError = ThisParameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble();
        mMeshConstant = ThisParameters["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble();
        mAnisRatio = ThisParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = ThisParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(ThisParameters["anisotropy_parameters"]["interpolation"].GetString());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, class TVarType>  
void ComputeHessianSolMetricProcess<TDim, TVarType>::Execute()
{
    // Iterate in the nodes
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    int num_nodes = nodes_array.end() - nodes_array.begin();
    
    CalculateAuxiliarHessian();
    
    #pragma omp parallel for 
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
    #ifdef KRATOS_DEBUG 
        if ( it_node->SolutionStepsDataHas( mVariable ) == false )
        {
            KRATOS_ERROR << "Missing variable on node " << it_node->Id() << std::endl;
        }
    #endif
        
        const double distance = it_node->FastGetSolutionStepValue(DISTANCE); // TODO: This should be changed for the varaible of interestin the future. This means that the value of the boundary value would be changed to a threshold value instead
        const Vector& hessian = it_node->GetValue(AUXILIAR_HESSIAN);

        const double nodal_h = it_node->FastGetSolutionStepValue(NODAL_H);            
        
        double element_min_size = mMinSize;
        if ((element_min_size > nodal_h) && (mEnforceCurrent == true)) element_min_size = nodal_h;
        double element_max_size = mMaxSize;
        if ((element_max_size > nodal_h) && (mEnforceCurrent == true)) element_max_size = nodal_h;

        const double ratio = CalculateAnisotropicRatio(distance, mAnisRatio, mBoundLayer, mInterpolation);
        
        // For postprocess pourposes
        it_node->SetValue(ANISOTROPIC_RATIO, ratio); 
        
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
            const Vector& new_metric = ComputeHessianMetricTensor(hessian, ratio, element_min_size, element_max_size);    
            
            metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
        }
        else
        {
            metric = ComputeHessianMetricTensor(hessian, ratio, element_min_size, element_max_size);    
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, class TVarType>  
Vector ComputeHessianSolMetricProcess<TDim, TVarType>::ComputeHessianMetricTensor(
    const Vector& Hessian,
    const double& AnisotropicRatio,
    const double& ElementMinSize, // This way we can impose as minimum as the previous size if we desire
    const double& ElementMaxSize // This way we can impose as maximum as the previous size if we desire
    )
{        
    // Calculating Metric parameters
    const double CEpsilon = mMeshConstant/mInterpError;
    const double min_ratio = 1.0/(ElementMinSize * ElementMinSize);
//         const double min_ratio = 1.0/(mMinSize * mMinSize);
    const double max_ratio = 1.0/(ElementMaxSize * ElementMaxSize);
//         const double max_ratio = 1.0/(mMaxSize * mMaxSize);
    
    typedef bounded_matrix<double, TDim, TDim> TempType;
    
    // Declaring the eigen system
    bounded_matrix<double, TDim, TDim> eigen_vector_matrix;
    bounded_matrix<double, TDim, TDim> eigen_values_matrix;

    // We first transform into a matrix
    const bounded_matrix<double, TDim, TDim> hessian_matrix = MetricsMathUtils<TDim>::VectorToTensor(Hessian);
    
    MathUtils<double>::EigenSystem<TDim>(hessian_matrix, eigen_vector_matrix, eigen_values_matrix, 1e-18, 20);
    
    // Recalculate the Metric eigen values
    for (unsigned int i = 0; i < TDim; ++i)
    {
        eigen_values_matrix(i, i) = MathUtils<double>::Min(MathUtils<double>::Max(CEpsilon * std::abs(eigen_values_matrix(i, i)), max_ratio), min_ratio);
    }
    
    // Considering anisotropic
    if (AnisotropicRatio < 1.0)
    {
        double eigen_max = eigen_values_matrix(0, 0);
        double eigen_min = eigen_values_matrix(1, 1);
        for (unsigned int i = 1; i < TDim - 1; ++i)
        {
            eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
            eigen_min = MathUtils<double>::Min(eigen_max, eigen_values_matrix(i, i));
        }
        
        const double eigen_radius = std::abs(eigen_max - eigen_min) * (1.0 - AnisotropicRatio);
        const double relative_eigen_radius = std::abs(eigen_max - eigen_radius);
        
        for (unsigned int i = 0; i < TDim; ++i)
        {
            eigen_values_matrix(i, i) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(i, i), eigen_max), relative_eigen_radius);
        }
    }
    else // NOTE: For isotropic we should consider the maximum of the eigenvalues
    {
        double eigen_max = eigen_values_matrix(0, 0);
        for (unsigned int i = 1; i < TDim - 1; ++i)
        {
            eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
        }
        for (unsigned int i = 0; i < TDim; ++i)
        {
            eigen_values_matrix(i, i) = eigen_max;
        }
        eigen_vector_matrix = IdentityMatrix(TDim, TDim);
    }
        
    // We compute the product
    const bounded_matrix<double, TDim, TDim>& matric_matrix =  prod(trans(eigen_vector_matrix), prod<TempType>(eigen_values_matrix, eigen_vector_matrix));
    
    // Finally we transform to a vector
    const Vector& metric = MetricsMathUtils<TDim>::TensorToVector(matric_matrix);
    
    return metric;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, class TVarType>  
void ComputeHessianSolMetricProcess<TDim, TVarType>::CalculateAuxiliarHessian()
{
    // Iterate in the nodes
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    int num_nodes = nodes_array.end() - nodes_array.begin();
    
    // Declaring auxiliar vector
    const Vector aux_zero_vector = ZeroVector(3 * (TDim - 1));
    
    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        
        it_node->SetValue(AUXILIAR_HESSIAN, aux_zero_vector);  
    }
    
    // Compute auxiliar gradient
    ComputeNodalGradientProcess<TDim, TVarType, NonHistorical> gradient_process = ComputeNodalGradientProcess<TDim, TVarType, NonHistorical>(mThisModelPart, mVariable, AUXILIAR_GRADIENT, NODAL_AREA);
    gradient_process.Execute();
    
    // Iterate in the conditions
    ElementsArrayType& elements_array = mThisModelPart.Elements();
    int num_elements = elements_array.end() - elements_array.begin();
    
    #pragma omp parallel for
    for(int i = 0; i < num_elements; ++i) 
    {
        auto it_elem = elements_array.begin() + i;
        
        Element::GeometryType& geom = it_elem->GetGeometry();

        double Volume;
        if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
        {
            bounded_matrix<double,3, 2> DN_DX;
            array_1d<double, 3> N;

            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            
            bounded_matrix<double,3, 2> values;
            for(unsigned int i_node = 0; i_node < 3; ++i_node)
            {
                const array_1d<double, 3> aux_grad = geom[i_node].GetValue(AUXILIAR_GRADIENT);
                values(i_node, 0) = aux_grad[0];
                values(i_node, 1) = aux_grad[1];
            }
            
            const bounded_matrix<double,2, 2>& hessian = prod(trans(DN_DX), values); 
            const Vector& hessian_cond = MetricsMathUtils<2>::TensorToVector(hessian);
            
            for(unsigned int i_node = 0; i_node < geom.size(); ++i_node)
            {
                for(unsigned int k = 0; k < 3; k++)
                {
                    double& val = geom[i_node].GetValue(AUXILIAR_HESSIAN)[k];
                    
                    #pragma omp atomic
                    val += N[i_node] * Volume * hessian_cond[k];
                }
            }
        }
        else if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4)
        {
            bounded_matrix<double,4,  3> DN_DX;
            array_1d<double, 4> N;
            
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
            
            bounded_matrix<double,4, 3> values;
            for(unsigned int i_node = 0; i_node < 4; ++i_node)
            {
                const array_1d<double, 3> aux_grad = geom[i_node].GetValue(AUXILIAR_GRADIENT);
                values(i_node, 0) = aux_grad[0];
                values(i_node, 1) = aux_grad[1];
                values(i_node, 2) = aux_grad[2];
            }
            
            const bounded_matrix<double, 3, 3> hessian = prod(trans(DN_DX), values); 
            const Vector& hessian_cond = MetricsMathUtils<3>::TensorToVector(hessian);
            
            for(unsigned int i_node = 0; i_node < geom.size(); ++i_node)
            {
                for(unsigned int k = 0; k < 6; k++)
                {
                    double& val = geom[i_node].GetValue(AUXILIAR_HESSIAN)[k];
                    
                    #pragma omp atomic
                    val += N[i_node] * Volume * hessian_cond[k];
                }
            }
        }
        else
        {
            KRATOS_ERROR << "WARNING: YOU CAN USE JUST 2D TRIANGLES OR 3D TETRAEDRA RIGHT NOW IN THE GEOMETRY UTILS: " << geom.size() << std::endl;
        }
    }
        
    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) 
    {
        auto it_node = nodes_array.begin() + i;
        it_node->GetValue(AUXILIAR_HESSIAN) /= it_node->FastGetSolutionStepValue(NODAL_AREA);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, class TVarType>  
Interpolation ComputeHessianSolMetricProcess<TDim, TVarType>::ConvertInter(const std::string& str)
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

template<unsigned int TDim, class TVarType>  
double ComputeHessianSolMetricProcess<TDim, TVarType>::CalculateAnisotropicRatio(
    const double& distance,
    const double& rAnisRatio,
    const double& rBoundLayer,
    const Interpolation& rInterpolation
    )
{
    const double tolerance = 1.0e-12;
    double ratio = 1.0; // NOTE: Isotropic mesh
    if (rAnisRatio < 1.0)
    {                           
        if (std::abs(distance) <= rBoundLayer)
        {
            if (rInterpolation == Constant)
            {
                ratio = rAnisRatio;
            }
            else if (rInterpolation == Linear)
            {
                ratio = rAnisRatio + (std::abs(distance)/rBoundLayer) * (1.0 - rAnisRatio);
            }
            else if (rInterpolation == Exponential)
            {
                ratio = - std::log(std::abs(distance)/rBoundLayer) * rAnisRatio + tolerance;
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

template class ComputeHessianSolMetricProcess<2, Variable<double>>;
template class ComputeHessianSolMetricProcess<3, Variable<double>>;
template class ComputeHessianSolMetricProcess<2, ComponentType>;
template class ComputeHessianSolMetricProcess<3, ComponentType>;

};// namespace Kratos.
