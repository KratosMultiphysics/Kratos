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

#if !defined(KRATOS_METRICS_UTILITY)
#define KRATOS_METRICS_UTILITY

// Project includes
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "processes/compute_nodal_gradient_process.h" // TODO: Not prism or quadrilaterals implemented yet

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

//// This class is can be used to compute

template<unsigned int TDim>  
class MetricsUtility
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;
    typedef Node <3>                                                   NodeType;
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart: The model part to be computed
     * @param rMinSize: The min size of element
     * @param rAnisRatio: The anisotropic ratio
     * @param rBoundLayer: The boundary layer limit
     * @param rInterpolation: The interpolation type
     */
    
    MetricsUtility(
        const double rMinSize,
        const double rAnisRatio = 1.0,
        const double rBoundLayer =  1.0,
        const std::string rInterpolation = "Linear"
        )
        :mMinSize(rMinSize),
        mAnisRatio(rAnisRatio),
        mBoundLayer(rBoundLayer),
        mInterpolation(rInterpolation)
    {       
        std::cout << "Initializing metric utility" << std::endl;
    }
    
    /// Destructor.
    ~MetricsUtility() {}
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * We initialize the metrics of the MMG sol using a level set approach
     * @param rVariableGradient: The gradient variable to compute
     */
    
    void ComputeLevelSetSolMetric(
        ModelPart& rThisModelPart,
        const Variable<array_1d<double,3>> rVariableGradient = DISTANCE_GRADIENT
        )
    {
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if ( itNode->SolutionStepsDataHas( rVariableGradient ) == false )
            {
                KRATOS_THROW_ERROR( std::invalid_argument, "Missing gradient variable on node ", itNode->Id() )
            }
            
            const double distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            array_1d<double, 3> gradient_value = itNode->FastGetSolutionStepValue(rVariableGradient, 0);
            
            double element_size = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (element_size > mMinSize)
            {
                element_size = mMinSize;
            }
            
            const double ratio = CalculateAnisotropicRatio(distance);
            
            // For postprocess pourposes
            double& anisotropic_ratio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            anisotropic_ratio = ratio;
            
            const double tolerance = 1.0e-12;
            const double norm = norm_2(gradient_value);
            if (norm > tolerance)
            {
                gradient_value /= norm;
            }
            
            // We compute the metric
            #ifdef KRATOS_DEBUG 
            if( itNode->Has(MMG_METRIC) == false) 
            {
                KRATOS_ERROR <<  " MMG_METRIC not defined for node " << itNode>Id();
            }
            #endif     
            Vector& metric = itNode->GetValue(MMG_METRIC);
            
            #ifdef KRATOS_DEBUG 
            if(metric.size() != TDim * 3 - 3) 
            {
                KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << itNode>Id() << " size is " << metric.size() << " expected size was " << TDim * 3 - 3;
            }
            #endif
            
            const double normmetric = norm_2(metric);
            if (normmetric > 0.0) // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            {
                const Vector old_metric = itNode->GetValue(MMG_METRIC);
                const Vector new_metric = ComputeLevelSetMetricTensor(gradient_value, ratio, element_size);
                
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
    
    /**
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     * @param rMaxSize: The maximal size of the elements
     * @param rVariable: The variable considered for remeshing
     * @param rInterpError: The interpolation error assumed
     * @param rMeshConstant: The constant that appears in papers an none knows where comes from, where...
     */
    
    void ComputeHessianMetric(
        ModelPart& rThisModelPart,
        const double& rMaxSize,
        Variable<double> & rVariable,
        const double rInterpError = 1.0e-6,
        const double rMeshConstant = 0.28125 // TODO: This is for tetrahedron, look in literature for more
        )
    {
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
        CalculateAuxiliarHessian(rThisModelPart, rVariable);
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            if ( itNode->SolutionStepsDataHas( rVariable ) == false )
            {
                KRATOS_THROW_ERROR( std::invalid_argument, "Missing variable on node ", itNode->Id() )
            }
            
            const double distance = itNode->FastGetSolutionStepValue(DISTANCE, 0);
            const Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);
            
            double element_size = itNode->FastGetSolutionStepValue(NODAL_H, 0);
            if (element_size > mMinSize)
            {
                element_size = mMinSize;
            }
            
            const double ratio = CalculateAnisotropicRatio(distance);
            
            // For postprocess pourposes
            double& anisotropic_ratio = itNode->FastGetSolutionStepValue(ANISOTROPIC_RATIO, 0); 
            anisotropic_ratio = ratio;
            
            // We compute the metric
            #ifdef KRATOS_DEBUG 
            if( itNode->Has(MMG_METRIC) == false) 
            {
                KRATOS_ERROR <<  " MMG_METRIC not defined for node " << itNode>Id();
            }
            #endif     
            Vector& metric = itNode->GetValue(MMG_METRIC);
            
            #ifdef KRATOS_DEBUG 
            if(metric.size() != TDim * 3 - 3) 
            {
                KRATOS_ERROR << "Wrong size of vector MMG_METRIC found for node " << itNode>Id() << " size is " << metric.size() << " expected size was " << TDim * 3 - 3;
            }
            #endif
            
            const double normmetric = norm_2(metric);
            if (normmetric > 0.0) // NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            {
                const Vector old_metric = itNode->GetValue(MMG_METRIC);
                const Vector new_metric = ComputeHessianMetricTensor(hessian, element_size, rMaxSize, rInterpError, rMeshConstant, ratio);
//                 const Vector new_metric = ComputeHessianMetricTensor(hessian, mMinSize, rMaxSize, rInterpError, rMeshConstant, ratio);        
                
                metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
            }
            else
            {
                metric = ComputeHessianMetricTensor(hessian, element_size, rMaxSize, rInterpError, rMeshConstant, ratio);
//                 metric = ComputeHessianMetricTensor(hessian, mMinSize, rMaxSize, rInterpError, rMeshConstant, ratio);        
            }
        }
    }
       
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{
    
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    double mMinSize;            // The minimal size of the elements
    double mAnisRatio;          // The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;         // The boundary layer limit distance
    std::string mInterpolation; // The interpolation type
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * It calculates the tensor of the scalar, necessary to get the solution before remeshing
     * @param gradient_value: The gradient of the scalar to remesh
     * @param ratio: The alpha parameter used to remesh
     * @param element_size: The minimum size of the elements
     * @param node_id: The id of the node
     */
        
    Vector ComputeLevelSetMetricTensor(
        const array_1d<double, 3>& gradient_value,
        const double& ratio,
        const double& element_size
    )
    {
        Vector metric;
        metric.resize(3 * (1 + TDim), false);
        
        const double coeff0 = 1.0/(element_size * element_size);
        
        if (TDim == 2) // 2D: The order of the metric is m11,m12,m22
        {
            const double coeff1 = coeff0/(ratio * ratio);
            
            const double v0v0 = gradient_value[0]*gradient_value[0];
            const double v0v1 = gradient_value[0]*gradient_value[1];
            const double v1v1 = gradient_value[1]*gradient_value[1];
            
            metric[0] = coeff0*(1.0 - v0v0) + coeff1*v0v0;
            metric[1] = coeff0*(    - v0v1) + coeff1*v0v1;  
            metric[2] = coeff0*(1.0 - v1v1) + coeff1*v1v1;

        }
        else // 3D: The order of the metric is m11,m12,m13,m22,m23,m33
        {
            const double coeff1 = coeff0/(ratio * ratio);
            
            const double v0v0 = gradient_value[0]*gradient_value[0];
            const double v0v1 = gradient_value[0]*gradient_value[1];
            const double v0v2 = gradient_value[0]*gradient_value[2];
            const double v1v1 = gradient_value[1]*gradient_value[1];
            const double v1v2 = gradient_value[1]*gradient_value[2];
            const double v2v2 = gradient_value[2]*gradient_value[2];
            
            metric[0] = coeff0*(1.0 - v0v0) + coeff1*v0v0;
            metric[1] = coeff0*(    - v0v1) + coeff1*v0v1; 
            metric[2] = coeff0*(    - v0v2) + coeff1*v0v2; 
            metric[3] = coeff0*(1.0 - v1v1) + coeff1*v1v1; 
            metric[4] = coeff0*(    - v1v2) + coeff1*v1v2; 
            metric[5] = coeff0*(1.0 - v2v2) + coeff1*v2v2;
        }
        
        return metric;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * This function is used to compute the Hessian metric tensor, note that when using the Hessian, more than one metric can be defined simultaneously, so in consecuence we need to define the elipsoid which defines the volume of maximal intersection
     * @param hessian: The hessian tensor condensed already computed
     * @param mMinSize: The minimum size of the elements
     * @param rMaxSize: The maximal size of the elements
     * @param rInterpError: The interpolation error assumed
     * @param rMeshConstant: The constant that appears in papers an none knows where comes from, where...
     */
        
    Vector ComputeHessianMetricTensor(
        const Vector& hessian,
        const double& mMinSize,
        const double& rMaxSize,
        const double& rInterpError,
        const double& rMeshConstant,
        const double& ratio
        )
    {        
        // Calculating metric parameters
        const double c_epsilon = rMeshConstant/rInterpError;
        const double min_ratio = 1.0/(mMinSize * mMinSize);
        const double max_ratio = 1.0/(rMaxSize * rMaxSize);
        
        typedef boost::numeric::ublas::bounded_matrix<double, TDim, TDim> temp_type;
        
        // Declaring the eigen system
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> eigen_vector_matrix;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> eigen_values_matrix;

        // We first transform into a matrix
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> hessian_matrix = MetricsMathUtils<TDim>::VectorToTensor(hessian);
        
        MetricsMathUtils<TDim>::EigenSystem(hessian_matrix, eigen_vector_matrix, eigen_values_matrix, 1e-18, 10);
            
        // Recalculate the metric eigen values
        for (unsigned int i = 0; i < TDim; i++)
        {
            eigen_values_matrix(i, i) = MathUtils<double>::Min(MathUtils<double>::Max(c_epsilon * std::abs(eigen_values_matrix(i, i)), max_ratio), min_ratio);
        }
            
        // Considering anisotropic
        if (ratio < 1.0)
        {
            double eigen_max = eigen_values_matrix(0, 0);
            double eigen_min = eigen_values_matrix(1, 1);
            for (unsigned int i = 1; i < TDim - 1; i++)
            {
                eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
                eigen_min = MathUtils<double>::Min(eigen_max, eigen_values_matrix(i, i));
            }
            const double eigen_radius = std::abs(eigen_max - eigen_min) * (1.0 - ratio);
            const double rel_eigen_radius = std::abs(eigen_max - eigen_radius);
            
            for (unsigned int i = 0; i < TDim; i++)
            {
                eigen_values_matrix(i, i) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(i, i), eigen_max), rel_eigen_radius);
            }
        }
            
        // We compute the product
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> metric_matrix =  prod(trans(eigen_vector_matrix), prod<temp_type>(eigen_values_matrix,eigen_vector_matrix));
        
        // Finally we transform to a vector
        const Vector metric = MetricsMathUtils<TDim>::TensorToVector(metric_matrix);
        
        return metric;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This calculates the auxiliar hessian needed for the metric
     * @param rThisModelPart: The original model part where we compute the hessian
     * @param rVariable: The variable to calculate the hessian
     */
    
    void CalculateAuxiliarHessian(
        ModelPart& rThisModelPart,
        Variable<double> & rVariable
        )
    {
        // Iterate in the nodes
        NodesArrayType& pNode = rThisModelPart.Nodes();
        auto numNodes = pNode.end() - pNode.begin();
        
//         #pragma omp parallel for // NOTE: Be careful with the parallel (MUST BE INITIALIZED TO BE THREAD SAFE)
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            
            Vector& hessian = itNode->GetValue(AUXILIAR_HESSIAN);  
            hessian = ZeroVector(3 * (TDim - 1));
        }
        
        // Compute auxiliar gradient
        ComputeNodalGradientProcess<TDim> GradientProcess = ComputeNodalGradientProcess<TDim>(rThisModelPart, rVariable, AUXILIAR_GRADIENT, NODAL_AREA);
        GradientProcess.Execute();
        
        // Iterate in the conditions
        ElementsArrayType& pElement = rThisModelPart.Elements();
        auto numElements = pElement.end() - pElement.begin();
        
        #pragma omp parallel for
        for(unsigned int i = 0; i < numElements; i++) 
        {
            auto itElem = pElement.begin() + i;
            
            Element::GeometryType& geom = itElem->GetGeometry();
            const unsigned int geom_size = geom.size();

            double Volume;
            if (TDim == 2)
            {
                if (geom_size == 3)
                {
                    boost::numeric::ublas::bounded_matrix<double,3, 2> DN_DX;
                    array_1d<double, 3> N;
       
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                    
                    boost::numeric::ublas::bounded_matrix<double,3, 2> values;
                    for(unsigned int i_node = 0; i_node < 3; i_node++)
                    {
                        const array_1d<double, 3> aux_grad = geom[i_node].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                        values(i_node, 0) = aux_grad[0];
                        values(i_node, 1) = aux_grad[1];
                    }
                    
                    const boost::numeric::ublas::bounded_matrix<double,2, 2> Hessian = prod(trans(DN_DX), values); 
                    const Vector HessianCond = MetricsMathUtils<TDim>::TensorToVector(Hessian);
                    
                    for(unsigned int i_node = 0; i_node < geom_size; i_node++)
                    {
                        for(unsigned int k = 0; k < 3; k++)
                        {
                            double& val = geom[i_node].FastGetSolutionStepValue(AUXILIAR_HESSIAN)[k];
                            
                            #pragma omp atomic
                            val += N[i_node] * Volume * HessianCond[k];
                        }
                    }
                }
                else
                {
                    KRATOS_THROW_ERROR( std::logic_error, "WARNING: YOU CAN USE JUST TRIANGLES RIGHT NOW IN THE GEOMETRY UTILS: ", geom_size );
                }
            }
            else
            {
                if (geom_size == 4)
                {
                    boost::numeric::ublas::bounded_matrix<double,4,  3> DN_DX;
                    array_1d<double, 4> N;
                    
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);
                    
                    boost::numeric::ublas::bounded_matrix<double,4, 3> values;
                    for(unsigned int i_node = 0; i_node < 4; i_node++)
                    {
                        const array_1d<double, 3> aux_grad = geom[i_node].FastGetSolutionStepValue(AUXILIAR_GRADIENT);
                        values(i_node, 0) = aux_grad[0];
                        values(i_node, 1) = aux_grad[1];
                        values(i_node, 2) = aux_grad[2];
                    }
                    
                    const boost::numeric::ublas::bounded_matrix<double, 3, 3> Hessian = prod(trans(DN_DX), values); 
                    const Vector HessianCond = MetricsMathUtils<TDim>::TensorToVector(Hessian);
                    
                    for(unsigned int i_node = 0; i_node < geom_size; i_node++)
                    {
                        for(unsigned int k = 0; k < 6; k++)
                        {
                            double& val = geom[i_node].FastGetSolutionStepValue(AUXILIAR_HESSIAN)[k];
                            
                            #pragma omp atomic
                            val += N[i_node] * Volume * HessianCond[k];
                        }
                    }
                }
                else
                {
                    KRATOS_THROW_ERROR( std::logic_error, "WARNING: YOU CAN USE JUST TETRAEDRA RIGHT NOW IN THE GEOMETRY UTILS: ", geom.size() );
                }
            }
        }
            
        #pragma omp parallel for
        for(unsigned int i = 0; i < numNodes; i++) 
        {
            auto itNode = pNode.begin() + i;
            itNode->GetValue(AUXILIAR_HESSIAN) /= itNode->FastGetSolutionStepValue(NODAL_AREA);
            
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * This calculates the anisotropic ratio
     * @param distance: Distance parameter
     */
    
    double CalculateAnisotropicRatio(const double& distance)
    {
        const double tolerance = 1.0e-12;
        double ratio = 1.0; // NOTE: Isotropic mesh
        if (mAnisRatio < 1.0)
        {                           
            if (std::abs(distance) <= mBoundLayer)
            {
                if (mInterpolation.find("Constant") != std::string::npos)
                {
                    ratio = mAnisRatio;
                }
                else if (mInterpolation.find("Linear") != std::string::npos)
                {
                    ratio = mAnisRatio + (std::abs(distance)/mBoundLayer) * (1.0 - mAnisRatio);
                }
                else if (mInterpolation.find("Exponential") != std::string::npos)
                {
                    ratio = - std::log(std::abs(distance)/mBoundLayer) * mAnisRatio + tolerance;
                    if (ratio > 1.0)
                    {
                        ratio = 1.0;
                    }
                }
                else
                {
                    std::cout << "No interpolation defined, considering linear" << std:: endl;
                    ratio = 
                    mAnisRatio + (std::abs(distance)/mBoundLayer) * (1.0 - mAnisRatio);
                }
            }
        }
        
        return ratio;
    }
    
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
};// class MetricsUtility
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}// namespace Kratos.
#endif /* KRATOS_METRICS_UTILITY defined */
