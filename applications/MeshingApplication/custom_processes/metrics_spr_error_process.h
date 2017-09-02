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

#if !defined(KRATOS_SPR_ERROR_METRICS_PROCESS)
#define KRATOS_SPR_ERROR_METRICS_PROCESS

// Project includes
#include "utilities/math_utils.h"
#include "custom_utilities/metrics_math_utils.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "meshing_application.h"
#include "processes/compute_nodal_gradient_process.h" // TODO: Not prism or quadrilaterals implemented yet
#include "processes/find_nodal_neighbours_process.h"


namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef ModelPart::NodesContainerType                                     NodesArrayType;
    typedef ModelPart::ElementsContainerType                               ElementsArrayType;
    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;
    typedef Node <3>                                                                NodeType;
    
///@}
///@name  Enum's
///@{
    
    #if !defined(INTERPOLATION_METRIC)
    #define INTERPOLATION_METRIC
        enum Interpolation {Constant = 0, Linear = 1, Exponential = 2};
    #endif
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

//// This class is can be used to compute the metrics of the model part with an Hessian approach

//template<unsigned int TDim, class TVarType>  
template<unsigned int TDim> 
class ComputeSPRErrorSolMetricProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ComputeSPRErrorSolMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeSPRErrorSolMetricProcess);
    
    ///@}
    ///@name Life Cycle
    ///@{
     
    // Constructor
    
    /**
     * This is the default constructor
     * @param rThisModelPart: The model part to be computed
     * @param ThisParameters: The input parameters
     */
    
    ComputeSPRErrorSolMetricProcess(
        ModelPart& rThisModelPart,
        //TVarType& rVariable,
        Parameters ThisParameters = Parameters(R"({})")
        )
        :mThisModelPart(rThisModelPart)
        //:mThisModelPart(rThisModelPart),
        //mVariable(rVariable)
    {               
      /*  Parameters DefaultParameters = Parameters(R"(
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
        })" );*/
        Parameters DefaultParameters = Parameters(R"(
        {
            "minimal_size"                        : 0.1,
            "maximal_size"                        : 10.0, 
            "error"                               : 0.05
        })" );
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
         
        mMinSize = ThisParameters["minimal_size"].GetDouble();
        mMaxSize = ThisParameters["maximal_size"].GetDouble();
        /*mEnforceCurrent = ThisParameters["enforce_current"].GetBool();
        
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
        }*/
    }
    
    /// Destructor.
    virtual ~ComputeSPRErrorSolMetricProcess() {}
    
    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * We initialize the metrics of the MMG sol using the Hessian metric matrix approach
     */
    
    virtual void Execute()
    {
        CalculateSuperconvergentPatchRecovery();

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
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ComputeSPRErrorSolMetricProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ComputeSPRErrorSolMetricProcess";
    }

    /// Print object"s data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
    
protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


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
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{
    
    ModelPart& mThisModelPart;               // The model part to compute
    //TVarType mVariable;            // The variable to calculate the hessian
    double mMinSize;                         // The minimal size of the elements
    double mMaxSize;                         // The maximal size of the elements
    bool mEnforceCurrent;                    // With this we choose if we inforce the current nodal size (NODAL_H)
    double mInterpError;                     // The error of interpolation allowed
    double mMeshConstant;                    // The mesh constant to remesh (depends of the element type)
    double mAnisRatio;                       // The minimal anisotropic ratio (0 < ratio < 1)
    double mBoundLayer;                      // The boundary layer limit distance
    Interpolation mInterpolation;            // The interpolation type
    

    void CalculateSuperconvergentPatchRecovery()
    {
        /************************************************************************
        --1-- calculate superconvergent stresses (at the nodes) --1--
        ************************************************************************/
        
        //std::vector<std::string> submodels;
        //submodels= mThisModelPart.GetSubModelPartNames();
        //for (std::vector<std::string>::const_iterator i = submodels.begin();i!=submodels.end();i++) 
        //    std::cout << *i<<std::endl; 
        FindNodalNeighboursProcess findNeighbours(mThisModelPart);
        findNeighbours.Execute();
        //std::vector<Vector> stress_vector(1);
        //std::vector<array_1d<double,3>> coordinates_vector(1);
        //ariable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
        //iteration over all nodes -- construction of patches
        ModelPart::NodesContainerType& rNodes = mThisModelPart.Nodes();
        for(ModelPart::NodesContainerType::iterator i_nodes = rNodes.begin(); i_nodes!=rNodes.end(); i_nodes++){
            int neighbour_size = i_nodes->GetValue(NEIGHBOUR_ELEMENTS).size();
            std::cout << "Node: " << i_nodes->Id() << " has " << neighbour_size << " neighbouring elements: " << std::endl;
            Vector sigma_recovered(3,0);
            if(neighbour_size>2){ 
                CalculatePatch(i_nodes,i_nodes,neighbour_size,sigma_recovered);
                i_nodes->SetValue(RECOVERED_STRESS,sigma_recovered);
                std::cout<<"recovered sigma"<<sigma_recovered<<std::endl;
            }
            else{
                for(WeakPointerVector< Node<3> >::iterator i_neighbour_nodes = i_nodes->GetValue(NEIGHBOUR_NODES).begin(); i_neighbour_nodes != i_nodes->GetValue(NEIGHBOUR_NODES).end(); i_neighbour_nodes++){
                    Vector sigma_recovered_i(3);
                    unsigned int count_i=0;
                    for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++){
                        if (i->Id() == i_neighbour_nodes->Id() && i->GetValue(NEIGHBOUR_ELEMENTS).size()>2){
                            CalculatePatch(i_nodes,i,neighbour_size,sigma_recovered_i);
                            count_i ++;
                        }
                    }
                    //average solution from different patches
                    if(count_i != 0)
                        sigma_recovered =sigma_recovered*(count_i-1)/count_i + sigma_recovered_i/count_i;
                }
                i_nodes->SetValue(RECOVERED_STRESS,sigma_recovered);
            }
       }
        /******************************************************************************
        --2-- calculate error estimation and new element size (for each element) --2--
        ******************************************************************************/
        //loop over all elements: 
        double error_overall_squared=0;
        double energy_norm_overall_squared=0;

        //compute the error estimate per element
        for(ModelPart::ElementsContainerType::iterator i_elements = mThisModelPart.Elements().begin() ; i_elements != mThisModelPart.Elements().end(); i_elements++) 
        {
            std::vector<double> error_integration_point;
            i_elements->GetValueOnIntegrationPoints(ERROR_INTEGRATION_POINT,error_integration_point,mThisModelPart.GetProcessInfo());
            double error_energy_norm=0;
            for(unsigned int i=0;i<error_integration_point.size();i++)
                error_energy_norm += error_integration_point[i];
            error_overall_squared += error_energy_norm;
            error_energy_norm= sqrt(error_energy_norm);
            i_elements->SetValue(ELEMENT_ERROR,error_energy_norm);
            std::cout<<"element_error:"<<error_energy_norm<<std::endl;


            std::vector<double> strain_energy;
            i_elements->GetValueOnIntegrationPoints(STRAIN_ENERGY,strain_energy,mThisModelPart.GetProcessInfo());
            double energy_norm=0;
            for(unsigned int i=0;i<strain_energy.size();i++)
                energy_norm += 2*strain_energy[i];
            energy_norm_overall_squared += energy_norm;
            energy_norm= sqrt(energy_norm);
            std::cout<<"energy norm:"<<energy_norm<<std::endl;
        }
        std::cout<<"overall error norm (squared):"<<error_overall_squared<<std::endl;
        std::cout<<"overall energy norm (squarde):"<<energy_norm_overall_squared<<std::endl;
        
        //compute new element size
        for(ModelPart::ElementsContainerType::iterator i_elements = mThisModelPart.Elements().begin() ; i_elements != mThisModelPart.Elements().end(); i_elements++) 
        {
            //compute the current element size h
            i_elements->CalculateElementSize();

            //compute new element size
            double new_element_size;
            new_element_size = i_elements->GetValue(ELEMENT_H)/i_elements->GetValue(ELEMENT_ERROR);
            new_element_size *= sqrt((energy_norm_overall_squared+error_overall_squared)/mThisModelPart.Elements().size())*0.05;
            std::cout<<"old element size: "<<i_elements->GetValue(ELEMENT_H)<<std::endl;
            i_elements->SetValue(ELEMENT_H,new_element_size);
            std::cout<<"new element size: "<<i_elements->GetValue(ELEMENT_H)<<std::endl;
        }

        /******************************************************************************
        --3-- calculate metric (for each node) --3--
        ******************************************************************************/

        for(ModelPart::NodesContainerType::iterator i_nodes = rNodes.begin(); i_nodes!=rNodes.end(); i_nodes++){
            // get maximal element size from neighboring elements
            double h_min=0;
            for(WeakPointerVector< Element >::iterator i_neighbour_elements = i_nodes->GetValue(NEIGHBOUR_ELEMENTS).begin(); i_neighbour_elements != i_nodes->GetValue(NEIGHBOUR_ELEMENTS).end(); i_neighbour_elements++){
                if(h_min==0||h_min>i_neighbour_elements->GetValue(ELEMENT_H))
                    h_min = i_neighbour_elements->GetValue(ELEMENT_H);
                
            }
            //std::cout<<"h_min: "<<h_min<<std::endl;
            // set metric
            Matrix metric_matrix(2,2,0);
            metric_matrix(0,0)=1/(h_min*h_min);
            metric_matrix(1,1)=1/(h_min*h_min);
            // transform metric matrix to a vector
            const Vector metric = MetricsMathUtils<TDim>::TensorToVector(metric_matrix);
            i_nodes->SetValue(MMG_METRIC,metric);


            std::cout<<"metric: "<<i_nodes->GetValue(MMG_METRIC)<<std::endl;
        }
    }
    //calculates the recovered stress at a node 
    // i_node: the node for which the recovered stress should be calculated
    // i_patch_node: the center node of the patch
    void CalculatePatch(
        ModelPart::NodesContainerType::iterator i_nodes,
        ModelPart::NodesContainerType::iterator i_patch_node,
        int neighbour_size,
        Vector& rsigma_recovered)
    {
        std::vector<Vector> stress_vector(1);
        std::vector<array_1d<double,3>> coordinates_vector(1);
        Variable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
        Variable<Vector> variable_stress = CAUCHY_STRESS_VECTOR;
        Matrix A(3,3,0);
        Matrix b(3,3,0); 
        Matrix p_k(1,3,0);
        for( WeakPointerVector< Element >::iterator i_elements = i_patch_node->GetValue(NEIGHBOUR_ELEMENTS).begin(); i_elements != i_patch_node->GetValue(NEIGHBOUR_ELEMENTS).end(); i_elements++) {
            std::cout << "\tElement: " << i_elements->Id() << std::endl;
            i_elements->GetValueOnIntegrationPoints(variable_stress,stress_vector,mThisModelPart.GetProcessInfo());
            i_elements->GetValueOnIntegrationPoints(variable_coordinates,coordinates_vector,mThisModelPart.GetProcessInfo());

            std::cout << "\tstress: " << stress_vector[0] << std::endl;
            std::cout << "\tx: " << coordinates_vector[0][0] << "\ty: " << coordinates_vector[0][1] << "\tz_coordinate: " << coordinates_vector[0][2] << std::endl;
            Matrix sigma(1,3);
            for(int j=0;j<3;j++)
                sigma(0,j)=stress_vector[0][j];
            p_k(0,0)=1;
            p_k(0,1)=coordinates_vector[0][0]-i_patch_node->X(); 
            p_k(0,2)=coordinates_vector[0][1]-i_patch_node->Y();   
            A+=prod(trans(p_k),p_k);
            b+=prod(trans(p_k),sigma);
        }
        Matrix invA(3,3);
        double det;
        MathUtils<double>::InvertMatrix(A,invA,det);
        //std::cout <<A<<std::endl;
        //std::cout <<invA<<std::endl;
        //std::cout << det<< std::endl;

        Matrix coeff(3,3);
        coeff = prod(invA,b);
        if(neighbour_size > 2)
            rsigma_recovered = MatrixRow(coeff,0);
        else{
            p_k(0,1)=i_nodes->X()-i_patch_node->X(); 
            p_k(0,2)=i_nodes->Y()-i_patch_node->Y();
            Matrix sigma(1,3);
            sigma = prod(p_k,coeff);
            rsigma_recovered = MatrixRow(sigma,0);
        }
    }

    /// Assignment operator.
    ComputeSPRErrorSolMetricProcess& operator=(ComputeSPRErrorSolMetricProcess const& rOther);

    /// Copy constructor.
    //ComputeSPRErrorSolMetricProcess(ComputeSPRErrorSolMetricProcess const& rOther);

};// class ComputeSPRErrorSolMetricProcess

};// namespace Kratos.
#endif /* KRATOS_SPR_ERROR_METRICS_PROCESS defined */
