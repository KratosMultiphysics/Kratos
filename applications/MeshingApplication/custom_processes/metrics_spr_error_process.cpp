// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Anna Rehr
//

// Project includes
#include "custom_processes/metrics_spr_error_process.h"

namespace Kratos
{
template<unsigned int TDim> 
SPRMetricProcess<TDim>::SPRMetricProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    )
    :mThisModelPart(rThisModelPart)
{               
    Parameters DefaultParameters = Parameters(R"(
    {
        "minimal_size"                        : 0.01,
        "maximal_size"                        : 10.0, 
        "error"                               : 0.1,
        "penalty_normal"                      : 10000.0,
        "penalty_tangential"                  : 10000.0,
        "echo_level"                          : 0,
        "set_number_of_elements"              : false,
        "number_of_elements"                  : 1000,
        "average_nodal_h"                     : false
    })" 
    );
    
    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
    
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mPenaltyNormal = ThisParameters["penalty_normal"].GetDouble();
    mPenaltyTangent = ThisParameters["penalty_tangential"].GetDouble();
    mEchoLevel = ThisParameters["echo_level"].GetInt();
    mSigmaSize = (TDim == 2) ? 3 : 6;
    mSetElementNumber = ThisParameters["set_number_of_elements"].GetBool();
    mElementNumber = ThisParameters["number_of_elements"].GetInt();
    mTargetError = ThisParameters["error"].GetDouble();
    mAverageNodalH = ThisParameters["average_nodal_h"].GetBool();
}
    
/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim> 
void SPRMetricProcess<TDim>::Execute()
{

    /************************************************************************
    --1-- calculate superconvergent stresses (at the nodes) --1--
    ************************************************************************/

    FindNodalNeighboursProcess find_neighbours(mThisModelPart);
    find_neighbours.Execute();

    // Iteration over all nodes -- construction of patches
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    int num_nodes = nodes_array.end() - nodes_array.begin();
    
    for(int i_node = 0; i_node < num_nodes; ++i_node) {
        auto it_node = nodes_array.begin() + i_node;
        
        const unsigned int neighbour_size = it_node->GetValue(NEIGHBOUR_ELEMENTS).size();

        Vector sigma_recovered(mSigmaSize, 0);
        
        if(neighbour_size > TDim){ 
            
            CalculatePatch(it_node, it_node, neighbour_size,sigma_recovered);
            it_node->SetValue(RECOVERED_STRESS,sigma_recovered);
            
            if(mEchoLevel>2)
                std::cout<<"recovered sigma"<<sigma_recovered<<std::endl;
        }
        else{
            auto& neigh_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            for(auto it_neighbour_nodes = neigh_nodes.begin(); it_neighbour_nodes != neigh_nodes.end(); it_neighbour_nodes++){
                
                Vector sigma_recovered_i(mSigmaSize,0);
                
                unsigned int count_i=0;
                for(int i_node_loop = 0; i_node_loop < num_nodes; ++i_node_loop) { // FIXME: Avoid thsi double loop, extreamily expensive
                    auto it_node_loop = nodes_array.begin() + i_node_loop;
                    const std::size_t size_elem_neigh = it_node_loop->GetValue(NEIGHBOUR_ELEMENTS).size();
                    if (it_node_loop->Id() == it_neighbour_nodes->Id() && size_elem_neigh > TDim){
                        CalculatePatch(it_node, it_node_loop, neighbour_size, sigma_recovered_i);
                        ++count_i;
                    }
                }
                
                // Average solution from different patches
                if(count_i != 0)
                    sigma_recovered = sigma_recovered*(count_i-1)/count_i + sigma_recovered_i/count_i;
            }
            
            it_node->SetValue(RECOVERED_STRESS,sigma_recovered);
            if(mEchoLevel > 2)
                std::cout << "Recovered sigma: " << sigma_recovered << std::endl;
        }
    }
    /******************************************************************************
    --2-- calculate error estimation and new element size (for each element) --2--
    ******************************************************************************/
    //loop over all elements: 
    double error_overall = 0.0;
    double energy_norm_overall = 0.0;

    ElementsArrayType& elements_array = mThisModelPart.Elements();
    int num_elem = elements_array.end() - elements_array.begin();
    
    // Compute the error estimate per element
    for(int i_elem = 0; i_elem < num_elem; ++i_elem){
        auto it_elem = elements_array.begin() + i_elem;
        
        std::vector<double> error_integration_point;
        auto& process_info = mThisModelPart.GetProcessInfo();
        it_elem->GetValueOnIntegrationPoints(ERROR_INTEGRATION_POINT, error_integration_point, process_info);
        //std::cout<<"Error:"<<error_integration_point<<std::endl;
        double error_energy_norm = 0.0;
        for(unsigned int i = 0;i < error_integration_point.size();++i)
            error_energy_norm += error_integration_point[i];
        error_overall += error_energy_norm;
        error_energy_norm = std::sqrt(error_energy_norm);
        it_elem->SetValue(ELEMENT_ERROR, error_energy_norm);
        
        if (mEchoLevel > 2)
            std::cout << "Element_error: " << error_energy_norm << std::endl;

        std::vector<double> strain_energy;
        it_elem->GetValueOnIntegrationPoints(STRAIN_ENERGY, strain_energy, process_info);
        
        double energy_norm = 0.0;
        for(unsigned int i = 0;i < strain_energy.size(); ++i)
            energy_norm += 2.0 * strain_energy[i];
        energy_norm_overall += energy_norm;
        energy_norm= std::sqrt(energy_norm);
        
        if(mEchoLevel > 2)
            std::cout << "Energy norm: " << energy_norm << std::endl;
    }
    
    error_overall = std::sqrt(error_overall);
    energy_norm_overall = std::sqrt(energy_norm_overall);
    double error_percentage = error_overall/std::sqrt((std::pow(error_overall, 2) + std::pow(energy_norm_overall, 2)));
    
    if(mEchoLevel>1){
        std::cout << "Overall error norm: " << error_overall << std::endl;
        std::cout << "Overall energy norm: "<< energy_norm_overall << std::endl;
        std::cout << "Error in percent: " << error_percentage << std::endl;}
    
    // Compute new element size
    for(int i_elem = 0; i_elem < num_elem; ++i_elem){
        auto it_elem = elements_array.begin() + i_elem;
        
        //Compute the current element size h
        //it_elem->CalculateElementSize();
        ComputeElementSize(it_elem);

        // Compute new element size
        double new_element_size;
        new_element_size = it_elem->GetValue(ELEMENT_H)/it_elem->GetValue(ELEMENT_ERROR);

        // if a target number for elements is given: use this, else: use current element number
        //if(mSetElementNumber == true && mElementNumber<mThisModelPart.Elements().size())
        if(mSetElementNumber == true)
        new_element_size *= std::sqrt((std::pow(energy_norm_overall, 2)+ std::pow(error_overall, 2))/mElementNumber) * mTargetError;
        else
        new_element_size *= std::sqrt((energy_norm_overall*energy_norm_overall+error_overall*error_overall)/mThisModelPart.Elements().size())*mTargetError;
        
        
        // Check if element sizes are in specified limits. If not, set them to the limit case
        if(new_element_size < mMinSize)
            new_element_size = mMinSize;
        if(new_element_size > mMaxSize)
            new_element_size = mMaxSize;

        
        it_elem->SetValue(ELEMENT_H, new_element_size);
    }

    /******************************************************************************
    --3-- calculate metric (for each node) --3--
    ******************************************************************************/

    for(int i_node = 0; i_node < num_nodes; ++i_node) {
        auto it_node = nodes_array.begin() + i_node;
        /**************************************************************************
        ** Determine nodal element size h:
        ** if average_nodal_h == true : the nodal element size is averaged from the element size of neighboring elements
        ** if average_nodal_h == false: the nodal element size is the minimum element size from neighboring elements
        */
        double h_min = 0.0;
        auto& neigh_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
        for(WeakElementItType i_neighbour_elements = neigh_elements.begin(); i_neighbour_elements != neigh_elements.end(); i_neighbour_elements++){
            const double element_h = i_neighbour_elements->GetValue(ELEMENT_H);
            if(mAverageNodalH == false){
                if(h_min == 0.0 || h_min > element_h) h_min = element_h;
            }
            else
                h_min += element_h;
        }
        if(mAverageNodalH == true)
        h_min = h_min/double(neigh_elements.size());

        // Set metric
        Matrix metric_matrix(TDim, TDim, 0.0);
        for(unsigned int i = 0;i < TDim; ++i)
            metric_matrix(i,i) = 1.0/std::pow(h_min, 2);

        // Transform metric matrix to a vector
        const Vector metric = MetricsMathUtils<TDim>::TensorToVector(metric_matrix);
        it_node->SetValue(MMG_METRIC, metric);

        if(mEchoLevel>2)
            std::cout<<"Node "<<it_node->Id()<<" has metric: "<<it_node->GetValue(MMG_METRIC)<<std::endl;
    }
    
    mThisModelPart.GetProcessInfo()[ERROR_ESTIMATE] = error_overall/std::pow((error_overall*error_overall+energy_norm_overall*energy_norm_overall),0.5);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim> 
void SPRMetricProcess<TDim>::CalculatePatch(
    NodeItType itNode,
    NodeItType itPatchNode,
    unsigned int NeighbourSize,
    Vector& rSigmaRecovered)
{
    // Determine if contact BC has to be regarded
    const bool regard_contact = (itNode->GetValue(CONTACT_PRESSURE) != 0.0);
    
    //regard_contact = itPatchNode->Has(CONTACT_PRESSURE);
    /*if(regard_contact == false)
    {
        for( auto it_neighbour_nodes = itPatchNode->GetValue(NEIGHBOUR_NODES).begin(); it_neighbour_nodes != itPatchNode->GetValue(NEIGHBOUR_NODES).end(); it_neighbour_nodes++) {
            if (it_neighbour_nodes->Has(CONTACT_PRESSURE))
            {
                regard_contact = true;
                break;
            }
        }
    }*/
    
    if (regard_contact == false)
        CalculatePatchStandard(itNode, itPatchNode, NeighbourSize, rSigmaRecovered);
    else
        CalculatePatchContact(itNode, itPatchNode, NeighbourSize, rSigmaRecovered);
}

/***********************************************************************************/
/***********************************************************************************/
    
template<unsigned int TDim> 
void SPRMetricProcess<TDim>::CalculatePatchStandard(
    NodeItType itNode,
    NodeItType itPatchNode,
    unsigned int NeighbourSize,
    Vector& rSigmaRecovered
    )
{
    std::vector<Vector> stress_vector(1);
    std::vector<array_1d<double,3>> coordinates_vector(1);
    Variable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
    Variable<Vector> variable_stress = CAUCHY_STRESS_VECTOR;
    Matrix A(TDim+1,TDim+1,0);
    Matrix b(TDim+1,mSigmaSize,0); 
    Matrix p_k(1,TDim+1,0);
    
    auto& neigh_elements = itPatchNode->GetValue(NEIGHBOUR_ELEMENTS);
    for( WeakElementItType it_elem = neigh_elements.begin(); it_elem != neigh_elements.end(); ++it_elem) {
        
        it_elem->GetValueOnIntegrationPoints(variable_stress,stress_vector,mThisModelPart.GetProcessInfo());
        it_elem->GetValueOnIntegrationPoints(variable_coordinates,coordinates_vector,mThisModelPart.GetProcessInfo());

        //std::cout << "\tstress: " << stress_vector[0] << std::endl;
        //std::cout << "\tx: " << coordinates_vector[0][0] << "\ty: " << coordinates_vector[0][1] << "\tz_coordinate: " << coordinates_vector[0][2] << std::endl;
        
        Matrix sigma(1,mSigmaSize);
        for(unsigned int j = 0; j < mSigmaSize; ++j)
            sigma(0,j)=stress_vector[0][j];
        p_k(0,0)=1;
        p_k(0,1)=coordinates_vector[0][0]-itPatchNode->X(); 
        p_k(0,2)=coordinates_vector[0][1]-itPatchNode->Y();
        if(TDim == 3)
            p_k(0,3)=coordinates_vector[0][2]-itPatchNode->Z();       
        
        A += prod(trans(p_k), p_k);
        b += prod(trans(p_k), sigma);
    }
    
    Matrix invA(TDim+1,TDim+1);
    double det;
    MathUtils<double>::InvertMatrix(A,invA,det);
    //std::cout <<A<<std::endl;
    //std::cout <<invA<<std::endl;
    //std::cout << det<< std::endl;
    if(det<1e-10){
        //std::cout<<A<<std::endl;
        for( int i=0; i<TDim+1;i++){
            for( int j=0; j<TDim+1; j++)
                A(i,j)+= 0.001;
        }
        MathUtils<double>::InvertMatrix(A,invA,det);
        std::cout <<"det: "<< det<< std::endl;
    }

    Matrix coeff(TDim+1,mSigmaSize);
    coeff = prod(invA,b);
    
    if(NeighbourSize > TDim)
        rSigmaRecovered = MatrixRow(coeff,0);
    else{
        p_k(0,1)=itNode->X()-itPatchNode->X(); 
        p_k(0,2)=itNode->Y()-itPatchNode->Y();
        if(TDim ==3)
            p_k(0,3)=itNode->Z()-itPatchNode->Z();
        Matrix sigma(1,mSigmaSize);
        sigma = prod(p_k,coeff);
        rSigmaRecovered = MatrixRow(sigma,0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim> 
void SPRMetricProcess<TDim>::CalculatePatchContact(
    NodeItType itNode,
    NodeItType itPatchNode,
    unsigned int NeighbourSize,
    Vector& rSigmaRecovered)
{

    std::vector<Vector> stress_vector(1);
    std::vector<array_1d<double,3>> coordinates_vector(1);
    Variable<array_1d<double,3>> variable_coordinates = INTEGRATION_COORDINATES;
    Variable<Vector> variable_stress = CAUCHY_STRESS_VECTOR;

    CompressedMatrix A((mSigmaSize*(TDim+1)),(mSigmaSize*(TDim+1)),0);
    Matrix b((mSigmaSize*(TDim+1)),1,0); 
    Matrix p_k(mSigmaSize,(mSigmaSize*(TDim+1)),0);
    Matrix N_k(1,mSigmaSize,0);
    Matrix T_k(1,mSigmaSize,0);
    Matrix T_k2(1,mSigmaSize,0);  // in case of 3D: second tangential vector
    Matrix sigma(mSigmaSize,1);
    
    /* Computation A and b */
    // PART 1: contributions from the neighboring elements
    auto& neigh_elements = itPatchNode->GetValue(NEIGHBOUR_ELEMENTS);
    for( WeakElementItType it_elem = neigh_elements.begin(); it_elem != neigh_elements.end(); ++it_elem) {
        
        //std::cout << "\tElement: " << it_elem->Id() << std::endl;
        auto& process_info = mThisModelPart.GetProcessInfo();
        it_elem->GetValueOnIntegrationPoints(variable_stress,stress_vector, process_info);
        it_elem->GetValueOnIntegrationPoints(variable_coordinates,coordinates_vector, process_info);

        //std::cout << "\tstress: " << stress_vector[0] << std::endl;
        //std::cout << "\tx: " << coordinates_vector[0][0] << "\ty: " << coordinates_vector[0][1] << "\tz_coordinate: " << coordinates_vector[0][2] << std::endl;
        
        for( unsigned int j = 0; j < mSigmaSize; ++j)
            sigma(j,0) = stress_vector[0][j];
        
        for ( unsigned int j = 0; j < mSigmaSize; ++j){
            p_k(j,j*(TDim+1))=1;
            p_k(j,j*(TDim+1)+1)=coordinates_vector[0][0]-itPatchNode->X(); 
            p_k(j,j*(TDim+1)+2)=coordinates_vector[0][1]-itPatchNode->Y();
            if(TDim == 3)
                p_k(j,j*(TDim+1)+3)=coordinates_vector[0][2]-itPatchNode->Z();
        }
        
        A += prod(trans(p_k),p_k);
        b += prod(trans(p_k),sigma);
    }
    
    // Computing A and b
    Matrix A1((mSigmaSize * (TDim + 1)),1,0), A2(1,(mSigmaSize *(TDim + 1)),0);
    for (unsigned int j = 0; j < mSigmaSize; ++j){
        p_k(j,j * (TDim + 1) + 1)= itNode->X() - itPatchNode->X();
        p_k(j,j * (TDim + 1) + 2)= itNode->Y() - itPatchNode->Y();
        if(TDim == 3)
            p_k(j,j * (TDim + 1) + 3)= itNode->Z() - itPatchNode->Z();
    }
    
    // Set the normal and tangential vectors in Voigt Notation
    std::vector<double> n(TDim);

    const array_1d<double, 3>& normal = itNode->GetValue(NORMAL);
    for( unsigned int j = 0; j < TDim; ++j) n[j] = normal[j];

    if(TDim ==2){
        N_k(0,0) = n[0]*n[0];
        N_k(0,1) = n[1]*n[1];
        N_k(0,2) = 2*n[0]*n[1];
        
        T_k(0,0) = n[0]*n[1];
        T_k(0,1) = -n[0]*n[1];
        T_k(0,2) = n[1]*n[1]-n[0]*n[0];
    }
    else if (TDim ==3){
        N_k(0,0) = n[0]*n[0];
        N_k(0,1) = n[1]*n[1];
        N_k(0,1) = n[2]*n[2];
        N_k(0,3) = 2*n[1]*n[2];
        N_k(0,4) = 2*n[2]*n[0];
        N_k(0,5) = 2*n[0]*n[1];

        // Set tangential vectors
        std::vector<double> t1(3),t2(3);
        if(n[0]!=0 || n[1] !=0){
            double norm = sqrt((t1[0]*t1[0]+t1[1]*t1[1]));
            t1[0] = n[1]/norm;
            t1[1] = -n[0]/norm;
            t1[2] = 0;

            t2[0] = -n[0]*n[2]/norm;
            t2[1] = -n[1]*n[2]/norm;
            t2[2] = n[0]*n[0]+n[1]*n[1]/norm;
        }
        else{
            t1[0] = 1;
            t1[1] = 0;
            t1[2] = 0;

            t2[0] = 0;
            t2[1] = 1;
            t2[2] = 0;
        }

        T_k(0,0) = n[0]*t1[0];
        T_k(0,1) = n[1]*t1[1];
        T_k(0,2) = n[2]*t1[2];
        T_k(0,3) = n[1]*t1[2]+n[2]*t1[1];
        T_k(0,4) = n[2]*t1[0]+n[0]*t1[2];
        T_k(0,5) = n[0]*t1[1]+n[1]*t1[0];
        
        T_k2(0,0) = n[0]*t2[0];
        T_k2(0,1) = n[1]*t2[1];
        T_k2(0,2) = n[2]*t2[2];
        T_k2(0,3) = n[1]*t2[2]+n[2]*t2[1];
        T_k2(0,4) = n[2]*t2[0]+n[0]*t2[2];
        T_k2(0,5) = n[0]*t2[1]+n[1]*t2[0];
    }
    
    A1 = prod(trans(p_k),trans(N_k));
    A2 = prod(N_k,p_k);
    A += mPenaltyNormal*prod(A1, A2);

    A1 = prod(trans(p_k),trans(T_k));
    A2 = prod(T_k,p_k);
    A += mPenaltyTangent*prod(A1, A2);

    b += mPenaltyNormal*prod(trans(p_k),trans(N_k))*itNode->GetValue(CONTACT_PRESSURE);
    /*
    //PART 2: contributions from contact nodes: regard all nodes from the patch which are in contact
    //patch center node:
    if (itPatchNode->Has(CONTACT_PRESSURE)){
        p_k(0,1)=0;
        p_k(0,2)=0;
        p_k(1,4)=0;
        p_k(1,5)=0;
        p_k(2,7)=0;
        p_k(2,8)=0;
        N_k(0,0) = itPatchNode->GetValue(NORMAL)[0]*itPatchNode->GetValue(NORMAL)[0];
        N_k(0,1) = itPatchNode->GetValue(NORMAL)[1]*itPatchNode->GetValue(NORMAL)[1];
        N_k(0,2) = 2*itPatchNode->GetValue(NORMAL)[0]*itPatchNode->GetValue(NORMAL)[1];
        T_k(0,0) = itPatchNode->GetValue(NORMAL)[0]*itPatchNode->GetValue(NORMAL)[1];
        T_k(0,1) = -itPatchNode->GetValue(NORMAL)[0]*itPatchNode->GetValue(NORMAL)[1];
        T_k(0,2) = itPatchNode->GetValue(NORMAL)[1]*itPatchNode->GetValue(NORMAL)[1]-itPatchNode->GetValue(NORMAL)[0]*itPatchNode->GetValue(NORMAL)[0];

        A1 = prod(trans(p_k),trans(N_k));
        A2 = prod(N_k,p_k);
        A+= mPenaltyNormal*prod(A1, A2);

        A1 = prod(trans(p_k),trans(T_k));
        A2 = prod(T_k,p_k);
        A+= mPenaltyTangent*prod(A1, A2);
        //A+= mPenaltyNormal*prod(prod(trans(p_k),trans(N_k)),prod(N_k,p_k));
        //A+= mPenaltyTangent*prod(prod(prod(trans(p_k),trans(T_k)),T_k),p_k);

        b-= mPenaltyNormal*prod(trans(p_k),trans(N_k))*itPatchNode->GetValue(CONTACT_PRESSURE);
    }

    //neighboring nodes:
    
    for( auto it_neighbour_nodes = itPatchNode->GetValue(NEIGHBOUR_NODES).begin(); it_neighbour_nodes != itPatchNode->GetValue(NEIGHBOUR_NODES).end(); it_neighbour_nodes++) {
        if (it_neighbour_nodes->Has(CONTACT_PRESSURE)){
            p_k(0,1)= it_neighbour_nodes->X()-itPatchNode->X();
            p_k(0,2)= it_neighbour_nodes->Y()-itPatchNode->Y();
            p_k(1,4)= it_neighbour_nodes->X()-itPatchNode->X();;
            p_k(1,5)= it_neighbour_nodes->Y()-itPatchNode->Y();
            p_k(2,7)= it_neighbour_nodes->X()-itPatchNode->X();;
            p_k(2,8)= it_neighbour_nodes->Y()-itPatchNode->Y();
            N_k(0,0) = it_neighbour_nodes->GetValue(NORMAL)[0]*it_neighbour_nodes->GetValue(NORMAL)[0];
            N_k(0,1) = it_neighbour_nodes->GetValue(NORMAL)[1]*it_neighbour_nodes->GetValue(NORMAL)[1];
            N_k(0,2) = 2*it_neighbour_nodes->GetValue(NORMAL)[0]*it_neighbour_nodes->GetValue(NORMAL)[1];
            T_k(0,0) = it_neighbour_nodes->GetValue(NORMAL)[0]*it_neighbour_nodes->GetValue(NORMAL)[1];
            T_k(0,1) = -it_neighbour_nodes->GetValue(NORMAL)[0]*it_neighbour_nodes->GetValue(NORMAL)[1];
            T_k(0,2) = it_neighbour_nodes->GetValue(NORMAL)[1]*it_neighbour_nodes->GetValue(NORMAL)[1]-it_neighbour_nodes->GetValue(NORMAL)[0]*it_neighbour_nodes->GetValue(NORMAL)[0];

            A1 = prod(trans(p_k),trans(N_k));
            A2 = prod(N_k,p_k);
            A+= mPenaltyNormal*prod(A1, A2);

            A1 = prod(trans(p_k),trans(T_k));
            A2 = prod(T_k,p_k);
            A+= mPenaltyTangent*prod(A1, A2);

            b+= mPenaltyNormal*prod(trans(p_k),trans(N_k))*i_neighbour_node->GetValue(CONTACT_PRESSURE);
        }
    }*/

    // Computing coefficients a: A*a=b
    //UblasSpace<double,CompressedMatrix,Vector> U1 = UblasSpace<double, Matrix,Vector>();
    //UblasSpace<double, Matrix,Vector> U2 = UblasSpace<double, Matrix,Vector>();
    SkylineLUFactorizationSolver< UblasSpace<double,CompressedMatrix,Vector>, UblasSpace<double,Matrix,Vector>> solver = SkylineLUFactorizationSolver< UblasSpace<double,CompressedMatrix,Vector>, UblasSpace<double,Matrix,Vector>>();
    //std::cout<<A<<std::endl;
    /*CompressedMatrix compA(9,9,81);// = A.sparseView();
    for (unsigned i = 0; i < compA.size1 (); ++ i){
        for (unsigned j = 0; j < compA.size2 (); ++ j)
            compA (i, j) = 9 * i + j;
    }*/
    
    Vector coeff(mSigmaSize*(TDim+1));
    Vector b_vector = MatrixColumn(b,0);
    solver.Solve(A,coeff,b_vector);

    for (unsigned int j = 0; j < mSigmaSize;++j){    
        p_k(j,j*(TDim + 1) + 1)= itNode->X() - itPatchNode->X();
        p_k(j,j*(TDim + 1) + 2)= itNode->Y() - itPatchNode->Y();
        if (TDim == 3)
            p_k(j,j*(TDim + 1) + 3)= itNode->Z() - itPatchNode->Z();
    }
    
    Matrix coeff_matrix(mSigmaSize*(TDim + 1), 1);
    for (unsigned int i=0; i<mSigmaSize*(TDim + 1); ++i)
        coeff_matrix(i,0)=coeff(i);
    
    sigma = prod(p_k,coeff_matrix);

    rSigmaRecovered = MatrixColumn(sigma,0);
    
    if(mEchoLevel > 1)
        std::cout<<" Recovered pressure: "<< prod(N_k,sigma) <<", LM: "<<itNode->GetValue(CONTACT_PRESSURE)<<std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim> 
void SPRMetricProcess<TDim>::ComputeElementSize(ElementItType itElement){

    auto& this_geometry = itElement->GetGeometry(); 
    
    // Triangular elements
    if (this_geometry.size()==3){
        itElement->SetValue(ELEMENT_H, 2.0 * this_geometry.Circumradius());
    }

    // Tetrahedral elements
    if(this_geometry.size() == 4){
        itElement->SetValue(ELEMENT_H,std::pow(12.0 * GeometryUtils::CalculateVolume3D(this_geometry)/std::sqrt(2.0), 1.0/3.0));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class SPRMetricProcess<2>;
template class SPRMetricProcess<3>;

};// namespace Kratos.
