//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED)
#define KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/linear_master_slave_constraint.h"
#include "utilities/binbased_fast_point_locator_conditions.h"


namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

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

/// Short class definition.
/** Detail class definition.
*/
class RVEPeriodicityUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RVEPeriodicityUtility
    KRATOS_CLASS_POINTER_DEFINITION(RVEPeriodicityUtility);

    typedef std::tuple< std::vector<unsigned int>, std::vector<double>, Vector > DataTuple;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RVEPeriodicityUtility(ModelPart& rDestinationModelPart) : mrModelPart(rDestinationModelPart){}

    /// Destructor.
    virtual ~RVEPeriodicityUtility(){}

    ///@}
    ///@name Operators
    ///@{
    void AssignPeriodicity(ModelPart& master,
                           ModelPart& slave,
                           const Matrix& strain,
                           const Vector& direction
                           )
    {
        if(master.Conditions().size() == 0)
            KRATOS_ERROR << "the master is expected to have conditions and it is empty" << std::endl;

        Vector translation = prod(strain,direction);

        BinBasedFastPointLocatorConditions<3> bin_based_point_locator(master);
        bin_based_point_locator.UpdateSearchDatabase();



        int max_search_results = 100;
        double search_tolerance = 1e-6;

        //construct auxiliary data structure to contain the master slave relation.
        //slave nodes must appear once, however a non-circular dependency is allowed between the masters
        for(unsigned int i=0; i<slave.Nodes().size(); ++i)
        {
            //search in which condition it falls
            auto i_node = slave.NodesBegin() + i;

            KRATOS_WATCH(i_node->Id());

            Condition::Pointer p_host_cond;
            Vector N;
            array_1d<double, 3 > transformed_slave_coordinates = i_node->Coordinates() - direction;


            // Finding the host element for this node
            const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(transformed_slave_coordinates, N, p_host_cond, max_search_results, search_tolerance);
            if(is_found)
            {
                const auto& geom = p_host_cond->GetGeometry();

                DataTuple aux_data;
                
                auto& T = std::get<2>(aux_data);
                T = translation;

                auto& master_ids = std::get<0>(aux_data);
                auto& weights = std::get<1>(aux_data);
                for(unsigned int j=0; j<geom.size(); ++j)
                {
                    master_ids.push_back(geom[j].Id());
                    weights.push_back(N[j]);
                }

                if(mAuxPairings.find(i_node->Id()) == mAuxPairings.end()) //this slave is not already present
                    mAuxPairings[i_node->Id()] = aux_data;
                else
                {
                    std::cout << "slave model part = " << slave << std::endl;
                    std::cout << "master model part = " << master << std::endl;
                    KRATOS_ERROR << "attempting to add twice the slave node with Id " << i_node->Id() << std::endl;
                }
            }
            else
            {
                KRATOS_ERROR << "counterpart not found for slave node " << i_node->Id() << std::endl;
            }
        }
    }


    void Finalize()
    {

        std::cout << "BEFORE COMPACTING" << std::endl;
        for(auto& data : mAuxPairings)
        {
            const unsigned int slave_id = data.first;
            auto& master_data = data.second;
            auto& master_ids = std::get<0>(master_data);
            auto& master_weights = std::get<1>(master_data);
            auto& T = std::get<2>(master_data);

            std::cout << slave_id << " - ";
            for(auto& master : master_ids)
                std::cout << master << " ";
            std::cout << " - ";
            for(auto& w : master_weights)
                std::cout << w << " ";  
            std::cout << " - " << T << std::endl; 
        }
       
        for(auto& data : mAuxPairings)
        {
            const unsigned int slave_id = data.first;
            auto& master_data = data.second;
            auto& master_ids = std::get<0>(master_data);
            auto& master_weights = std::get<1>(master_data);
            auto& T = std::get<2>(master_data);

            std::vector<unsigned int> final_master_ids;
            std::vector<double> final_master_weights;
            Vector final_T = T;

            for(unsigned int i = 0; i<master_ids.size(); ++i)
            {
                AppendIdsAndWeights( mAuxPairings, master_ids[i], master_weights[i],final_master_ids, final_master_weights, final_T);
            }

            //assign back the finalized pairings and weights to the data structure
            master_ids = final_master_ids;
            master_weights = final_master_weights;
            T = final_T;
        }

        //first assign master and slave all to false
        for(auto& node : mrModelPart.Nodes())
        {
            node.Set(SLAVE,false);
            node.Set(MASTER,false);
        }



        //compute the max id of the constraint
        int ConstraintId = 0;
        if(mrModelPart.MasterSlaveConstraints().size() != 0)
            ConstraintId = (mrModelPart.MasterSlaveConstraints().end()-1)->Id();
        ConstraintId+= 1;

        for(const auto& data : mAuxPairings)
        {
            const unsigned int slave_id = data.first;
            const auto& master_data = data.second;
            auto& master_ids = std::get<0>(master_data);
            auto& master_weights = std::get<1>(master_data);
            auto& T = std::get<2>(master_data);

            std::cout << slave_id << " - ";
            for(auto& master : master_ids)
                std::cout << master << " ";
            std::cout << " - ";
            for(auto& w : master_weights)
                std::cout << w << " ";  
            std::cout << " - " << T << std::endl;  

            //flag slave and master nodes
            mrModelPart.pGetNode(slave_id)->Set(SLAVE);
            for(auto& id : master_ids)
                mrModelPart.pGetNode(id)->Set(MASTER);
            
            //now construct the linear constraint
            LinearMasterSlaveConstraint::DofPointerVectorType slave_dofs, master_dofs;

            auto pslave_node = mrModelPart.pGetNode(slave_id);

            //define translation vector
            Vector xtranslation(1);
            Vector ytranslation(1);
            Vector ztranslation(1);
            xtranslation[0] = T[0];
            ytranslation[0] = T[1];
            ztranslation[0] = T[2];

            //define relation matrix (same for the different components)
            Matrix relation_matrix(1,master_weights.size());
            for(unsigned int i=0; i<relation_matrix.size2(); ++i)
                relation_matrix(0,i) = master_weights[i];

            //x-component
            slave_dofs.push_back(pslave_node->pGetDof(DISPLACEMENT_X));
            for(unsigned int i=0; i<master_ids.size(); ++i)
                master_dofs.push_back(mrModelPart.pGetNode(master_ids[i])->pGetDof(DISPLACEMENT_X));
            mrModelPart.AddMasterSlaveConstraint( 
                    Kratos::make_shared<LinearMasterSlaveConstraint>(ConstraintId,master_dofs,slave_dofs, relation_matrix, xtranslation) );
            ConstraintId+= 1;
            slave_dofs.clear();
            master_dofs.clear();

            //y-component
            slave_dofs.push_back(pslave_node->pGetDof(DISPLACEMENT_Y));
            for(unsigned int i=0; i<master_ids.size(); ++i)
                master_dofs.push_back(mrModelPart.pGetNode(master_ids[i])->pGetDof(DISPLACEMENT_Y));
            mrModelPart.AddMasterSlaveConstraint( 
                    Kratos::make_shared<LinearMasterSlaveConstraint>(ConstraintId,master_dofs,slave_dofs, relation_matrix, ytranslation) );
            ConstraintId+= 1;
            slave_dofs.clear();
            master_dofs.clear();

            //z-component
            slave_dofs.push_back(pslave_node->pGetDof(DISPLACEMENT_Z));
            for(unsigned int i=0; i<master_ids.size(); ++i)
                master_dofs.push_back(mrModelPart.pGetNode(master_ids[i])->pGetDof(DISPLACEMENT_Z));
            mrModelPart.AddMasterSlaveConstraint( 
                    Kratos::make_shared<LinearMasterSlaveConstraint>(ConstraintId,master_dofs,slave_dofs, relation_matrix, ztranslation) );
            ConstraintId+= 1;
            slave_dofs.clear();
            master_dofs.clear();
            
        }    
    }


    ///@}
    ///@name Operations
    ///@{


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
        std::stringstream buffer;
        buffer << "RVEPeriodicityUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "RVEPeriodicityUtility";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{
    void AppendIdsAndWeights(
                            std::map<unsigned int, DataTuple>& aux, 
                            const unsigned int master_id, 
                            const double& master_weight,
                            std::vector<unsigned int>& final_master_ids, 
                            std::vector<double>& final_master_weights, 
                            Vector& final_T)
    {
        if(std::abs(master_weight) > 1e-12  ) //discard nodes with negligible weight (note that weights sum to 1)
        {
            if(aux.find(master_id) == aux.end()) //master is NOT also a slave
            {
                final_master_ids.push_back(master_id);
                final_master_weights.push_back(master_weight);
            }
            else //master also happens to be a slave
            {
                const auto& other_data = aux[master_id];
                const auto& other_master_ids = std::get<0>(other_data);
                const auto& other_master_weights = std::get<1>(other_data);
                const auto& other_T = std::get<2>(other_data);
                for(unsigned int j = 0; j<other_master_ids.size(); ++j)
                {     
                    AppendIdsAndWeights(aux,other_master_ids[j], master_weight*other_master_weights[j], final_master_ids, final_master_weights, final_T );                   
                }

                final_T += master_weight*other_T;
            }
        }
    }


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
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    ModelPart& mrModelPart;

    std::map< unsigned int, DataTuple  > mAuxPairings;

    


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
    RVEPeriodicityUtility& operator=(RVEPeriodicityUtility const& rOther) = delete;

    /// Copy constructor.
    RVEPeriodicityUtility(RVEPeriodicityUtility const& rOther) = delete;

    ///@}

}; // Class RVEPeriodicityUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                RVEPeriodicityUtility& rThis){return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const RVEPeriodicityUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_RVE_PERIODICITY_UTILITY_H_INCLUDED  defined


