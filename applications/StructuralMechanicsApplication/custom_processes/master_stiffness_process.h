// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Sebastian Ares de Parga Regalado
//                   
//

#if !defined(KRATOS_MASTER_STIFFNESS_PROCESS)
#define KRATOS_MASTER_STIFFNESS_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_filesystem.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

///@}
///@name  Functions
///@{

/**
 * @class MasterStiffnessProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes the stiffness matrix for list of master nodes related to a list of slave surfaces
 * @details It takes into account all elements in the ModelPart
 *
 * @author Sebastian Ares de Parga Regalado
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) MasterStiffnessProcess : public Process
{
public:
    
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;

    KRATOS_CLASS_POINTER_DEFINITION(MasterStiffnessProcess);

///-----------------------------------------------------------------------------------------------------
    /// Constructor.
    MasterStiffnessProcess(ModelPart& rThisModelPart, Parameters rParameters) : Process() ,
            mrThisModelPart(rThisModelPart)
    { 
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_list": ["Surface_name_1","Surface_name_2"],
                "list_of_master_coordinates": [[0.0,0.0,0.0],[1.0,0.0,0.0]],
                "eps_perturbation": 1e-5
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        //Assigning the process' variables
        this->parameters = rParameters;
        this->eps_perturbation = this->parameters["eps_perturbation"].GetDouble();
        this->Value();
        this->inf_rot=1e-5; //Infinitesimal rotation
        this->master_coor = this->parameters["list_of_master_coordinates"].GetMatrix();
        this->slave_surface_name = this->parameters["model_part_list"].GetStringArray();
        this->number_master_nodes = this->slave_surface_name.size();
        this->dim = rThisModelPart.GetProcessInfo()[DOMAIN_SIZE];
        KRATOS_CATCH("");
    }

    /// Destructor.
    ~MasterStiffnessProcess() override {}

///-----------------------------------------------------------------------------------------------------

    /// This method is executed in order to initialize the current step
    void ExecuteInitializeSolutionStep() override;

    /// This method is executed in order to finalize the current step
    void ExecuteFinalizeSolutionStep() override;


    /// This function gives an infinitesimal displacement for the 3D model part
    void Value(){
        //Still need to make cheaper and decide the final criterion.
        KRATOS_TRY;

        double minimum = 1e30;
        double maximum = -1e30;
        NodesArrayType &r_nodes   = mrThisModelPart.Nodes();

        for(NodeType& node : r_nodes){
            double mini = std::min({node.X(),node.Y(),node.Z()});
            double maxi = std::max({node.X(),node.Y(),node.Z()});
            if(mini < minimum){
                minimum = mini;
            }
            if(maxi > maximum){
                maximum = maxi;
            }
        this->value = this->eps_perturbation*std::abs(maximum - minimum);
        }

        KRATOS_CATCH("");
    }

    ///Sets the displacements/rotations of the slave surfaces for the specific DoF of the step. 
    void ChangeVectorValues(){

        KRATOS_TRY;

        int current_step = mrThisModelPart.GetProcessInfo()[STEP]-1; 
        
        this->n = current_step/(2*this->dim);// n stands for the current salve surface.

        this->DoF = current_step-(2*this->dim*this->n);//Local surface DoF

        if(this->DoF==0){
            //Creating Master- Slave constraint for each interface (displacements)
            int k=0;
            if(this->n==0){
                this->nodes_per_master_part = Vector(this->number_master_nodes+1,0);//Contains the sum of the master nodes so later we can access when recreating the master-slave constraints
                this->Id_master = Vector(this->number_master_nodes);//Contains the the Id of the master nodes
            }
            int counter=0;
            for(int j=0;j<this->number_master_nodes;j++){//Creating the master nodes as fictitious nodes to create the master-slave constraints with the slaves surfaces.
                ModelPart &current_master_slave_part=mrThisModelPart.GetSubModelPart(this->slave_surface_name[j]);
                if(this->n==0){
                    current_master_slave_part.CreateNewNode(mrThisModelPart.NumberOfNodes()+1,this->master_coor(j,0), this->master_coor(j,1), this->master_coor(j,2));
                    VariableUtils().AddDofWithReaction(DISPLACEMENT_X,REACTION_X, current_master_slave_part);
                    VariableUtils().AddDofWithReaction(DISPLACEMENT_Y,REACTION_Y, current_master_slave_part);
                    VariableUtils().AddDofWithReaction(DISPLACEMENT_Z,REACTION_Z, current_master_slave_part);
                    int current_number_of_nodes = current_master_slave_part.NumberOfNodes()-1;
                    this->nodes_per_master_part(j+1)  =  this->dim*(counter+current_number_of_nodes);
                    this->Id_master(j) = mrThisModelPart.NumberOfNodes();
                    counter += current_number_of_nodes;
                }
                NodeType& master_node = current_master_slave_part.GetNode(int(this->Id_master(j)));
                NodesArrayType &slave_nodes = current_master_slave_part.Nodes();
                for(NodeType& slave_node : slave_nodes){
                    if(slave_node.Id()==master_node.Id()){
                        continue;
                    }
                    current_master_slave_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", k, master_node, DISPLACEMENT_X, slave_node, DISPLACEMENT_X, 1.0, 0);
                    k+=1;
                    current_master_slave_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", k, master_node, DISPLACEMENT_Y, slave_node, DISPLACEMENT_Y, 1.0, 0);
                    k+=1;
                    current_master_slave_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", k, master_node, DISPLACEMENT_Z, slave_node, DISPLACEMENT_Z, 1.0, 0);
                    k+=1;
                }
            }
            VariableUtils().SetFlag(TO_ERASE,false, mrThisModelPart.MasterSlaveConstraints());
        }
        if(this->DoF>2){
            ModelPart &current_master_slave_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[this->n]);
            VariableUtils().SetFlag(TO_ERASE, true, current_master_slave_part.MasterSlaveConstraints());
            mrThisModelPart.RemoveMasterSlaveConstraints(TO_ERASE);
            int k=this->nodes_per_master_part(this->n);
            Kratos::Vector v = Vector(3);
            Kratos::Matrix R = this->QuaternionRotationMatrix();
            NodeType &master_node = current_master_slave_part.GetNode(this->Id_master(this->n));
            NodesArrayType &slave_nodes = current_master_slave_part.Nodes();
            for(NodeType& slave_node : slave_nodes){//Loop on interface nodes to each constraint
                Kratos::Vector disp = Vector(this->dim,0);
                if(slave_node.Id()==master_node.Id()){ 
                    continue;   
                }
                v(0) = slave_node.X0() - this->master_coor(this->n,0);
                v(1) = slave_node.Y0() - this->master_coor(this->n,1);
                v(2) = slave_node.Z0() - this->master_coor(this->n,2);
                disp = prod(R,v)-v;
                current_master_slave_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", k, master_node, DISPLACEMENT_X, slave_node, DISPLACEMENT_X, disp(0), 0);
                k+=1;
                current_master_slave_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", k, master_node, DISPLACEMENT_Y, slave_node, DISPLACEMENT_Y, disp(1), 0);
                k+=1;
                current_master_slave_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", k, master_node, DISPLACEMENT_Z, slave_node, DISPLACEMENT_Z, disp(2), 0);
                k+=1;
            }
        }
        //Assign the current submodel slave part.
        ModelPart &slave_submodel_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[this->n]);
        NodeType& master_node = slave_submodel_part.GetNode(this->Id_master(this->n));
        
        //Setting the infinitesimal displacement/rotation boundary condition for current slave surface.
        if(this->DoF<3){//Displacement DoF(0-2)
            Kratos::Vector disp = Vector(3,0);
            disp(this->DoF) = this->value;
            master_node.FastGetSolutionStepValue(DISPLACEMENT) = disp;
            master_node.Fix(DISPLACEMENT_X);
            master_node.Fix(DISPLACEMENT_Y);
            master_node.Fix(DISPLACEMENT_Z);
        } else{
            Kratos::Vector disp = Vector(3,1);
            master_node.FastGetSolutionStepValue(DISPLACEMENT) = disp;
            master_node.Fix(DISPLACEMENT_X);
            master_node.Fix(DISPLACEMENT_Y);
            master_node.Fix(DISPLACEMENT_Z);
        }

        //All other surfaces are fixed.
        for(int i=0;i<this->number_master_nodes;i++){
            if (i!=this->n){
                ModelPart &fixed_submodel_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[i]);
                NodeType& master_node = fixed_submodel_part.GetNode(this->Id_master(i));
                Kratos::Vector disp = Vector(this->dim,0);
                master_node.FastGetSolutionStepValue(DISPLACEMENT) = disp;
                master_node.Fix(DISPLACEMENT_X);
                master_node.Fix(DISPLACEMENT_Y);
                master_node.Fix(DISPLACEMENT_Z);
            }
        }
        
        KRATOS_CATCH("");
    }

    //This function computes the rotation matrix for the current DoF by using quaternions.
    Kratos::Matrix QuaternionRotationMatrix(){
        Kratos::Vector u = Vector(3,0);
        u(this->DoF-3) = 1;
        Quaternion<double> q = Quaternion<double>::FromAxisAngle(u(0),u(1),u(2), this->inf_rot);
        q.normalize(); //Checks if the normalization of a Quaternion (make it a unit quaternion) is being calculated correctly.
        Kratos::Matrix rotation_matrix = Matrix(3,3);
        q.ToRotationMatrix(rotation_matrix);
        return rotation_matrix;
    }

    //This function computes the rotation matrix for the current DoF.
    Kratos::Matrix ComputeRotationMatrix(){
        Kratos::Matrix R = Matrix(this->dim,this->dim);
        Kratos::Vector u = Vector(this->dim,0);
        u(this->DoF-this->dim) = 1;
        R(0,0) = cos(this->inf_rot)+pow(u(0),2)*(1-cos(this->inf_rot));
        R(0,1) = u(0)*u(1)*(1-cos(this->inf_rot))-u(2)*sin(this->inf_rot);
        R(0,2) = u(0)*u(2)*(1-cos(this->inf_rot))+u(1)*sin(this->inf_rot);
        R(1,0) = u(1)*u(0)*(1-cos(this->inf_rot))+u(2)*sin(this->inf_rot);
        R(1,1) = cos(this->inf_rot)+pow(u(1),2)*(1-cos(this->inf_rot));
        R(1,2) = u(1)*u(2)*(1-cos(this->inf_rot))-u(0)*sin(this->inf_rot);
        R(2,0) = u(2)*u(0)*(1-cos(this->inf_rot))-u(1)*sin(this->inf_rot);
        R(2,1) = u(2)*u(1)*(1-cos(this->inf_rot))+u(0)*sin(this->inf_rot);
        R(2,2) = cos(this->inf_rot)+pow(u(2),2)*(1-cos(this->inf_rot));
        return R;
    }

    //Obtains the Stiffness column Vector of the current DoF of the n surface.
    void MasterStiffnessVector(){
        KRATOS_TRY;

        int current_step = mrThisModelPart.GetProcessInfo()[STEP]-1; 

        if(this->DoF==0 && this->n==0){//Only initialize once the master stiffness matrix 
            int matrix_size = this->number_master_nodes*2*this->dim;
            this->master_stiffness = Matrix(matrix_size,matrix_size,0);//Initialize Master Stiffness Matrix
        }

        Kratos::Matrix resultant = Matrix(this->number_master_nodes,2*this->dim,0);
        for(int k=0;k<this->number_master_nodes;k++){
            ModelPart &surface_model_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[k]);
            //Initializaing results Matrices an Vectors
            int surface_num_nodes = surface_model_part.NumberOfNodes();
            Kratos::Matrix surface_undeformed_coordinates = Matrix(surface_num_nodes,this->dim);
            Kratos::Matrix surface_reaction = Matrix(surface_num_nodes,this->dim);
            int counter = 0;

            NodesArrayType &surface_nodes   = surface_model_part.Nodes();

            for(NodeType& node : surface_nodes){
                if(counter==surface_num_nodes){
                    continue;
                }
                Kratos::Vector undeformed_coordinates = Vector(this->dim);//-----May be declared outside the loop
                undeformed_coordinates(0) = node.X0();
                undeformed_coordinates(1) = node.Y0();
                undeformed_coordinates(2) = node.Z0();
                for(int i=0;i<this->dim;i++){
                    surface_undeformed_coordinates(counter,i) = undeformed_coordinates(i); //Matrix of undeformed coordinates
                    surface_reaction(counter,i) = node.FastGetSolutionStepValue(REACTION)[i]; //Matrix of reactions

                }
                counter += 1;
            }

            Kratos::Matrix r = Matrix(surface_num_nodes,this->dim); // Initialize Position vector r

            for(int i=0; i<surface_num_nodes; i++){
                for(int j=0; j<this->dim; j++){
                    r(i,j) = surface_undeformed_coordinates(i,j)-this->master_coor(k,j);// Assign position vector r
                }
            }

            //Calculating the moments about the master node
            Kratos::Matrix moments = Cross_product(r,surface_reaction,surface_num_nodes);
            

            // Resultant assembling
            for(int i=0; i<surface_num_nodes; i++){
                for(int j=0; j<this->dim; j++){
                    resultant(k,j) += surface_reaction(i,j); // Resultant forces assembling
                    resultant(k,j+this->dim) += moments(i,j); // Resultant moments assembling
                }
            }

        }
        

        float constraint;
        if(this->DoF<3){
            constraint = this->value;
        }else{
            constraint = this->inf_rot;
        }

        //Stiffness Calculation
        for(int k=0; k<this->number_master_nodes;k++){
            for(int i=0; i<2*this->dim;i++){
                this->master_stiffness((2*k*this->dim)+i,current_step) = resultant(k,i)/constraint;
            }
        }
        
        
        if(this->DoF==5){
            VariableUtils().SetFlag(TO_ERASE, true, mrThisModelPart.MasterSlaveConstraints());
            mrThisModelPart.RemoveMasterSlaveConstraints(TO_ERASE);
            if(this->n==this->number_master_nodes-1){
                std::cout<<"Master Stiffness: "<<std::endl;
                std::cout<<this->master_stiffness<<std::endl;
                this->json_parameters.AddEmptyValue("Stiffness Matrix");
                this->json_parameters["Stiffness Matrix"].SetMatrix(this->master_stiffness);
                this->json_parameters.AddEmptyValue("Displacement");
                this->json_parameters["Displacement"].SetDouble(this->value);
                this->json_parameters.AddEmptyValue("Rotation");
                this->json_parameters["Rotation"].SetDouble(this->inf_rot);
                this->json_parameters.AddEmptyValue("Interface Frame Origins");
                this->json_parameters["Interface Frame Origins"].SetMatrix(this->master_coor);
                this->CreateJSONfile(); 
            }
        }

        KRATOS_CATCH("");
    }

    //This function creates a .json file containing the master stiffness matrices of the slave surfaces.
    void CreateJSONfile(){
        KRATOS_TRY;

        const std::string &r_json_text = this->json_parameters.PrettyPrintJsonString();
        std::filebuf buffer;
        buffer.open(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "Stiffness_Matrix.json"}), std::ios::out);
        std::ostream os(&buffer);
        os << r_json_text;
        buffer.close();
        
        KRATOS_CATCH("");
    }

    // Obtains the dot product (used to obtain the moments).
    Kratos::Matrix Cross_product(Kratos::Matrix a, Kratos::Matrix b, int n){
        Kratos::Matrix c = Matrix(n,this->dim);
        for(int i=0;i<n;i++){
            c(i,0) = a(i,1)*b(i,2)-a(i,2)*b(i,1); //Moments x
            c(i,1) = a(i,2)*b(i,0)-a(i,0)*b(i,2); //Moments y
            c(i,2) = a(i,0)*b(i,1)-a(i,1)*b(i,0); //Moments z
        }
        return c;
    }

    
    /// Turn back information as a string.
    std::string Info() const override
    {
        return "MasterStiffnessProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MasterStiffnessProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///-----------------------------------------------------------------------------------------------------

protected:
    Parameters parameters, json_parameters;
    float eps_perturbation,value,inf_rot;
    int number_master_nodes,dim,n,DoF;
    std::vector<std::string> slave_surface_name;
    Kratos::Matrix master_coor,master_stiffness;
    Kratos::Vector nodes_per_master_part,Id_master;
    
///-----------------------------------------------------------------------------------------------------

private:

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrThisModelPart; // The main model part

    /// Assignment operator.
    MasterStiffnessProcess& operator=(MasterStiffnessProcess const& rOther) = delete;

    /// Copy constructor.
    MasterStiffnessProcess(MasterStiffnessProcess const& rOther) = delete;


    ///@}

}; // Class MasterStiffnessProcess

///-----------------------------------------------------------------------------------------------------

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MasterStiffnessProcess& rThis);

/// output stream function
 inline std::ostream& operator << (std::ostream& rOStream,
                                   const MasterStiffnessProcess& rThis)
 {
     rThis.PrintInfo(rOStream);
     rOStream << std::endl;
     rThis.PrintData(rOStream);

     return rOStream;
 }

} //namesapce Kratos.
#endif /* MASTER_STIFFNESS_PROCESS defined */
