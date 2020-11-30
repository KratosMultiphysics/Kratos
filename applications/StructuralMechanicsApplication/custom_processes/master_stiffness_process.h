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
            double mini=std::min({node.X(),node.Y(),node.Z()});
            double maxi=std::max({node.X(),node.Y(),node.Z()});
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

        //Order the vector of strings to make the current slave surface be the first of the list.
        if(this->DoF==0){
            std::reverse(std::begin(this->slave_surface_name),std::begin(this->slave_surface_name)+this->n+1);
        }

        //Assign the current submodel slave part.
        ModelPart &slave_submodel_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[0]);
        NodesArrayType &slave_nodes   = slave_submodel_part.Nodes();
        //Setting the infinitesimal displacement/rotation boundary condition for current slave surface.
        if(this->DoF<3){//Displacement DoF(0-2)
            for(NodeType& node : slave_nodes){//Slave surface
                Kratos::Vector disp = Vector(this->dim,0);
                disp(this->DoF) = this->value;
                node.FastGetSolutionStepValue(DISPLACEMENT) = disp;
                node.Fix(DISPLACEMENT_X);
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
            }
        } else {//Rotation DoF(3-5)
            Kratos::Matrix R = this->ComputeRotationMatrix();
            Kratos::Vector v = Vector(this->dim);
            for(NodeType& node : slave_nodes){//Slave surface
                Kratos::Vector disp = Vector(3,0);
                v(0) = node.X0() - this->master_coor(this->n,0);
                v(1) = node.Y0() - this->master_coor(this->n,1);
                v(2) = node.Z0() - this->master_coor(this->n,2);
                disp = prod(R,v)-v; 
                node.FastGetSolutionStepValue(DISPLACEMENT) = disp;
                node.Fix(DISPLACEMENT_X);
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
            }
        }

        //All other surfaces are fixed.
        for(int i=1;i<this->number_master_nodes;i++){
            ModelPart &fixed_submodel_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[i]);
            NodesArrayType &fixed_nodes   = fixed_submodel_part.Nodes();
            //Setting the infinitesimal fixed boundary condition.
            for(NodeType& node : fixed_nodes){//Fixed surface
                Kratos::Vector disp = Vector(this->dim,0);
                node.FastGetSolutionStepValue(DISPLACEMENT) = disp;
                node.Fix(DISPLACEMENT_X);
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
            }
        }
        
        KRATOS_CATCH("");
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

        ModelPart &slave_model_part = mrThisModelPart.GetSubModelPart(this->slave_surface_name[0]);

        if(this->DoF==0){//Only initialize once the master stiffness matrix for each slave surface
            this->master_stiffness = Matrix(2*this->dim,2*this->dim,0);//Initialize Master Stiffness Matrix
        }

        //Initializaing results Matrices an Vectors
        int slave_num_nodes = slave_model_part.NumberOfNodes();
        Kratos::Matrix slave_undeformed_coordinates = Matrix(slave_num_nodes,this->dim);
        Kratos::Matrix slave_reaction = Matrix(slave_num_nodes,this->dim);
        int counter = 0;

        NodesArrayType &slave_nodes   = slave_model_part.Nodes();

        for(NodeType& node : slave_nodes){
            Kratos::Vector undeformed_coordinates = Vector(this->dim);
            undeformed_coordinates(0) = node.X0();
            undeformed_coordinates(1) = node.Y0();
            undeformed_coordinates(2) = node.Z0();
            for(int i=0;i<this->dim;i++){
                slave_undeformed_coordinates(counter,i) = undeformed_coordinates(i); //Matrix of undeformed coordinates
                slave_reaction(counter,i) = node.FastGetSolutionStepValue(REACTION)[i]; //Matirx of reactions

            }
            counter += 1;
        }

        Kratos::Matrix r = Matrix(slave_num_nodes,this->dim); // Initialize Position vector r

        for(int i=0; i<slave_num_nodes; i++){
            for(int j=0; j<this->dim; j++){
                r(i,j) = slave_undeformed_coordinates(i,j)-this->master_coor(this->n,j);// Assign position vector r
            }
        }

        //Calculating the moments about the master node
        Kratos::Matrix moments = Cross_product(r,slave_reaction,slave_num_nodes);
        
        Kratos::Vector resultant = Vector(2*this->dim,0);

        // Resultant assembling
        for(int i=0; i<slave_num_nodes; i++){
            for(int j=0; j<this->dim; j++){
                resultant(j) += slave_reaction(i,j); // Resultant forces assembling
                resultant(j+this->dim) += moments(i,j); // Resultant moments assembling
            }
        }

        float constraint;
        if(this->DoF<3){
            constraint = this->value;
        }else{
            constraint = this->inf_rot;
        }

        //Stiffness Calculation
        for(int i=0; i<this->dim;i++){
            master_stiffness(i,this->DoF) = resultant(i)/constraint;
            master_stiffness(i+this->dim,this->DoF) = resultant(i+this->dim)/constraint;
        }

        if(this->DoF==5){
            std::cout<<"Master Stiffness: "<<this->slave_surface_name[0]<<std::endl;
            std::cout<<this->master_stiffness<<std::endl;
            this->json_parameters.AddEmptyValue(this->slave_surface_name[0]);//Creates a value on the json_parameters object for the current slave's stiffnes matrix
            this->json_parameters[this->slave_surface_name[0]].SetMatrix(this->master_stiffness);//Sets the current slave's stiffness matrix
            //Reorders the list to its original form.
            std::reverse(std::begin(this->slave_surface_name),std::begin(this->slave_surface_name)+this->n+1);
            if(this->n==this->number_master_nodes-1){
                std::cout<<"-------------------------------------------------------------------------------------------------"<<std::endl;
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
        buffer.open(FilesystemExtensions::JoinPaths({FilesystemExtensions::CurrentWorkingDirectory(), "Master_Stiffness_Matrices.json"}), std::ios::out);
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
