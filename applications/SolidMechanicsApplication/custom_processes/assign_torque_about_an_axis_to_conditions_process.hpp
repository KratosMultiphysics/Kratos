//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:             January 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_ASSIGN_TORQUE_ABOUT_AN_AXIS_TO_CONDITIONS_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_TORQUE_ABOUT_AN_AXIS_TO_CONDITIONS_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/beam_math_utilities.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignTorqueAboutAnAxisToConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignTorqueAboutAnAxisToConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignTorqueAboutAnAxisToConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{


    AssignTorqueAboutAnAxisToConditionsProcess(ModelPart& model_part) : Process(Flags()) , mrModelPart(model_part) {}



    AssignTorqueAboutAnAxisToConditionsProcess(ModelPart& model_part,
					       Parameters rParameters
 				              ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "modulus" : 1.0,
                "direction" : [],
                "center" : []
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();


	if( KratosComponents< Variable<array_1d<double, 3> > >::Has( mvariable_name ) ) //case of array_1d variable
        {

	    mvalue = rParameters["modulus"].GetDouble();

	    for( unsigned int i=0; i<3; i++)
	      {
		mdirection[i] = rParameters["direction"][i].GetDouble();
		mcenter[i] = rParameters["center"][i].GetDouble();
	      }

	    double norm = norm_2(mdirection);
	    if(norm!=0)
		mdirection/=norm;

	}
	else //case of other variable type
        {
	  KRATOS_ERROR << "trying to set a variable that is not in the model_part - variable name is " << mvariable_name <<std::endl;
	}

        KRATOS_CATCH("");
    }


    /// Destructor.
    ~AssignTorqueAboutAnAxisToConditionsProcess() override {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the AssignTorqueAboutAnAxisToConditionsProcess algorithms.
    void Execute() override
    {

        KRATOS_TRY;

	this->AssignTorqueAboutAnAxis(KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name));

        KRATOS_CATCH("");

    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
    {
	array_1d<double,3> vector_value;
	vector_value.clear();
	InternalAssignValue(KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name), vector_value);
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
    std::string Info() const override
    {
        return "AssignTorqueAboutAnAxisToConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignTorqueAboutAnAxisToConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


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

    ModelPart& mrModelPart;
    std::string mvariable_name;
    double mvalue;
    array_1d<double,3> mdirection;
    array_1d<double,3> mcenter;

    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    AssignTorqueAboutAnAxisToConditionsProcess(AssignTorqueAboutAnAxisToConditionsProcess const& rOther);

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
    ///@}
    ///@name Private Operators
    ///@{


    void InternalAssignValue(const Variable<array_1d<double,3> >& rVariable,
			     const array_1d<double,3>& rvector_value)
    {
        const int nconditions = mrModelPart.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

            #pragma omp parallel for
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                it->SetValue(rVariable, rvector_value);
            }
        }
    }

    void AssignTorqueAboutAnAxis(const Variable<array_1d<double,3> >& rVariable)
    {
	KRATOS_TRY

	const int nconditions = mrModelPart.GetMesh().Conditions().size();

        if(nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = mrModelPart.GetMesh().ConditionsBegin();

	    std::vector<array_1d<double,3> > Couples(nconditions);
	    std::vector<array_1d<double,3> > Forces(nconditions);

	    Matrix rotation_matrix;
	    array_1d<double,3> radius;
	    array_1d<double,3> distance;
	    array_1d<double,3> force;

            #pragma omp parallel for private(rotation_matrix,radius,distance,force)
            for(int i = 0; i<nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

		Geometry< Node<3> >& rGeometry = it->GetGeometry();

		//Get geometry size
		unsigned int size  = rGeometry.size();
		array_1d<double,3> couple;
		couple.clear();
		array_1d<double,3> moment;
		moment.clear();
		for ( unsigned int j = 0; j < size; j++ )
		{

		    noalias(distance) = rGeometry[j].GetInitialPosition() - mcenter;

		    noalias(radius)  = distance-inner_prod(distance,mdirection) * mdirection,

		    //compute the skewsymmmetric tensor for the torque axis
		    BeamMathUtils<double>::VectorToSkewSymmetricTensor(mdirection, rotation_matrix);

		    double norm_radius = norm_2(radius);
		    if(norm_radius!=0)
			radius/=norm_radius;

		    noalias(force) = prod(rotation_matrix, radius);

		    noalias(couple) += force*norm_radius;

		    double norm = norm_2(force);
		    if(norm!=0)
			force/=norm;

		    //compute the skewsymmmetric tensor for the radius
		    BeamMathUtils<double>::VectorToSkewSymmetricTensor(radius, rotation_matrix);

		    noalias(moment) += norm_radius * prod(rotation_matrix, force);
		}

		const unsigned int dimension = rGeometry.WorkingSpaceDimension();
		double domain_size = 1.0;
		if(dimension==3)
		    domain_size = rGeometry.Area();
		if(dimension==2)
		    domain_size = rGeometry.Length();

		Couples[i] = moment * domain_size * (1.0/double(size));
		Forces[i]  = couple * (1.0/double(size));

            }

	    double total_size = 1.0;
	    array_1d<double,3> torque;
	    for(int i = 0; i<nconditions; i++)
	    {
		ModelPart::ConditionsContainerType::iterator it = it_begin + i;
		Geometry< Node<3> >& rGeometry = it->GetGeometry();
		const unsigned int dimension = rGeometry.WorkingSpaceDimension();
		double domain_size = 0.0;
		if(dimension==3)
		    domain_size = rGeometry.Area();
		if(dimension==2)
		    domain_size = rGeometry.Length();

		total_size += domain_size;
		torque += Couples[i];
	    }

	    torque /=total_size;

	    //solve distributed couple
	    double value = 0;
	    for(int i = 0; i<3; i++)
	    {
		if( torque[i] != 0 )
		    value = mdirection[i]*mvalue/torque[i];
		if( value != 0 )
		    break;
	    }

	    array_1d<double,3> load;
            #pragma omp parallel for private(torque)
	    for(int i = 0; i<nconditions; i++)
	    {
		ModelPart::ConditionsContainerType::iterator it = it_begin + i;
		load = value * Forces[i];
		it->SetValue(rVariable, load);
	    }

	}

      KRATOS_CATCH( "" )

    }


    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignTorqueAboutAnAxisToConditionsProcess& operator=(AssignTorqueAboutAnAxisToConditionsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignTorqueAboutAnAxisToConditionsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignTorqueAboutAnAxisToConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignTorqueAboutAnAxisToConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_TORQUE_ABOUT_AN_AXIS_TO_CONDITIONS_PROCESS_H_INCLUDED  defined
