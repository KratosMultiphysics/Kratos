//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:        July 2016 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS )
#define  KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "dam_application_variables.h"

namespace Kratos
{

class DamUpliftConditionLoadProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(DamUpliftConditionLoadProcess);
    
    typedef Table<double,double> TableType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamUpliftConditionLoadProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY
			 
        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed"                                              : false,
                "Modify"                                                : true,
                "Gravity_Direction"                                     : "Y",
                "Reservoir_Bottom_Coordinate_in_Gravity_Direction"      : 0.0,
                "Upstream_Coordinate"                                   : [0.0,0.0,0.0],
                "Downstream_Coordinate"                                 : [0.0,0.0,0.0],
                "Upstream_Longitudinal_Coordinate"                      : [0.0,0.0,0.0],
                "Spe_weight"                                            : 0.0,
                "Water_level"                                           : 10.0,
                "Drains"                                                : false,
                "Height_drain"                                          : 0.0,
                "Distance"                                              : 0.0,
                "Effectiveness"                                         : 0.0,
                "table"                                                 : 0
            }  )" );
            
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"];
        rParameters["Upstream_Coordinate"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
        
        mmesh_id = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();
        mis_fixed = rParameters["is_fixed"].GetBool();
        mgravity_direction = rParameters["Gravity_Direction"].GetString();
        mreference_coordinate = rParameters["Reservoir_Bottom_Coordinate_in_Gravity_Direction"].GetDouble();
        mspecific = rParameters["Spe_weight"].GetDouble();
        
        // Getting the values of the coordinates (reference value)
        mx_0.resize(3,false);
        mx_0[0] = rParameters["Upstream_Coordinate"][0].GetDouble();
        mx_0[1] = rParameters["Upstream_Coordinate"][1].GetDouble();
        mx_0[2] = rParameters["Upstream_Coordinate"][2].GetDouble();
        
        mx_1.resize(3,false);
        mx_1[0] = rParameters["Downstream_Coordinate"][0].GetDouble();
        mx_1[1] = rParameters["Downstream_Coordinate"][1].GetDouble();
        mx_1[2] = rParameters["Downstream_Coordinate"][2].GetDouble();
        
        mx_2.resize(3,false);
        mx_2[0] = rParameters["Upstream_Longitudinal_Coordinate"][0].GetDouble();
        mx_2[1] = rParameters["Upstream_Longitudinal_Coordinate"][1].GetDouble();
        mx_2[2] = rParameters["Upstream_Longitudinal_Coordinate"][2].GetDouble();        
        
        // Drains
        mdrain = rParameters["Drains"].GetBool();
        mheight_drain = rParameters["Height_drain"].GetDouble();
        mdistance_drain = rParameters["Distance"].GetDouble();
        meffectiveness_drain = rParameters["Effectiveness"].GetDouble();
        mwater_level = rParameters["Water_level"].GetDouble();
        
        mtime_unit_converter = mr_model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        mTableId = rParameters["table"].GetInt();
        
        if(mTableId != 0)
            mpTable = model_part.pGetTable(mTableId);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    virtual ~DamUpliftConditionLoadProcess() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute()
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize()
{
        
        KRATOS_TRY;
        
        //Defining necessary variables
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        boost::numeric::ublas::bounded_matrix<double,3,3> RotationMatrix;
                
        // Computing the rotation matrix accoding with the introduced points by the user
        this->CalculateRotationMatrix(RotationMatrix);
        array_1d<double,3> newCoordinate;
        array_1d<double,3> auxiliar_vector;
        array_1d<double,3> reference_vector;

        // Gravity direction for computing the hydrostatic pressure
        int direction;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else if( mgravity_direction == "Z")
            direction = 3;
            
        // Computing the reference vector (coordinates)    
        reference_vector = prod(RotationMatrix,mx_0);
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
            
            double ref_coord = mreference_coordinate + mwater_level;
            
            if( mdrain == true)
            {
				double coefficient_effectiveness = 1.0 - meffectiveness_drain;
				double aux_drain = coefficient_effectiveness *(mwater_level - mheight_drain)* ((mbase_dam-mdistance_drain)/mbase_dam) + mheight_drain;

				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = it->Coordinate(1);
                    auxiliar_vector[1] = it->Coordinate(2);
                    auxiliar_vector[2] = it->Coordinate(3);
                    
                    // Computing the new coordinates                                        
                    newCoordinate = prod(RotationMatrix,auxiliar_vector);

					if(mis_fixed)
					{
						it->Fix(var);
					}
		
                    // We compute the first part of the uplift law 
                    muplift_pressure = (mspecific*((ref_coord-aux_drain)- (it->Coordinate(direction))))*(1.0-((1.0/(mdistance_drain))*(fabs( (newCoordinate(0)) - reference_vector(0))))) + mspecific*aux_drain;
                    
                    // If uplift pressure is greater than the limit we compute the second part and we update the value
                        if(muplift_pressure <= mspecific*aux_drain)
                        {
                            muplift_pressure = (mspecific*((mreference_coordinate+aux_drain)- (it->Coordinate(direction))))*(1.0-((1.0/(mbase_dam - mdistance_drain))*(fabs( (newCoordinate(0)) - (reference_vector(0)+mdistance_drain)))));
                        }
                        
					if(muplift_pressure<0.0)
					{
						it->FastGetSolutionStepValue(var)=0.0;
					}
					else
					{
						it->FastGetSolutionStepValue(var) = muplift_pressure;
					}
				}
					
			}
			else
			{
				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = it->Coordinate(1);
                    auxiliar_vector[1] = it->Coordinate(2);
                    auxiliar_vector[2] = it->Coordinate(3);
                    
                    newCoordinate = prod(RotationMatrix,auxiliar_vector);

					if(mis_fixed)
					{
						it->Fix(var);
					}
				
					muplift_pressure = (mspecific*(ref_coord- (it->Coordinate(direction))))*(1.0-((1.0/mbase_dam)*(fabs(newCoordinate(0)-reference_vector(0)))));
                    
					if(muplift_pressure<0.0)
					{
						it->FastGetSolutionStepValue(var)=0.0;
					}
					else
					{
						it->FastGetSolutionStepValue(var) = muplift_pressure;
					}
				}
			}            
        }
        
        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep()
    {
        
        KRATOS_TRY;
        
        //Defining necessary variables
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        boost::numeric::ublas::bounded_matrix<double,3,3> RotationMatrix;
        
        // Getting the values of table in case that it exist        
        if(mTableId != 0)
        { 
            double time = mr_model_part.GetProcessInfo()[TIME];
            time = time/mtime_unit_converter;
            mwater_level = mpTable->GetValue(time);
        }
        
        // && mr_model_part.GetProcessInfo()[TIME] >= 1.0
        
        // Computing the rotation matrix accoding with the introduced points by the user
        this->CalculateRotationMatrix(RotationMatrix);
        array_1d<double,3> newCoordinate;
        array_1d<double,3> auxiliar_vector;
        array_1d<double,3> reference_vector;

        // Gravity direction for computing the hydrostatic pressure
        int direction;
        
        if( mgravity_direction == "X")
            direction = 1;
        else if( mgravity_direction == "Y")
            direction = 2;
        else if( mgravity_direction == "Z")
            direction = 3;
            
        // Computing the reference vector (coordinates)    
        reference_vector = prod(RotationMatrix,mx_0);
                
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
            
            double ref_coord = mreference_coordinate + mwater_level;
            
            if( mdrain == true)
            {
				double coefficient_effectiveness = 1.0 - meffectiveness_drain;
				double aux_drain = coefficient_effectiveness *(mwater_level - mheight_drain)* ((mbase_dam-mdistance_drain)/mbase_dam) + mheight_drain;

				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = it->Coordinate(1);
                    auxiliar_vector[1] = it->Coordinate(2);
                    auxiliar_vector[2] = it->Coordinate(3);
                    
                    // Computing the new coordinates                                        
                    newCoordinate = prod(RotationMatrix,auxiliar_vector);

					if(mis_fixed)
					{
						it->Fix(var);
					}
		
                    // We compute the first part of the uplift law 
                    muplift_pressure = (mspecific*((ref_coord-aux_drain)- (it->Coordinate(direction))))*(1.0-((1.0/(mdistance_drain))*(fabs( (newCoordinate(0)) - reference_vector(0))))) + mspecific*aux_drain;
                    
                    // If uplift pressure is greater than the limit we compute the second part and we update the value
                        if(muplift_pressure <= mspecific*aux_drain)
                        {
                            muplift_pressure = (mspecific*((mreference_coordinate+aux_drain)- (it->Coordinate(direction))))*(1.0-((1.0/(mbase_dam - mdistance_drain))*(fabs( (newCoordinate(0)) - (reference_vector(0)+mdistance_drain)))));
                        }
                        
					if(muplift_pressure<0.0)
					{
						it->FastGetSolutionStepValue(var)=0.0;
					}
					else
					{
						it->FastGetSolutionStepValue(var) = muplift_pressure;
					}
				}
					
			}
			else
			{
				#pragma omp parallel for
				for(int i = 0; i<nnodes; i++)
				{
					ModelPart::NodesContainerType::iterator it = it_begin + i;
                    
                    auxiliar_vector.resize(3,false);                    
                    auxiliar_vector[0] = it->Coordinate(1);
                    auxiliar_vector[1] = it->Coordinate(2);
                    auxiliar_vector[2] = it->Coordinate(3);
                    
                    newCoordinate = prod(RotationMatrix,auxiliar_vector);

					if(mis_fixed)
					{
						it->Fix(var);
					}
				
					muplift_pressure = (mspecific*(ref_coord- (it->Coordinate(direction))))*(1.0-((1.0/mbase_dam)*(fabs(newCoordinate(0)-reference_vector(0)))));
                    
					if(muplift_pressure<0.0)
					{
						it->FastGetSolutionStepValue(var)=0.0;
					}
					else
					{
						it->FastGetSolutionStepValue(var) = muplift_pressure;
					}
				}
			}            
        }
        
        KRATOS_CATCH("");
    }
    
    void CalculateRotationMatrix(boost::numeric::ublas::bounded_matrix<double,3,3>& rRotationMatrix)
    {
        KRATOS_TRY;
        
        //Unitary vector in uplift direction
        array_1d<double,3> V_uplift;      
        V_uplift = (mx_1 - mx_0);
        mbase_dam= norm_2(V_uplift);
        double inv_norm_uplift = 1.0/norm_2(V_uplift);
        V_uplift[0] *= inv_norm_uplift;
        V_uplift[1] *= inv_norm_uplift;
        V_uplift[2] *= inv_norm_uplift;
        
        //Unitary vector in longitudinal direction
        array_1d<double,3> V_longitudinal;
        V_longitudinal = (mx_2 - mx_0);
        double inv_norm_longitudinal = 1.0/norm_2(V_longitudinal);
        V_longitudinal[0] *= inv_norm_longitudinal;
        V_longitudinal[1] *= inv_norm_longitudinal;
        V_longitudinal[2] *= inv_norm_longitudinal;
        
        //Unitary vector in local z direction
        array_1d<double,3> V_normal;
        MathUtils<double>::CrossProduct( V_normal, V_uplift, V_longitudinal);
                
        //Rotation Matrix
        rRotationMatrix(0,0) = V_uplift[0];
        rRotationMatrix(0,1) = V_uplift[1];
        rRotationMatrix(0,2) = V_uplift[2];
        
        rRotationMatrix(1,0) = V_longitudinal[0];
        rRotationMatrix(1,1) = V_longitudinal[1];
        rRotationMatrix(1,2) = V_longitudinal[2];
        
        rRotationMatrix(2,0) = V_normal[0];
        rRotationMatrix(2,1) = V_normal[1];
        rRotationMatrix(2,2) = V_normal[2];
        
        KRATOS_CATCH( "" )
        
    }

    /// Turn back information as a string.
    std::string Info() const
    {
        return "DamUpliftConditionLoadProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DamUpliftConditionLoadProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    std::size_t mmesh_id;
    std::string mvariable_name;
    std::string mgravity_direction;
    bool mis_fixed;
    double mreference_coordinate;
    double mspecific;
    double mbase_dam;
    double mwater_level;
    bool mdrain;
    double mheight_drain;
    double mdistance_drain;
    double meffectiveness_drain;
    double muplift_pressure;
    Vector mx_0;
    Vector mx_1;
    Vector mx_2;
    double mtime_unit_converter;
    TableType::Pointer mpTable;
    int mTableId;   
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamUpliftConditionLoadProcess& operator=(DamUpliftConditionLoadProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  DamUpliftConditionLoadProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamUpliftConditionLoadProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_UPLIFT_CONDITION_LOAD_PROCESS defined */
