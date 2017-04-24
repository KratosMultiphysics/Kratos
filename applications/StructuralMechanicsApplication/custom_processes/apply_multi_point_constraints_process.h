//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala 
//
//


#ifndef APPLY_MULTI_POINT_CONSTRAINTS_PROCESS_H
#define APPLY_MULTI_POINT_CONSTRAINTS_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{

class ApplyMultipointConstraintsProcess : public Process
{
public:

    /// Pointer definition of MoveRotorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyMultipointConstraintsProcess);

    typedef MpcData::Pointer MpcDataPointerType;
    typedef Dof<double>* DofPointerType;
    typedef Dof<double> DofType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef ProcessInfo      ProcessInfoType;
    typedef ProcessInfo::Pointer      ProcessInfoPointerType;

    /// Constructor.
    ApplyMultipointConstraintsProcess(  ModelPart& model_part,
                                        Parameters rParameters
                                        ) : Process(Flags()) , mr_model_part(model_part), mpcDataMap()
    {

         Parameters default_parameters( R"(
            {
                "master_model_part_name":"default_master",
                "slave_model_part_name":"default_slave",
                "constraint_sets":["default"],
                "interpolation_type":"nearest_element"            
            }  )" );

        rParameters.ValidateAndAssignDefaults(default_parameters);
        //mrMpcData(model_part.GetValue(KratosComponents< Variable<MpcData> >::Get( "MPC_DATA" ))
        mpcDataMap["default"] = MpcDataPointerType( new MpcData() );
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        info->SetValue(MPC_POINTER, mpcDataMap["default"]);
    }

    ApplyMultipointConstraintsProcess(  ModelPart& model_part
                                        ) : Process(Flags()) , mr_model_part(model_part), mpcDataMap()
    {

        //mrMpcData(model_part.GetValue(KratosComponents< Variable<MpcData> >::Get( "MPC_DATA" ))
        mpcDataMap["default"] = MpcDataPointerType( new MpcData() );
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        info->SetValue(MPC_POINTER, mpcDataMap["default"]);

    }    

    void AddMasterSlaveRelation(Node<3> &MasterNode, VariableComponentType& MasterVariable, Node<3> &SlaveNode, VariableComponentType& SlaveVariable, double weight, int PartitionId=0)
    {
        SlaveNode.Set(SLAVE);        
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
    	DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, PartitionId);
    }

    void AddMasterSlaveRelationWithDofs(DofType slaveDOF, DofType masterDOF, double masterWeight, int PartitionId=0 )
    {
        ProcessInfoType info = mr_model_part.GetProcessInfo();
        MpcDataPointerType pMpc = info[MPC_POINTER];
        pMpc->AddConstraint(slaveDOF, masterDOF,  masterWeight, PartitionId);
    }


    /// Destructor.
    virtual ~ApplyMultipointConstraintsProcess(){
        /*for(auto mpcDataMapElem : mpcDataMap){
            delete mpcDataMapElem.second;
        }*/
    }


    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;



        KRATOS_CATCH("");
    }


    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ApplyMultipointConstraintsProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ApplyMultipointConstraintsProcess";}

    /// Print object's data.
    void PrintData() {
        std::cout<<"Number of slave nodes :: "<<std::endl;
        mpcDataMap["default"]->GetInfo();
    }



protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    ModelPart&                                 mr_model_part;
    MpcDataMapType                             mpcDataMap;

private:

    /// Assignment operator.
    ApplyMultipointConstraintsProcess& operator=(ApplyMultipointConstraintsProcess const& rOther){return *this;}

}; // Class MoveRotorProcess

};  // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
