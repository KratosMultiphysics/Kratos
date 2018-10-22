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
//                   Pooyan Dadvand
//


#if !defined(KRATOS_MODEL_H_INCLUDED )
#define  KRATOS_MODEL_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@addtogroup KratosCore
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

/**
* @class Model
* @ingroup KratosCore
* @brief This class aims to manage different model parts across multi-physics simulations
* @details The class behaves as a manager of the different model parts. It uses unordered_maps of the variables and the model parts for that purpose
* @author Riccardo Rossi
*/
class KRATOS_API(KRATOS_CORE) Model
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the index type
    typedef ModelPart::IndexType IndexType;

    /// Pointer definition of Model
    KRATOS_CLASS_POINTER_DEFINITION(Model);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Model(){};

    /// Destructor.
    virtual ~Model()
    {
        mRootModelPartMap.clear();
        //mListOfVariablesLists.clear(); //this has to be done AFTER clearing the RootModelParts
    }

    Model & operator=(const Model&) = delete;
    Model(const Model&) = delete;


    ///@}
    ///@name Operators
    ///@{
    void Reset();

    ModelPart& CreateModelPart( const std::string ModelPartName, IndexType NewBufferSize=1 );

    void DeleteModelPart( const std::string ModelPartName );

    void RenameModelPart( const std::string OldName, const std::string NewName );

    ModelPart& GetModelPart(const std::string& rFullModelPartName);

    bool HasModelPart(const std::string& rFullModelPartName);

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
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


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
    std::map< std::string, std::unique_ptr<ModelPart> > mRootModelPartMap;

    std::set< std::unique_ptr<VariablesList> >& GetListOfVariableLists() const
    {
    static std::set< std::unique_ptr<VariablesList> > mListOfVariablesLists;
    return mListOfVariablesLists;
    }
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        //we construct auxiliary arrays to avoid having to serialize sets and maps of unique_ptrs
        std::vector<VariablesList* > aux_var_lists;
        std::vector<std::string> aux_names;
        std::vector<ModelPart* > aux_model_part_pointers;
        aux_var_lists.reserve(GetListOfVariableLists().size());
        aux_names.reserve(mRootModelPartMap.size());
        aux_model_part_pointers.reserve(mRootModelPartMap.size());

        for(auto it = mRootModelPartMap.begin(); it!=mRootModelPartMap.end(); ++it)
        {
            aux_names.push_back(it->first);
            aux_model_part_pointers.push_back((it->second).get());
        }

        for(auto it = GetListOfVariableLists().begin(); it!=GetListOfVariableLists().end(); ++it)
            aux_var_lists.push_back(it->get());

        rSerializer.save("ListOfVariablesLists", aux_var_lists);
        rSerializer.save("ModelPartNames", aux_names);
        rSerializer.save("ModelPartPointers", aux_model_part_pointers);
    }

    void load(Serializer& rSerializer)
    {
        //we construct auxiliary arrays to avoid having to serialize sets and maps of unique_ptrs
        std::vector<VariablesList* > aux_var_lists;
        std::vector<std::string> aux_names;
        std::vector<ModelPart* > aux_model_part_pointers;

        rSerializer.load("ListOfVariablesLists", aux_var_lists);
        rSerializer.load("ModelPartNames", aux_names);
        rSerializer.load("ModelPartPointers", aux_model_part_pointers);

        for(IndexType i=0; i<aux_var_lists.size(); ++i) {
            auto p_aux_list = std::unique_ptr<VariablesList>(aux_var_lists[i]);
            GetListOfVariableLists().insert(std::move(p_aux_list)); //NOTE: the ordering may be changed since the pointers are changed, however it should not matter
        }

        for(IndexType i=0; i<aux_names.size(); ++i) {
            mRootModelPartMap.insert(std::make_pair(aux_names[i],std::unique_ptr<ModelPart>(aux_model_part_pointers[i])));
        }


    }


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    ModelPart* RecursiveSearchByName(const std::string& ModelPartName, ModelPart* pModelPart);

    std::vector<std::string> GetSubPartsList(const std::string& rFullModelPartName);


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
//       Model& operator=(Model const& rOther);

    /// Copy constructor.
//       Model(Model const& rOther);


    ///@}

}; // Class Model

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                Model& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const Model& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODEL_H_INCLUDED  defined
