//  KratosEmpireSurfaceMappingApplication
//
//  License:		 BSD License
//					 license: KratosEmpireSurfaceMappingApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_EXTERNAL_MODEL_IMPORT_FUNCTION)
#define KRATOS_EXTERNAL_MODEL_IMPORT_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"

// OpenCASCADE includes
#include <gp_Vec.hxx>
#include <gp_Trsf.hxx>
#include <gp_Pnt.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopExp_Explorer.hxx>

#include <TopOpeBRepBuild_HBuilder.hxx>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>

#include <BRepAlgo_Cut.hxx>
#include <BRepAlgo.hxx>

#include <TDF_Data.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelMap.hxx>
#include <TDF_ChildIterator.hxx>
#include <TDF_MapIteratorOfLabelMap.hxx>

#include <TNaming_NamedShape.hxx>
#include <TNaming_Selector.hxx>
#include <TNaming_Tool.hxx>
#include <TNaming_Builder.hxx>
#include <TNaming.hxx>


#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <boost/algorithm/string.hpp>

// Application includes

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A base class for response functions.
class ExternalModelImport
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ExternalModelImport);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    ExternalModelImport( Parameters& rParameters,  std::string FileTypes)
    {
        KRATOS_TRY();

        Parameters default_params(R"(
            {
                "input_type": "PLEASE_SPECIFY_INPUT_MODEL_FILE_TYPE",
                "input_filename": "PLEASE_SPECIFY_INPUT_FILENAME",
                "echo_level": 0
            })");
        
        rParameters.ValidateAndAssignDefaults(default_params);

        std::vector<std::string> file_types;
        boost::split(file_types, FileTypes, ",");

        bool valid_file_type = false;
        for (std::string file_type : file_types)
            if (file_type.compare(rParameters["input_type"].GetString()) == 0)
            {
                valid_file_type = true;
                break;
            }
        
        if (!valid_file_type)
        {
            std::stringstream msg;                
            msg << "this module only supports "<<rFileTypes<<", provided input_type = "<<rParameters["input_type"].GetString()<<std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," invalid input_type: ",msg.str());
        }

        mFilename = rParameters["input_filename"].GetString();
        mEchoLevel = rParameters["echo_level"].GetInteger();
        mFileType = rParameters["input_type"].GetString();

        KRATOS_CATCH("");

    }

    /// Destructor.
    virtual ~ExternalModelImport()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    void Execute()
    {
        KRATOS_TRY;

        this->ReadModelPart();
        this->ConvertModelPartToNURBS();
        this->SewModelPart();

        KRATOS_CATCH("");
    }

    virtual void ReadModelPart() 
    {
        KRATOS_TRY;

        KRATOS_THROW_ERROR(std::invalid_argument, " not implemented: ", "ReadModelPart()");

        KRATOS_CATCH("")
    }

    virtual void WriteModelPart( std::string FileName ) 
    {
        KRATOS_TRY;

        KRATOS_THROW_ERROR(std::invalid_argument, " not implemented: ", "WriteModelPart()");

        KRATOS_CATCH("")
    }

    void ConvertModelPartToNURBS()
    {
        KRATOS_TRY;

        BRepBuilderAPI_NurbsConvert shape_converter();

        for (IndexType i=0; i< this->mCompoundsList.size(); i++)
        {
            if (echo_level > 0)
                std::cout<<"Converting shape: "<<mCompoundsList[i]<<" to NURBS"<<std::endl;
            
            shape_converter.Perform(mCompoundsList[i]);
            mCompoundNurbsList.push(shape_converter.Shape());
        }

        std::cout<<"Converted "<<mCompoundNurbsList.size()<<" shapes to NURBS."<<std::endl;

        KRATOS_CATCH("");
    }

    void SewModelPart()
    {
        KRATOS_TRY;
        
        for (IndexType i = 0; i<mCompoundNurbsList.size(); i++)
        {
            mSewedCompound.Add(mCompoundNurbsList[i]);
        }

        if (echo_level > 0)
        {
            std::cout<<"> --- Shape Summary Before Sewing ---"<<std::endl;
            mSewedCompound.Dump();
        }

        mSewedCompound.Perform();

        if (echo_level > 0)
        {
            std::cout<<"> --- Shape Summary After Sewing ---"<<std::endl;
            mSewedCompound.Dump();
        }

        std::cout<<"Sewed "<<mCompoundNurbsList.size()<<" shapes."<<std::endl;
        
        KRATOS_CATCH("");        
    }

    const std::vector<TopoDS_Shape>& GetCompoundsList() const
    {
      return mCompoundsList;
    }

    const std::vector<TopoDS_Shape>& GetCompoundNURBSList() const
    {
      return mCompoundNurbsList;
    }

    const BRepBuilderAPI_Sewing& GetSewedShape() const
    {
      return mSewedCompound;
    }
    ///@}

protected:
    ///@name Protected member Variables
    ///@{
    int mEchoLevel;
    std::string mFilename;
    std::string mFileType;
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}

private:
    ///@name Member Variables
    ///@{
    std::vector<TopoDS_Shape> mCompoundsList;
    std::vector<TopoDS_Shape> mCompoundNurbsList;
    BRepBuilderAPI_Sewing mSewedCompound();
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_RESPONSE_FUNCTION defined */
