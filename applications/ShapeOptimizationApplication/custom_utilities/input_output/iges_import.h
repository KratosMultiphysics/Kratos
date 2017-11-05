//  KratosEmpireSurfaceMappingApplication
//
//  License:		 BSD License
//					 license: KratosEmpireSurfaceMappingApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#if !defined(KRATOS_IGES_EXTERNAL_MODEL_IMPORT_FUNCTION)
#define KRATOS_IGES_EXTERNAL_MODEL_IMPORT_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/external_model_import.h"

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <IGESControl_Reader.hxx>

namespace Kratos
{
///@addtogroup KratosEmpireSurfaceMappingApplication
///@{

///@name Kratos Classes
///@{

/// A response function for drag.
template <unsigned int TDim>
class IGESExternalModelImport : public ExternalModelImport
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(IGESExternalModelImport);

    typedef ExternalModelImport BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    IGESExternalModelImport(Parameters& rParameters)
      : ExternalModelImport(rParameters, "iges,igs")
    {
    }

    /// Destructor.
    ~IGESExternalModelImport() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void ReadModelPart() override
    {
        KRATOS_TRY;

        IGESControl_Reader reader();
        reader.ReadFile(mFilename);
        reader.TransferRoots();
        
        for (IndexType i=0; i<reader.NbShapes(); i++)
        {
            this->mCompoundsList.push(reader.Shape(i+1));
        }

        std::cout<<" Successfully read "<<reader.NbShapes()<<" shapes."<<std::endl;

        KRATOS_CATCH("");
    }

    void WriteModelPart( std::string FileName ) override
    {
        KRATOS_TRY;

        

        KRATOS_CATCH("");
    }
    ///@}

protected:
    ///@name Protected member Variables
    ///@{

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

#endif /* KRATOS_DRAG_RESPONSE_FUNCTION defined */
