//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_FILENAME_H_INCLUDED )
#define  KRATOS_FILENAME_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/io.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// This class reads from STL file format and creates triangular elements in given model_part
/** The current version only reads triangles from the STL and not higher order polygons
 * The nodes corresponging to given vertices are not collapsed
 * A SubModelPart for each additional solid block will be created
 * For definition STL format please check https://en.wikipedia.org/wiki/STL_(file_format)
 * A sample file format with 3 triangles:
   solid 3 triangles
   facet normal  1.000000 0.000000 0.000000 
       outer loop 
          vertex 0.1 -2.56114e-08 0.1
          vertex 0.1 -0.499156 -0.0352136
          vertex 0.1 -0.473406 -0.0446259
       endloop 
   endfacet 
   facet normal  1.000000 -0.000000 0.000000 
       outer loop 
          vertex 0.1 -0.473406 -0.0446259
          vertex 0.1 -0.447464 -0.0534931
          vertex 0.1 -2.56114e-08 0.1
       endloop 
   endfacet 
   facet normal  1.000000 0.000000 0.000000 
       outer loop 
          vertex 0.1 -0.6 0.1
          vertex 0.1 -0.524702 -0.0252604
          vertex 0.1 -0.499156 -0.0352136
       endloop 
   endfacet 
   endsolid 3 triangles

*/
class KRATOS_API(KRATOS_CORE) StlIO : public IO
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StlIO
    KRATOS_CLASS_POINTER_DEFINITION(StlIO);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with filename.
    StlIO(std::string const& Filename);

    /// Constructor with stream.
    StlIO(std::iostream* pInputStream);

    /// Destructor.
    virtual ~StlIO(){
        delete mpInputStream;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ReadModelPart(ModelPart & rThisModelPart) override;

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
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;


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

    std::iostream* mpInputStream;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ReadSolid(ModelPart & rThisModelPart);
    
    void ReadFacet(ModelPart & rThisModelPart);

    void ReadLoop(ModelPart & rThisModelPart);

    Point ReadPoint();

    void ReadKeyword(std::string const& Keyword);

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
    StlIO& operator=(StlIO const& rOther);

    /// Copy constructor.
    StlIO(StlIO const& rOther);


    ///@}

}; // Class StlIO

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                StlIO& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const StlIO& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
