//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_GID_EIGEN_IO_H_INCLUDED )
#define  KRATOS_GID_EIGEN_IO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/gid_io.h"


namespace Kratos
{
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

  /// GidIO extended for writting Eigenvalue Results

  /** The main functionality of this class is to write the custom format
  * where the result label is customizable, i.e. it can contain e.g.
  * the Eigenvalue or -frequency of the solution
  */
class GidEigenIO : public GidIO<>
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GidEigenIO
    KRATOS_CLASS_POINTER_DEFINITION(GidEigenIO);

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GidEigenIO( const std::string& rDatafilename,
                GiD_PostMode Mode,
                MultiFileFlag use_multiple_files_flag,
                WriteDeformedMeshFlag write_deformed_flag,
                WriteConditionsFlag write_conditions_flag) :
            GidIO<>(rDatafilename,
                    Mode,
                    use_multiple_files_flag,
                    write_deformed_flag,
                    write_conditions_flag) { }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * Write the post-processed eigensolver-results
    * for Variable of type "double"
    * The label is can be e.g. the Eigenvalue or -frequency
    */
    void WriteEigenResults( ModelPart& rModelPart,
                            const Variable<double>& rVariable,
                            std::string Label,
                            const SizeType NumberOfAnimationStep )
    {
        Label += "_" + rVariable.Name();
        GiD_fBeginResult( mResultFile, (char*)Label.c_str() , "EigenVector_Animation",
                          NumberOfAnimationStep, GiD_Scalar,
                          GiD_OnNodes, NULL, NULL, 0, NULL );

        for (const auto& r_node : rModelPart.Nodes())
        {
            const double& nodal_result = r_node.FastGetSolutionStepValue(rVariable);
            GiD_fWriteScalar( mResultFile, r_node.Id(), nodal_result );
        }

        GiD_fEndResult(mResultFile);
    }

    /**
    * Write the post-processed eigensolver-results
    * for Variable of type "array_1d<double, 3>""
    * The label is can be e.g. the Eigenvalue or -frequency
    */
    void WriteEigenResults( ModelPart& rModelPart,
                            const Variable<array_1d<double, 3>>& rVariable,
                            std::string Label,
                            const SizeType NumberOfAnimationStep)
    {
        Label += "_" + rVariable.Name();
        GiD_fBeginResult( mResultFile, (char*)Label.c_str() , "EigenVector_Animation",
                          NumberOfAnimationStep, GiD_Vector,
                          GiD_OnNodes, NULL, NULL, 0, NULL );

        for (auto& r_node : rModelPart.Nodes())
        {
            const array_1d<double, 3>& nodal_result = r_node.FastGetSolutionStepValue(rVariable);
            GiD_fWriteVector(mResultFile, r_node.Id(), nodal_result[0], nodal_result[1], nodal_result[2]);
        }

        GiD_fEndResult(mResultFile);
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
        std::stringstream buffer;
        buffer << "GidEigenIO" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "GidEigenIO";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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


    ///@}

}; // Class GidEigenIO

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GID_EIGEN_IO_H_INCLUDED  defined

