//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

#if !defined(KRATOS_METIS_DIVIDE_SUBMODELPARTS_HETEROGENEOUS_H_INCLUDED )
#define  KRATOS_METIS_DIVIDE_SUBMODELPARTS_HETEROGENEOUS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <unordered_set>

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

// Application includes
#include "metis_divide_heterogeneous_input_process.h"


namespace Kratos {

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

/// Short class definition.
class MetisDivideSubModelPartsHeterogeneousInputProcess : public MetisDivideHeterogeneousInputProcess {
public:

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Pointer Definitions
  /// Pointer definition of MetisDivideSubModelPartsHeterogeneousInputProcess
  KRATOS_CLASS_POINTER_DEFINITION(MetisDivideSubModelPartsHeterogeneousInputProcess);

  ///@}
  ///@name Life Cycle
  ///@{

  MetisDivideSubModelPartsHeterogeneousInputProcess(IO& rIO, Parameters Settings, SizeType NumberOfPartitions, int Dimension = 3, int Verbosity = 0, bool SynchronizeConditions = false)
    : MetisDivideHeterogeneousInputProcess(rIO, NumberOfPartitions, Dimension, Verbosity, SynchronizeConditions),
      mSettings(Settings) {}

  /// Destructor.
  ~MetisDivideSubModelPartsHeterogeneousInputProcess() override {
  }

  ///@}
  ///@name Operators
  ///@{

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
  std::string Info() const override {
    std::stringstream buffer;
    buffer << "MetisDivideSubModelPartsHeterogeneousInputProcess" ;
    return buffer.str();

  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override {
      rOStream << "MetisDivideSubModelPartsHeterogeneousInputProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override {}

  ///@}
  ///@name Friends
  ///@{

  ///@}

protected:

  ///@name Protected static Member Variables
  ///@{
    Parameters mSettings;

  ///@}
  ///@name Protected member Variables
  ///@{

  ///@}
  ///@name Protected Operators
  ///@{

  ///@}
  ///@name Protected Operations
  ///@{

    void GetNodesPartitions(std::vector<idxtype> &rNodePartition, SizeType &rNumNodes) override;

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

  /// Assignment operator.
  MetisDivideSubModelPartsHeterogeneousInputProcess& operator=(MetisDivideSubModelPartsHeterogeneousInputProcess const& rOther);
  ///@}

}; // Class MetisDivideSubModelPartsHeterogeneousInputProcess


///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif //KRATOS_METIS_DIVIDE_SUBMODELPARTS_HETEROGENEOUS_H_INCLUDED  defined
