//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//  Main authors:    Minas Apostolakis

#if !defined(KRATOS_COMBINE_SOLID_SHELL_MODEL_PARTS )
#define  KRATOS_COMBINE_SOLID_SHELL_MODEL_PARTS

// System includes

// External includes

// Project includes
#include "containers/model.h"

#include "processes/process.h"

#include "geometries/coupling_geometry.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/* @class CombineSolidShellModelPartsProcess
 * @ingroup IgaApplication
 * @brief This class outputs the location of the quadrature points within the local space of the containing geometry. */
class KRATOS_API(IGA_APPLICATION) CombineSolidShellModelPartsProcess 
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CombineSolidShellModelPartsProcess
    KRATOS_CLASS_POINTER_DEFINITION(CombineSolidShellModelPartsProcess);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

   CombineSolidShellModelPartsProcess(Model& mModel) : 
       Process(),
       _Model(mModel){};


    /// Destructor.
    ~CombineSolidShellModelPartsProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    // Minas Declare here the operations you want to use

    void ExecuteInitialize();
    void RecursiveAddEntities(ModelPart& TwinModelPart, ModelPart& ModelPart) ;
    void ReorderGeometryIds(size_t& CoupleIdCounter, ModelPart& mModelPart);
    void ReorderNodeIds(size_t& CoupleIdCounter, ModelPart& mModelPart);
    void ReorderElementIds(size_t& CoupleIdCounter, ModelPart& mModelPart);
    void ReorderConditionIds(size_t& CoupleIdCounter, ModelPart& mModelPart);
    void ReorderMasterSlaveConstraintIds(size_t& CoupleIdCounter, ModelPart& mModelPart);
    void ReorderPropertyIds(size_t& CoupleIdCounter, ModelPart& mModelPart);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CombineSolidShellModelPartsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CombineSolidShellModelPartsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    ///@name Member Variables
    ///@{

    //Model part and different settings
    Model& _Model;

    ///@}

}; // Class CombineSolidShellModelPartsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  CombineSolidShellModelPartsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CombineSolidShellModelPartsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_OUTPUT_QUADRATURE_DOMAIN_PROCESS_H_INCLUDED  defined
