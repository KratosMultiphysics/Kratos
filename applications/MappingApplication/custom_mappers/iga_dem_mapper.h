//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_IGA_DEM_MAPPER_H_INCLUDED )
#define  KRATOS_IGA_DEM_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"


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

/// Nearest Neighbor Mapper
/** This class implements the Nearest Neighbor Mapping technique.
* Each node on the destination side gets assigned is's closest neighbor on the other side of the interface.
* In the mapping phase every node gets assigned the value of it's neighbor
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/

class IGADEMMapper : public Mapper
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IGADEMMapper
    KRATOS_CLASS_POINTER_DEFINITION(IGADEMMapper);

    ///@}
    ///@name Life Cycle
    ///@{

    IGADEMMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination) : Mapper(
                         rModelPartOrigin, rModelPartDestination) {}

    IGADEMMapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
                          Parameters JsonParameters) : Mapper(
                                  rModelPartOrigin, rModelPartDestination, JsonParameters)
    {
        KRATOS_INFO("IGADEMMapper") << "In Ctor of " << std::endl;
        mpMapperCommunicator->InitializeOrigin(MapperUtilities::Meshless_Point); // IGA
        mpMapperCommunicator->InitializeDestination(MapperUtilities::Element_Center); // DEM
        mpMapperCommunicator->Initialize();


        mpInverseMapper.reset(); // explicitly specified to be safe
    }

    /// Destructor.
    virtual ~IGADEMMapper() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override
    {
        mpMapperCommunicator->UpdateInterface(MappingOptions, SearchRadius);
        if (mpInverseMapper)
        {
            mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);
        }

        if (MappingOptions.Is(MapperFlags::REMESHED))
        {
            ComputeNumberOfNodesAndConditions();
        }

        KRATOS_INFO("IGADEMMapper") << "In UpdateInterface" << std::endl;
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This is not implemented!" << std::endl;
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This is not implemented!" << std::endl;
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This is not implemented!" << std::endl;
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        KRATOS_ERROR << "This is not implemented!" << std::endl;
    }

    Mapper::Pointer Clone(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters) override
    {
        return Kratos::make_shared<IGADEMMapper>(rModelPartOrigin,
                                                          rModelPartDestination,
                                                          JsonParameters);
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
    virtual std::string Info() const override
    {
        return "IGADEMMapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IGADEMMapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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

    Mapper::Pointer mpInverseMapper;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    template <typename T>
    static T GetValueOfNode(InterfaceObject::Pointer pInterfaceObject, //TODO const
                            const Variable< T >& rVariable,
                            const Kratos::Flags& rOptions,
                            const std::vector<double>& rShapeFunctionValues)
    {
        Node<3>* p_base_node = pInterfaceObject->pGetBaseNode();

        return p_base_node->FastGetSolutionStepValue(rVariable);
    }


    template <typename T>
    static void SetValueOfNode(InterfaceObject::Pointer pInterfaceObject,
                               const T& rValue,
                               const Variable< T >& rVariable,
                               const Kratos::Flags& rOptions,
                               const double Factor)
    {
        Node<3>* p_base_node = pInterfaceObject->pGetBaseNode();

        if (rOptions.Is(MapperFlags::ADD_VALUES))
        {
            p_base_node->FastGetSolutionStepValue(rVariable) += rValue * Factor;
        }
        else
        {
            p_base_node->FastGetSolutionStepValue(rVariable) = rValue * Factor;
        }
    }

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
    IGADEMMapper& operator=(IGADEMMapper const& rOther);

    /// Copy constructor.
    //IGADEMMapper(IGADEMMapper const& rOther);

    ///@}

}; // Class IGADEMMapper

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


inline std::istream & operator >>(std::istream& rIStream,
                                  IGADEMMapper& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const IGADEMMapper& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_IGA_DEM_MAPPER_H_INCLUDED  defined