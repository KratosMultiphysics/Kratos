//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//                   Jordi Cotela
//

#if !defined(KRATOS_TWO_FLUID_VMS_LINEARIZED_DARCY_H_INCLUDED )
#define  KRATOS_TWO_FLUID_VMS_LINEARIZED_DARCY_H_INCLUDED
// System includes
#include <string>
#include <iostream>

#include "custom_elements/two_fluid_vms.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Version of TwoFluidVMS that treats the Darcy term explicitly.
/** It should work as the TwoFluidVMS element, but using the old
 *  velocity to linearize the non-linear Darcy losses. This
 *  improves solution quality when using a single non-linear iteration.
 */
template< unsigned int TDim,
          unsigned int TNumNodes = TDim + 1 >
class TwoFluidVMSLinearizedDarcy : public TwoFluidVMS<TDim,TNumNodes> {
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of TwoFluidVMS
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidVMSLinearizedDarcy);

    using ElementBaseType = TwoFluidVMS<TDim,TNumNodes>;

    ///definition of node type (default is: Node<3>)
    typedef Node < 3 > NodeType;
    
    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    ///definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    ///definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    TwoFluidVMSLinearizedDarcy(IndexType NewId = 0) : ElementBaseType(NewId) {}

    ///Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    TwoFluidVMSLinearizedDarcy(IndexType NewId, const NodesArrayType& ThisNodes)
        : ElementBaseType(NewId, ThisNodes) {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    TwoFluidVMSLinearizedDarcy(IndexType NewId, GeometryType::Pointer pGeometry)
        : ElementBaseType(NewId, pGeometry) {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    TwoFluidVMSLinearizedDarcy(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : ElementBaseType(NewId, pGeometry, pProperties) {}

    /// Destructor.
    ~TwoFluidVMSLinearizedDarcy() override {}

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
     * Returns a pointer to a new TwoFluidVMSLinearizedDarcy element, created using given input
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override {

        return Kratos::make_shared<TwoFluidVMSLinearizedDarcy>(
            NewId, (this->GetGeometry()).Create(ThisNodes), pProperties);
    }

    /// Create a new element of this type
    /**
     * Returns a pointer to a new TwoFluidVMSLinearizedDarcy element, created using given input
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override {

        return Kratos::make_shared<TwoFluidVMSLinearizedDarcy>(
            NewId, pGeom, pProperties);
    }

    ///@}

protected:

    ///@name Protected Operations
    ///@{

    double CalculateDarcyTerm(
        const double Density,
        const double DynamicViscosity,
        const double LinearCoefficient,
        const double NonlinearCoefficient,
        const array_1d<double, TNumNodes>& rShapefunctions) override {

        array_1d<double,3> old_velocity;
        this->GetAdvectiveVel(old_velocity, rShapefunctions,1);
        const double old_velocity_norm = MathUtils<double>::Norm3(old_velocity);

        return DynamicViscosity * LinearCoefficient + Density * NonlinearCoefficient*old_velocity_norm;
    }
    ///@}

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ElementBaseType);
    }

    void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ElementBaseType);
    }

    ///@}
};

///@}
///@name Input and output
///@{

/// input stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::istream & operator >>(std::istream& rIStream,
                                  TwoFluidVMSLinearizedDarcy<TDim, TNumNodes>& rThis)
{
    return rIStream;
}
/// output stream function
template< unsigned int TDim,
          unsigned int TNumNodes >
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const TwoFluidVMSLinearizedDarcy<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif //KRATOS_TWO_FLUID_VMS_LINEARIZED_DARCY_H_INCLUDED