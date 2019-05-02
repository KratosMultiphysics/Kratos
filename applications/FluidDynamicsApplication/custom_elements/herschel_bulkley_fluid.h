#ifndef KRATOS_HERSCHEL_BULKLEY_FLUID_H
#define KRATOS_HERSCHEL_BULKLEY_FLUID_H

#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"

namespace Kratos {

///@addtogroup FluidDynamicsApplication
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

template< class TBaseElement,
          class TShapeFunctionValues = typename TBaseElement::ShapeFunctionsType,
          class TShapeFunctionGradients = typename TBaseElement::ShapeFunctionDerivativesType >
class HerschelBulkleyFluid : public TBaseElement
{
public:

    ///@name Type Definitions
    ///@{

    // Pointer types for HerschelBulkleyFluid
    KRATOS_CLASS_POINTER_DEFINITION(HerschelBulkleyFluid);

    /// Node type (default is: Node<3>)
    typedef Node <3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Properties PropertiesType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    HerschelBulkleyFluid(IndexType NewId = 0) :
        TBaseElement(NewId)
    {}

    /// Constructor using an array of nodes.
    /**
     * @param NewId Index of the new element
     * @param ThisNodes An array containing the nodes of the new element
     */
    HerschelBulkleyFluid(IndexType NewId, const NodesArrayType& ThisNodes) :
        TBaseElement(NewId, ThisNodes)
    {}

    /// Constructor using a geometry object.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     */
    HerschelBulkleyFluid(IndexType NewId, GeometryType::Pointer pGeometry) :
        TBaseElement(NewId, pGeometry)
    {}

    /// Constuctor using geometry and properties.
    /**
     * @param NewId Index of the new element
     * @param pGeometry Pointer to a geometry object
     * @param pProperties Pointer to the element's properties
     */
    HerschelBulkleyFluid(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
        TBaseElement(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~HerschelBulkleyFluid() override
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type.
    /**
     * Returns a pointer to a new element, created using given input
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared<HerschelBulkleyFluid<TBaseElement> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new element of this type.
	/**
	 @param NewId Index of the new element
     @param pGeom A pointer to the geometry of the new element
	 @param pProperties Pointer to the element's properties
	 */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< HerschelBulkleyFluid<TBaseElement> >(NewId,pGeom,pProperties);
    }


    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int Error = 0;

        // Check that any required model parameters are defined


        // Call the underlying element's check routine
        Error = TBaseElement::Check(rCurrentProcessInfo);

        return Error;

        KRATOS_CATCH("");
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
        buffer << "HerschelBulkleyFluid " ;
        buffer << TBaseElement::Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HerschelBulkleyFluid ";
        TBaseElement::PrintInfo(rOStream);
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


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief EffectiveViscosity Calculate the effective viscosity at given integration point using the Herschel-Bulkley constitutive model.
     *
     *
     * @note This assumes that the underlying fluid element provides EquivalentStrainRate(rDN_DX) and
     * EvaluateInPoint(value,VARIABLE,rN) methods.
     *
     * @param Density The fluid's density at the integration point (unused for Herschel-Bulkley fluids).
     * @param rN Nodal shape functions evaluated at the integration points (area coordinates for the point).
     * @param rDN_DX Shape function derivatives at the integration point.
     * @param ElemSize Representative length of the element (unused for Herschel-Bulkley fluids).
     * @param rProcessInfo ProcessInfo instance passed from the ModelPart, containing additional data
     * @return The effective viscosity, in dynamic units (Pa*s or equivalent).
     */
    double EffectiveViscosity(double Density,
                                      const TShapeFunctionValues &rN,
                                      const TShapeFunctionGradients &rDN_DX,
                                      double ElemSize,
                                      const ProcessInfo &rProcessInfo) override
    {
        double GammaDot = this->EquivalentStrainRate(rDN_DX);

        double YieldStress = rProcessInfo[YIELD_STRESS];
        double m = rProcessInfo[REGULARIZATION_COEFFICIENT];
        double k = rProcessInfo[POWER_LAW_K];
        double n = rProcessInfo[POWER_LAW_N];

        double DynamicViscosity;

        if (GammaDot > 1e-12) // Normal behaviour
        {
            DynamicViscosity = k*std::pow(GammaDot,n-1);
            double Regularization = 1.0 - std::exp(-m*GammaDot);
            DynamicViscosity += Regularization * YieldStress / GammaDot;
        }
        else // fallback to avoid division by zero
        {
            // Fall back to Newtonian behaviour
            DynamicViscosity = k;
        }

        return DynamicViscosity;
    }


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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement );
    }

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
    HerschelBulkleyFluid& operator=(HerschelBulkleyFluid const& rOther){}

    /// Copy constructor.
    HerschelBulkleyFluid(HerschelBulkleyFluid const& rOther){}


    ///@}

  }; // Class ClassName

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TBaseElement >
inline std::istream& operator >> (std::istream& rIStream,
                                  HerschelBulkleyFluid<TBaseElement>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TBaseElement >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const HerschelBulkleyFluid<TBaseElement>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_HERSCHEL_BULKLEY_FLUID_H
