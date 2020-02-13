
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla (based on Elisa Magliozzi previous work)
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED)
#define  KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

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

// TODO: UPDATE THIS INFORMATION
/**this element is a 3D stokes element, stabilized by employing an ASGS stabilization
* formulation is described in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwZ2Zxd09YUmlPZ28/view?usp=sharing
* symbolic implementation is defined in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwaXRKRUpDbmx4VXM/view?usp=sharing
*/
template< unsigned int TDim, unsigned int TBlockSize = TDim + 2, unsigned int TNumNodes = TDim + 1 >
class CompressibleNavierStokesExplicit : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressibleNavierStokesExplicit);

    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, TBlockSize> U;
        BoundedMatrix<double, TNumNodes, TDim> f_ext;
        array_1d<double, TNumNodes> r; // At the moment considering all parameters as constant in the domain (mu, nu, etc...)
        array_1d<double, TDim> f_gauss;
        double r_gauss;

        array_1d<double, TNumNodes > N;
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;

        double h;           // Element size
        double volume;      // In 2D: element area. In 3D: element volume
        double mu;          // Dynamic viscosity
        double nu;          // Kinematic viscosity
        double nu_sc;       // Kinematic viscosity (shock capturing)
        double lambda;      // Heat conductivity
        double lambda_sc;   // Heat conductivity (shock capturing)
        double c_v;         // Heat capacity at constant volume
        double gamma;       // Heat capacity ratio
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    CompressibleNavierStokesExplicit(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    CompressibleNavierStokesExplicit(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~CompressibleNavierStokesExplicit() override {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< CompressibleNavierStokesExplicit < TDim, TBlockSize, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< CompressibleNavierStokesExplicit < TDim, TBlockSize, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling the CalculateLocalSystem() method for the explicit compressible Navier-Stokes element.";

        KRATOS_CATCH("")
    }

    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo) override;

    void AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo) override;

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        if (rVariable == ERROR_RATIO) {
            // rOutput = this->SubscaleErrorEstimate(data);
            this->SetValue(ERROR_RATIO, rOutput);
        }

        KRATOS_CATCH("")
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
        return "CompressibleNavierStokesExplicit #";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    // virtual void PrintData(std::ostream& rOStream) const override

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static member variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    void GetDofList(
        DofsVectorType &ElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const override;

    // double ShockCapturingViscosity(const ElementDataStruct& rData) const;

    // double ShockCapturingConductivity(const ElementDataStruct& rData) const;

    void CalculateShockCapturingValues(ElementDataStruct &rData) const;

    void ComputeGaussPointRHSContribution(
        array_1d<double, TNumNodes * TBlockSize>& rhs,
        const ElementDataStruct& rData);

    double SubscaleErrorEstimate(const ElementDataStruct& rData);

    ///@}
    ///@name Protected Operators
    ///@{

    CompressibleNavierStokesExplicit() = default;

    ///@}
    ///@name Protected Operations
    ///@{

    // Auxiliar function to fill the element data structure
    void FillElementData(
        ElementDataStruct& rData,
        const ProcessInfo& rCurrentProcessInfo);

    double ComputeH(BoundedMatrix<double,TNumNodes, TDim>& DN_DX);

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
                                    Fluid2DASGS& rThis);
 */
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
                                    const Fluid2DASGS& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);
      return rOStream;
    }*/
///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_H_INCLUDED  defined
