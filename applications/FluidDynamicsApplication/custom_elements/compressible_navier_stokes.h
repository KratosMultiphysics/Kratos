
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Elisa Magliozzi
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES)
#define  KRATOS_COMPRESSIBLE_NAVIER_STOKES

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
template< unsigned int TDim, unsigned int BlockSize = TDim+2, unsigned int TNumNodes = TDim + 1 >
class CompressibleNavierStokes : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokes);
    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, BlockSize> U, Un, Unn;
        BoundedMatrix<double, TNumNodes, TDim> f_ext;
        array_1d<double,TNumNodes> r; // At the moment considering all parameters as constant in the domain (mu, nu, etc...)
        array_1d<double, TDim> f_gauss;
        double r_gauss;

        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;

        double bdf0;
        double bdf1;
        double bdf2;
        double h;             // Element size
        double volume;        // In 2D: element area. In 3D: element volume
        double mu;
        double nu;
        double lambda;
        double c_v;
        double gamma;               //gamma

        //double c;               // TO DO : temporarily use for testing
        //double time;               // TO DO: used for manufactured solution
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    CompressibleNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    CompressibleNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~CompressibleNavierStokes() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_shared< CompressibleNavierStokes < TDim,BlockSize, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_shared< CompressibleNavierStokes < TDim,BlockSize, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        constexpr unsigned int MatrixSize = TNumNodes*(BlockSize);

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Gauss point position
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);

        for(unsigned int igauss = 0; igauss<Ncontainer.size2(); igauss++)
        {
            noalias(data.N) = row(Ncontainer, igauss);

            double v_sc = ShockCapturingViscosity(data);
            double k_sc = ShockCapturingConductivity(data);

            ComputeGaussPointRHSContribution(rhs_local, data, v_sc, k_sc);
            ComputeGaussPointLHSContribution(lhs_local, data, v_sc, k_sc);

            // here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }

        rLeftHandSideMatrix  *= data.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);
        /*
        std::cout<<this->Id()<<" "<<rRightHandSideVector<<std::endl;
        for(unsigned int i=0; i<GetGeometry().size(); i++){
            std::cout<<GetGeometry()[i].Id()<<" "<<GetGeometry()[i].Coordinates()<<std::endl;
            //std::cout<<GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM)<<std::endl;
            //std::cout<<GetGeometry()[i].FastGetSolutionStepValue(DENSITY)<<std::endl;
            //std::cout<<GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY)<<std::endl;
        }
        */
        KRATOS_CATCH("Error in Compressible Navier Stokes Element Symbolic")
    }


    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        constexpr unsigned int MatrixSize = TNumNodes*(BlockSize);

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        // Gauss point position
        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);

        // Loop on gauss point
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<Ncontainer.size2(); igauss++)
        {
            noalias(data.N) = row(Ncontainer, igauss);

            double v_sc = ShockCapturingViscosity(data);
            double k_sc = ShockCapturingConductivity(data);

            ComputeGaussPointRHSContribution(rhs_local, data,v_sc,k_sc);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rRightHandSideVector) += rhs_local;
        }

        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("")

    }

    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Check that all required variables have been registered
        if(MOMENTUM.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"MOMENTUM Key is 0. Check if the application was correctly registered.","");
        if(TOTAL_ENERGY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"TOTAL_ENERGY Key is 0. Check if the application was correctly registered.","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
        if(DYNAMIC_VISCOSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_VISCOSITY Key is 0. Check if the application was correctly registered.","");
        if(KINEMATIC_VISCOSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"KINEMATIC_VISCOSITY Key is 0. Check if the application was correctly registered.","");
        if(CONDUCTIVITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONDUCTIVITY Key is 0. Check if the application was correctly registered.","");
        if(SPECIFIC_HEAT.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"SPECIFIC_HEAT Key is 0. Check if the application was correctly registered.","");
        if(HEAT_CAPACITY_RATIO.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"HEAT_CAPACITY_RATIO Key is 0. Check if the application was correctly registered.","");



        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing MOMENTUM variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing TOTAL_ENERGY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) == false ||
               this->GetGeometry()[i].HasDofFor(MOMENTUM_Y) == false ||
               this->GetGeometry()[i].HasDofFor(MOMENTUM_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing MOMENTUM component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing TOTAL_ENERGY component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(DENSITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing DENSITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

       return 0;

        KRATOS_CATCH("");
    }


    void Calculate(const Variable<double>& rVariable,
                           double& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        if (rVariable == ERROR_RATIO)
        {
//             rOutput = this->SubscaleErrorEstimate(data);
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
        return "CompressibleNavierStokes #";
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

    // Symbolic function implementing the element
    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    double ShockCapturingViscosity(const ElementDataStruct& data);
    double ShockCapturingConductivity(const ElementDataStruct& data);
    void ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(BlockSize),TNumNodes*(BlockSize)>& lhs, const ElementDataStruct& data, double v_sc, double k_sc);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(BlockSize)>& rhs, const ElementDataStruct& data,double v_sc, double k_sc);

    double SubscaleErrorEstimate(const ElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    CompressibleNavierStokes() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    // Auxiliar function to fill the element data structure
    void FillElementData(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    {
        // Getting data for the given geometry
        // double Volume; // In 2D cases Volume variable contains the element area
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), rData.DN_DX, rData.N, rData.volume);

        // Compute element size
        rData.h = ComputeH(rData.DN_DX);

        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        rData.bdf0 = BDFVector[0];
        rData.bdf1 = BDFVector[1];
        rData.bdf2 = BDFVector[2];
        //rData.time = rCurrentProcessInfo[TIME];

        Properties& r_properties = this->GetProperties();
        rData.nu = r_properties.GetValue(KINEMATIC_VISCOSITY);
        rData.mu =  r_properties.GetValue(DYNAMIC_VISCOSITY);
        rData.lambda = r_properties.GetValue(CONDUCTIVITY);
        rData.c_v = r_properties.GetValue(SPECIFIC_HEAT);
        rData.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& moment = this->GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double,3>& moment_n = this->GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM,1);
            const array_1d<double,3>& moment_nn = this->GetGeometry()[i].FastGetSolutionStepValue(MOMENTUM,2);

            for(unsigned int k=0; k<TDim; k++)
            {
                rData.U(i,k+1)   = moment[k];
                rData.Un(i,k+1)  = moment_n[k];
                rData.Unn(i,k+1) = moment_nn[k];
                rData.f_ext(i,k)   = body_force[k];
            }
            rData.U(i,0)= this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            rData.Un(i,0)= this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY,1);
            rData.Unn(i,0)= this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY,2);

            rData.U(i,TDim+1) = this->GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY);
            rData.Un(i,TDim+1) = this->GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY,1);
            rData.Unn(i,TDim+1) = this->GetGeometry()[i].FastGetSolutionStepValue(TOTAL_ENERGY,2);

            rData.r(i) = this->GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);
         }

    }

    //~ template< unsigned int TDim, unsigned int TNumNodes=TDim+1>
    double ComputeH(BoundedMatrix<double,TNumNodes, TDim>& DN_DX)
    {
        double h=0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h += 1.0/h_inv;
        }
        h = sqrt(h)/static_cast<double>(TNumNodes);
        return h;
    }

    // 3D tetrahedra shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,4,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    // 2D triangle shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    }

    // 3D tetrahedra shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(BoundedMatrix<double,1,4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.25; Ncontainer(0,1) = 0.25; Ncontainer(0,2) = 0.25; Ncontainer(0,3) = 0.25;
    }

    // 2D triangle shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(BoundedMatrix<double,1,3>& Ncontainer)
    {
        Ncontainer(0,0) = 1.0/3.0; Ncontainer(0,1) = 1.0/3.0; Ncontainer(0,2) = 1.0/3.0;
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

#endif // KRATOS_STOKES_ELEMENT_SYMBOLIC_3D_INCLUDED  defined
