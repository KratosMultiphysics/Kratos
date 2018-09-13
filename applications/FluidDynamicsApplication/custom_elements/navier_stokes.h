//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_NAVIER_STOKES)
#define  KRATOS_NAVIER_STOKES

// System includes

// External includes
#include "boost/smart_ptr.hpp"

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
template< unsigned int TDim, unsigned int TNumNodes = TDim + 1 >
class NavierStokes : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(NavierStokes);

    struct ElementDataStruct
    {
        BoundedMatrix<double, TNumNodes, TDim> v, vn, vnn, vmesh, f;
        array_1d<double,TNumNodes> p, pn, pnn, rho, mu;

        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;

        Matrix C;
        Vector stress;
        Vector strain;

        double bdf0;
        double bdf1;
        double bdf2;
        double c;             // Wave velocity (needed if artificial compressibility is considered)
        double h;             // Element size
        double volume;        // In 2D: element area. In 3D: element volume
        double dt;            // Time increment
        double dyn_tau;       // Dynamic tau considered in ASGS stabilization coefficients
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    NavierStokes(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    NavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~NavierStokes() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_shared< NavierStokes < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_shared< NavierStokes < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

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

            ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);
            ComputeGaussPointLHSContribution(lhs_local, data);

            // here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rLeftHandSideMatrix) += lhs_local;
            noalias(rRightHandSideVector) += rhs_local;
        }

        rLeftHandSideMatrix  *= data.volume/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Stokes Element Symbolic")
    }


    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

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

            ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
            noalias(rRightHandSideVector) += rhs_local;
        }

        // rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
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
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
        if(PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
        if(DYNAMIC_TAU.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check if the application was correctly registered.","");
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if the application was correctly registered.","");
        if(SOUND_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"SOUND_VELOCITY Key is 0. Check if the application was correctly registered.","");

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

        // Check constitutive law
        if(mpConstitutiveLaw == nullptr)
            KRATOS_ERROR << "The constitutive law was not set. Cannot proceed. Call the navier_stokes.h Initialize() method needs to be called.";

        mpConstitutiveLaw->Check(GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

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
            rOutput = this->SubscaleErrorEstimate(data);
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
        return "NavierStokes #";
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

    // Constitutive law pointer
    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

    // Symbolic function implementing the element
    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override;
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void ComputeGaussPointLHSContribution(BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs, const ElementDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs, const ElementDataStruct& data);

    double SubscaleErrorEstimate(const ElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    NavierStokes() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

    // Element initialization (constitutive law)
    void Initialize() override
    {
        KRATOS_TRY

        mpConstitutiveLaw = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mpConstitutiveLaw->InitializeMaterial( GetProperties(), this->GetGeometry(), row( this->GetGeometry().ShapeFunctionsValues(), 0 ) );

        KRATOS_CATCH( "" )
    }

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

        rData.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];  // Only, needed if the temporal dependent term is considered in the subscales
        rData.dt = rCurrentProcessInfo[DELTA_TIME];         // Only, needed if the temporal dependent term is considered in the subscales

        rData.c = rCurrentProcessInfo[SOUND_VELOCITY];           // Wave velocity

        for (unsigned int i = 0; i < TNumNodes; i++)
        {

            const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& vel_n = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);
            const array_1d<double,3>& vel_mesh = this->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);

            for(unsigned int k=0; k<TDim; k++)
            {
                rData.v(i,k)   = vel[k];
                rData.vn(i,k)  = vel_n[k];
                rData.vnn(i,k) = vel_nn[k];
                rData.vmesh(i,k) = vel_mesh[k];
                rData.f(i,k)   = body_force[k];
            }

            rData.p[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            rData.pn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,1);
            rData.pnn[i] = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,2);
            rData.rho[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
            rData.mu[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
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

    // Computes the strain rate in Voigt notation
    void ComputeStrain(ElementDataStruct& rData, const unsigned int& strain_size)
    {
        const BoundedMatrix<double, TNumNodes, TDim>& v = rData.v;
        const BoundedMatrix<double, TNumNodes, TDim>& DN = rData.DN_DX;

        // Compute strain (B*v)
        // 3D strain computation
        if (strain_size == 6)
        {
            rData.strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
            rData.strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
            rData.strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
            rData.strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
            rData.strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
            rData.strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
        }
        // 2D strain computation
        else if (strain_size == 3)
        {
            rData.strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
            rData.strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
            rData.strain[2] = DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
        }
    }

    // Call the constitutive law to get the stress value
    virtual void ComputeConstitutiveResponse(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int strain_size = (TDim*3)-3;

        if(rData.C.size1() != strain_size)
            rData.C.resize(strain_size,strain_size,false);
        if(rData.stress.size() != strain_size)
            rData.stress.resize(strain_size,false);
        if(rData.strain.size() != strain_size)
            rData.strain.resize(strain_size,false);

        ComputeStrain(rData, strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(this->GetGeometry(), GetProperties(), rCurrentProcessInfo);

        const Vector Nvec(rData.N);
        Values.SetShapeFunctionsValues(Nvec);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        Values.SetStrainVector(rData.strain);       //this is the input parameter
        Values.SetStressVector(rData.stress);       //this is an ouput parameter
        Values.SetConstitutiveMatrix(rData.C);      //this is an ouput parameter

        //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
        //this is ok under the hypothesis that no history dependent behaviour is employed
        mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    }

    virtual double ComputeEffectiveViscosity(const ElementDataStruct& rData)
    {
        // Computes effective viscosity as sigma = mu_eff*eps -> sigma*sigma = mu_eff*sigma*eps -> sigma*sigma = mu_eff*(mu_eff*eps)*eps
        // return sqrt(inner_prod(rStress, rStress)/inner_prod(rStrain, rStrain)); // TODO: Might yield zero in some cases, think a more suitable manner

        // Computes the effective viscosity as the average of the lower diagonal constitutive tensor
        double mu_eff = 0.0;
        const unsigned int strain_size = (TDim-1)*3;
        for (unsigned int i=TDim; i<strain_size; ++i){mu_eff += rData.C(i,i);}
        mu_eff /= (strain_size - TDim);

        return mu_eff;
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
