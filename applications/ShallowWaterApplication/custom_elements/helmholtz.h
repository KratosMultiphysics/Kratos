//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_HELMHOLTZ)
#define  KRATOS_HELMHOLTZ

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
// #include "includes/cfd_variables.h"     // TODO: all the includes, like serializer, are need ??

// Application includes
#include "shallow_water_application_variables.h"

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
/**this element is a 2D Helmholtz element
* formulation is described in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwZ2Zxd09YUmlPZ28/view?usp=sharing
* symbolic implementation is defined in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwaXRKRUpDbmx4VXM/view?usp=sharing
*/
template< unsigned int TDim, unsigned int TNumVar = 1, unsigned int TNumNodes = TDim + 1 >  // TODO: TNumVar is OK? (See helmholtz_cpp_template.cpp)
class Helmholtz : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(Helmholtz);

    struct ElementDataStruct
    {
        //~ bounded_matrix<double, TNumNodes, TDim> v, vn, vnn, vmesh, f;   // TODO: ??
        bounded_matrix<double, TNumNodes, TDim> f;
        array_1d<double,TNumNodes> eta, etan, etann, H, h;              // Elevation, bathymetry and depth. n denotes previous time step.

        bounded_matrix<double, TNumNodes, TDim > DN_DX;   // Shape function spatial derivatives
        array_1d<double, TNumNodes > N;                   // Shape functions

        Matrix C;         // TODO: ??
//        Vector stress;

        double bdf0;
        double bdf1;
        double bdf2;
        double g;             // Gravity
        double l;             // Element size
        double area;          // In 2D: element area. In 3D: element volume
        double dt;            // Time increment
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    Helmholtz(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    Helmholtz(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~Helmholtz() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< Helmholtz < TDim, TNumNodes > >(NewId, this->GetGeometry().Create(rThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return boost::make_shared< Helmholtz < TDim, TNumNodes > >(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*TNumVar;

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        // Loop on gauss points
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

        // Gauss point position
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
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

        rLeftHandSideMatrix  *= data.area/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.area/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Stokes Element Symbolic")
    }


    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        constexpr unsigned int MatrixSize = TNumNodes*TNumVar;

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        // Gauss point position
        bounded_matrix<double,TNumNodes, TNumNodes> Ncontainer;
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

        // rRightHandSideVector *= area/static_cast<double>(TNumNodes);
        rRightHandSideVector *= data.area/static_cast<double>(TNumNodes);

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
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        // Perform basic element checks
        int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
        if(ErrorCode != 0) return ErrorCode;

        // Check that all required variables have been registered
        if(ELEVATION.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"ELEVATION Key is 0. Check if the application was correctly registered.",""); 
        if(BATHYMETRY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"BATHYMETRY Key is 0. Check if the application was correctly registered.","");
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if the application was correctly registered.","");

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(ELEVATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing ELEVATION variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(ELEVATION) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"Missing ELEVATION component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

//        // Check constitutive law          // TODO: Is constitutive law needed?
//        if(mpConstitutiveLaw == nullptr)
//            KRATOS_ERROR << "The constitutive law was not set. Cannot proceed. Call the navier_stokes.h Initialize() method needs to be called.";
//
//        mpConstitutiveLaw->Check(GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("");
    }


    // TODO: Check this Calculate function
    virtual void Calculate(const Variable<double>& rVariable,
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
        // if(rVariable == HEAT_FLUX) //compute the heat flux per unit volume induced by the shearing
        // {
        //     const unsigned int strain_size = (TDim*3)-3;
        //
        //     //struct to pass around the data
        //     ElementDataStruct data;
        //
        //     //getting data for the given geometry
        //     double Volume;
        //     GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, Volume);
        //
        //     //~ for (unsigned int i = 0; i < NumNodes; i++)
        //     for (unsigned int i = 0; i < TNumNodes; i++)
        //     {
        //         const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        //
        //         //~ for(unsigned int k=0; k<Dim; k++)
        //         for(unsigned int k=0; k<TDim; k++)
        //         {
        //             data.v(i,k)   = vel[k];
        //         }
        //     }
        //
        //     if (data.stress.size() != strain_size) data.stress.resize(strain_size,false);
        //
        //     //compute strain
        //     Vector strain(strain_size);
        //     ComputeStrain(data, strain_size, strain);
        //
        //     //create constitutive law parameters:
        //     ConstitutiveLaw::Parameters Values(this->GetGeometry(),GetProperties(),rCurrentProcessInfo);
        //
        //     //set constitutive law flags:
        //     Flags& ConstitutiveLawOptions = Values.GetOptions();
        //     ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        //     ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        //
        //     //this is to pass the shape functions. Unfortunately it is needed to make a copy to a flexible size vector
        //     const Vector Nvec(data.N);
        //     Values.SetShapeFunctionsValues(Nvec);
        //
        //     Values.SetStrainVector(strain); //this is the input parameter
        //     Values.SetStressVector(data.stress); //this is an ouput parameter
        //     // Values.SetConstitutiveMatrix(data.C);      //this is an ouput parameter
        //
        //     //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
        //     //this is ok under the hypothesis that no history dependent behaviour is employed
        //     mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
        //
        //     Output = inner_prod(data.stress, strain);
        // }

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
    virtual std::string Info() const override
    {
        return "Helmholtz #";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
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

    void ComputeGaussPointLHSContribution(bounded_matrix<double,TNumNodes*TNumVar,TNumNodes*TNumVar>& lhs, const ElementDataStruct& data);
    void ComputeGaussPointRHSContribution(array_1d<double,TNumNodes*TNumVar>& rhs, const ElementDataStruct& data);

    double SubscaleErrorEstimate(const ElementDataStruct& data);

    ///@}
    ///@name Protected Operators
    ///@{

    Helmholtz() : Element()
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
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), rData.DN_DX, rData.N, rData.area);

        // Compute element size
        rData.l = ComputeL(rData.DN_DX);

        // Database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        rData.bdf0 = BDFVector[0];
        rData.bdf1 = BDFVector[1];
        rData.bdf2 = BDFVector[2];

//        rData.dyn_tau = rCurrentProcessInfo[DYNAMIC_TAU];  // Only, needed if the temporal dependent term is considered in the subscales
        rData.dt = rCurrentProcessInfo[DELTA_TIME];         // Only, needed if the temporal dependent term is considered in the subscales

//        rData.c = rCurrentProcessInfo[SOUND_VELOCITY];           // Wave velocity

        for (unsigned int i = 0; i < TNumNodes; i++)
        {

            //~ const array_1d<double,3>& body_force = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            //~ const array_1d<double,3>& vel = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            //~ const array_1d<double,3>& vel_n = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            //~ const array_1d<double,3>& vel_nn = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);
            //~ const array_1d<double,3>& vel_mesh = this->GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
//~ 
            //~ for(unsigned int k=0; k<TDim; k++)
            //~ {
                //~ rData.v(i,k)   = vel[k];
                //~ rData.vn(i,k)  = vel_n[k];
                //~ rData.vnn(i,k) = vel_nn[k];
                //~ rData.vmesh(i,k) = vel_mesh[k];
                //~ rData.f(i,k)   = body_force[k];
            //~ }

            rData.eta[i] = this->GetGeometry()[i].FastGetSolutionStepValue(ELEVATION);
            rData.etan[i] = this->GetGeometry()[i].FastGetSolutionStepValue(ELEVATION,1);
            rData.etann[i] = this->GetGeometry()[i].FastGetSolutionStepValue(ELEVATION,2);
            rData.H[i] = this->GetGeometry()[i].FastGetSolutionStepValue(BATHYMETRY);
            rData.h[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DEPTH);
        }

    }

    // Compute element size
    double ComputeL(boost::numeric::ublas::bounded_matrix<double,TNumNodes, TDim>& DN_DX)
    {
        double l=0.0;
        for(unsigned int i=0; i<TNumNodes; i++)
        {
            double l_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                l_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            l += 1.0/l_inv;
        }
        l = sqrt(l)/static_cast<double>(TNumNodes);
        return l;
    }


    // 2D triangle shape functions values at Gauss points
    void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    }


    // 2D triangle shape functions values at centered Gauss point
    void GetShapeFunctionsOnUniqueGauss(boost::numeric::ublas::bounded_matrix<double,1,3>& Ncontainer)
    {
        Ncontainer(0,0) = 1.0/3.0; Ncontainer(0,1) = 1.0/3.0; Ncontainer(0,2) = 1.0/3.0;
    }

 
    //~ // Call the constitutive law to get the stress value
    //~ virtual void ComputeConstitutiveResponse(ElementDataStruct& rData, const ProcessInfo& rCurrentProcessInfo)
    //~ {
        //~ const unsigned int strain_size = (TDim*3)-3;
//~ 
        //~ if(rData.C.size1() != strain_size)
            //~ rData.C.resize(strain_size,strain_size,false);
        //~ if(rData.stress.size() != strain_size)
            //~ rData.stress.resize(strain_size,false);
//~ 
        //~ Vector strain(strain_size);
        //~ ComputeStrain(rData, strain_size, strain);
//~ 
        //~ // Create constitutive law parameters:
        //~ ConstitutiveLaw::Parameters Values(this->GetGeometry(), GetProperties(), rCurrentProcessInfo);
//~ 
        //~ const Vector Nvec(rData.N);
        //~ Values.SetShapeFunctionsValues(Nvec);
//~ 
        //~ // Set constitutive law flags:
        //~ Flags& ConstitutiveLawOptions=Values.GetOptions();
        //~ ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        //~ ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
//~ 
        //~ Values.SetStrainVector(strain);             //this is the input parameter
        //~ Values.SetStressVector(rData.stress);       //this is an ouput parameter
        //~ Values.SetConstitutiveMatrix(rData.C);      //this is an ouput parameter
//~ 
        //~ //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
        //~ //this is ok under the hypothesis that no history dependent behaviour is employed
        //~ mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);
//~ 
    //~ }


    // Compute total depth as h = eta + H
    virtual double ComputeTotalDepth(const double& rEta, const double& rH)
    {
        return rEta + rH;
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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    virtual void load(Serializer& rSerializer) override
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

///@}


} // namespace Kratos.

#endif // KRATOS_HELMHOLTZ_ELEMENT_SYMBOLIC_3D_INCLUDED  defined
