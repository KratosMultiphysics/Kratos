//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//



#if !defined(KRATOS_STOKES_ELEMENT_SYMBOLIC_3D_INCLUDED )
#define  KRATOS_STOKES_ELEMENT_SYMBOLIC_3D_INCLUDED


// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"

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

/**this element is a 3D stokes element, stabilized by employing an ASGS stabilization
* formulation is described in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwZ2Zxd09YUmlPZ28/view?usp=sharing
* symbolic implementation is defined in the file:
*    https://drive.google.com/file/d/0B_gRLnSH5vCwaXRKRUpDbmx4VXM/view?usp=sharing
*/
class Stokes3D
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(Stokes3D);

    template <unsigned int TNumNodes, unsigned int TDim>
    struct element_data
    {
        BoundedMatrix<double,TNumNodes, TDim> v, vn, vnn, f;
        array_1d<double,TNumNodes> p, rho;

        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;

        Matrix C;
        Vector stress;

        double bdf0;
        double bdf1;
        double bdf2;
        double h;
        double dyn_tau_coeff;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    Stokes3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    Stokes3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~Stokes3D() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_shared< Stokes3D >(NewId, GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }
    Element::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< Stokes3D >(NewId, pGeom, pProperties);
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const unsigned int NumNodes = 4;
        const unsigned int Dim = 3;
        const int ndofs = Dim + 1;
        const unsigned int MatrixSize = NumNodes*ndofs;

        if (rLeftHandSideMatrix.size1() != MatrixSize)
            rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        //struct to pass around the data
        element_data<NumNodes,Dim> data;

        //getting data for the given geometry

        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);

        //compute element size
//         data.h = ComputeH<4,3>(data.DN_DX, Volume);

        //gauss point position
        BoundedMatrix<double,NumNodes, NumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);

        //database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        data.bdf0 = BDFVector[0];
        data.bdf1 = BDFVector[1];
        data.bdf2 = BDFVector[2];
        data.dyn_tau_coeff = rCurrentProcessInfo[DYNAMIC_TAU] * data.bdf0;


        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel_n = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);

            for(unsigned int k=0; k<Dim; k++)
            {
                data.v(i,k)   = vel[k];
                data.vn(i,k)  = vel_n[k];
                data.vnn(i,k) = vel_nn[k];
                data.f(i,k)   = body_force[k];
            }

            data.p[i] = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            data.rho[i] = GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
        }

        //allocate memory needed
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        //loop on gauss points
//         noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
//         noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<1; igauss++) //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
        {
//             noalias(data.N) = row(Ncontainer, igauss); //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!

            ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);
            ComputeGaussPointLHSContribution(lhs_local, data);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
            noalias(rLeftHandSideMatrix) = lhs_local; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
            noalias(rRightHandSideVector) = rhs_local; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
        }

        rLeftHandSideMatrix  *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
        rRightHandSideVector *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!

//         rLeftHandSideMatrix  *= Volume/static_cast<double>(NumNodes);
//         rRightHandSideVector *= Volume/static_cast<double>(NumNodes);



        KRATOS_CATCH("Error in Stokes Element Symbolic")
    }





    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const unsigned int NumNodes = 4;
        const unsigned int Dim = 3;
        const int ndofs = Dim + 1;
        const unsigned int MatrixSize = NumNodes*ndofs;

        if (rRightHandSideVector.size() != MatrixSize)
            rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

        //struct to pass around the data
        element_data<NumNodes,Dim> data;

        //getting data for the given geometry
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);

        //gauss point position
        BoundedMatrix<double,NumNodes, NumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);

        //database access to all of the variables needed
        const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
        data.bdf0 = BDFVector[0];
        data.bdf1 = BDFVector[1];
        data.bdf2 = BDFVector[2];
        data.dyn_tau_coeff = rCurrentProcessInfo[DYNAMIC_TAU] * data.bdf0;


        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d<double,3>& vel_n = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
            const array_1d<double,3>& vel_nn = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);

            for(unsigned int k=0; k<Dim; k++)
            {
                data.v(i,k)   = vel[k];
                data.vn(i,k)  = vel_n[k];
                data.vnn(i,k) = vel_nn[k];
                data.f(i,k)   = body_force[k];
            }


            data.p[i] = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
            data.rho[i] = GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
        }

        //allocate memory needed
        array_1d<double,MatrixSize> rhs_local;

        //loop on gauss points - ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<1; igauss++) //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
        {
//             noalias(data.N) = row(Ncontainer, igauss);

            ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
            noalias(rRightHandSideVector) += rhs_local;
        }

        rRightHandSideVector *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//         rRightHandSideVector *= Volume/static_cast<double>(NumNodes);

        KRATOS_CATCH("")

    }





    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const unsigned int NumNodes = 4;
        const int Dim = 3;

        if (rResult.size() != NumNodes*(Dim+1))
            rResult.resize(NumNodes*(Dim+1), false);

        for(unsigned int i=0; i<NumNodes; i++)
        {
            rResult[i*(Dim+1)  ]  =  GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
            rResult[i*(Dim+1)+1]  =  GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
            rResult[i*(Dim+1)+2]  =  GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
            rResult[i*(Dim+1)+3]  =  GetGeometry()[i].GetDof(PRESSURE).EquationId();
        }


        KRATOS_CATCH("")

    }





    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const unsigned int NumNodes = 4;
        const int Dim = 3;

        if (ElementalDofList.size() != NumNodes*(Dim+1))
            ElementalDofList.resize(NumNodes*(Dim+1));

        for(unsigned int i=0; i<NumNodes; i++)
        {
            ElementalDofList[i*(Dim+1)  ]  =  GetGeometry()[i].pGetDof(VELOCITY_X);
            ElementalDofList[i*(Dim+1)+1]  =  GetGeometry()[i].pGetDof(VELOCITY_Y);
            ElementalDofList[i*(Dim+1)+2]  =  GetGeometry()[i].pGetDof(VELOCITY_Z);
            ElementalDofList[i*(Dim+1)+3]  =  GetGeometry()[i].pGetDof(PRESSURE);
        }

        KRATOS_CATCH("");
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

        // Checks on nodes

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
        if(mp_constitutive_law == nullptr)
            KRATOS_ERROR << "The constitutive law was not set. Cannot proceed.";

        mp_constitutive_law->Check(GetProperties(),GetGeometry(),rCurrentProcessInfo);

        return 0;

        KRATOS_CATCH("");
    }

    void Calculate(const Variable<double>& rVariable,
                           double& Output,
                           const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if(rVariable == HEAT_FLUX) //compute the heat flux per unit volume induced by the shearing
        {
            const unsigned int NumNodes = 4;
            const unsigned int Dim = 3;
            const unsigned int strain_size = 6;

            //struct to pass around the data
            element_data<NumNodes,Dim> data;

            //getting data for the given geometry
            double Volume;
            GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);


            for (unsigned int i = 0; i < NumNodes; i++)
            {
                const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);

                for(unsigned int k=0; k<Dim; k++)
                {
                    data.v(i,k)   = vel[k];
                }
            }

            if (data.stress.size() != strain_size) data.stress.resize(strain_size,false);

            const BoundedMatrix<double,NumNodes,Dim>& v = data.v;
            const BoundedMatrix<double,NumNodes,Dim>& DN = data.DN_DX;

            //compute strain
            Vector strain(strain_size);
            strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
            strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
            strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
            strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
            strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
            strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);

            //create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

            //set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            //this is to pass the shape functions. Unfortunately it is needed to make a copy to a flexible size vector
            const Vector Nvec(data.N);
            Values.SetShapeFunctionsValues(Nvec);

            Values.SetStrainVector(strain); //this is the input parameter
            Values.SetStressVector(data.stress); //this is an ouput parameter
//             Values.SetConstitutiveMatrix(data.C);      //this is an ouput parameter

            //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
            //this is ok under the hypothesis that no history dependent behaviour is employed
            mp_constitutive_law->CalculateMaterialResponseCauchy(Values);

            Output = inner_prod(data.stress, strain);
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
        return "Stokes3D #";
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;


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

    //this is the symbolic function implementing the element
    void ComputeGaussPointLHSContribution(BoundedMatrix<double,16,16>& lhs, const element_data<4,3>& data);
    void ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data<4,3>& data);

    ///@}
    ///@name Protected Operators
    ///@{

    Stokes3D() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{
    template< unsigned int TNumNodes, unsigned int TDim>
    double ComputeH(BoundedMatrix<double,TNumNodes, TDim>& DN_DX, const double Volume)
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

    //gauss points for the 3D case
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,4, 4>& Ncontainer)
    {
        Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
        Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;
        Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
        Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    }

    //gauss points for the 2D case
    void GetShapeFunctionsOnGauss(BoundedMatrix<double,3,3>& Ncontainer)
    {
        const double one_sixt = 1.0/6.0;
        const double two_third = 2.0/3.0;
        Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third;
        Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt;
        Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt;
    }

    virtual void ComputeConstitutiveResponse(element_data<4,3>& data, ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int nnodes = 4;
        const unsigned int dim = 3;
        const unsigned int strain_size = 6;

        if(data.C.size1() != strain_size)
            data.C.resize(strain_size,strain_size,false);
        if(data.stress.size() != strain_size)
            data.stress.resize(strain_size,false);

        const BoundedMatrix<double,nnodes,dim>& v = data.v;
        const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;


//         noalias(data.C) = ZeroMatrix(strain_size,strain_size);
//         const double nu = GetProperties()[VISCOSITY];
//         const double rho = inner_prod(data.rho, data.N);

        //compute strain
        Vector strain(strain_size);
        strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
        strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
        strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
        strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
        strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
        strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);

        //here we shall call the constitutive law
//         data.C(0,0) = 2.0*nu*rho;
//         data.C(1,1) = 2.0*nu*rho;
//         data.C(2,2) = 2.0*nu*rho;
//         data.C(3,3) = nu*rho;
//         data.C(4,4) = nu*rho;
//         data.C(5,5) = nu*rho;
//
//         noalias(data.stress) = prod(data.C,strain);


        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        const Vector Nvec(data.N);
        Values.SetShapeFunctionsValues(Nvec);
        //set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        Values.SetStrainVector(strain); //this is the input parameter
        Values.SetStressVector(data.stress); //this is an ouput parameter
        Values.SetConstitutiveMatrix(data.C);      //this is an ouput parameter

        //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
        //this is ok under the hypothesis that no history dependent behaviour is employed
        mp_constitutive_law->CalculateMaterialResponseCauchy(Values);

    }


    void Initialize() override
    {
        KRATOS_TRY

        mp_constitutive_law = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mp_constitutive_law->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), 0 ) );

        KRATOS_CATCH( "" )
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


protected:
    ConstitutiveLaw::Pointer mp_constitutive_law = nullptr;

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


