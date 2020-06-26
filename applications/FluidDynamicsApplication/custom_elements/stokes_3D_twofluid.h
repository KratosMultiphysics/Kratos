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



#if !defined(KRATOS_STOKES_ELEMENT_TWOFLUID_3D_INCLUDED )
#define  KRATOS_STOKES_ELEMENT_TWOFLUID_3D_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/cfd_variables.h"
#include "utilities/split_tetrahedra.h"
#include "custom_elements/stokes_3D.h"

#include "utilities/enrichment_utilities_duplicate_dofs.h"

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
class Stokes3DTwoFluid
    : public Stokes3D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Stokes3DTwoFluid);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    Stokes3DTwoFluid(IndexType NewId, GeometryType::Pointer pGeometry)
        : Stokes3D(NewId, pGeometry)
    {}

    Stokes3DTwoFluid(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Stokes3D(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    ~Stokes3DTwoFluid() override {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive< Stokes3DTwoFluid >(NewId, GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    Element::Pointer Create(IndexType NewId,
                           GeometryType::Pointer pGeom,
                           PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< Stokes3DTwoFluid >(NewId, pGeom, pProperties);
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

        array_1d<double, NumNodes> distances;
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
            distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }

        //allocate memory needed
        BoundedMatrix<double,MatrixSize, MatrixSize> lhs_local;
        array_1d<double,MatrixSize> rhs_local;

        unsigned int npos=0, nneg=0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                npos++;
            else
                nneg++;
        }


        //here we decide if the element is all FLUID/AIR/MIXED
        if(npos == NumNodes) //all AIR
        {
            ComputeElementAsAIR<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer, rCurrentProcessInfo);
        }
        else if (nneg == NumNodes) //all FLUID
        {
            ComputeElementAsFLUID<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer, rCurrentProcessInfo);
        }
        else //element includes both FLUID and AIR
        {
            ComputeElementAsMIXED<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo, distances);
        }


        KRATOS_CATCH("Error in StokesTwoFluid Element Symbolic")
    }





    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_THROW_ERROR(std::logic_error,"method not implemented yet","");

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
        if(DISTANCE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DISTANCE Key is 0. Check if the application was correctly registered.","");

        if(PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
        if(DENSITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
        if(DYNAMIC_TAU.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check if the application was correctly registered.","");
        if(DELTA_TIME.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if the application was correctly registered.","");

        // Checks on nodes

        //check Properties
        if(GetProperties().Has(DENSITY_AIR) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY_AIR is not set","");
        if(GetProperties().Has(CONSTITUTIVE_LAW) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW is not set","");

        //check constitutive CONSTITUTIVE_LAW
        GetProperties().GetValue(CONSTITUTIVE_LAW)->Check(GetProperties(),GetGeometry(),rCurrentProcessInfo);

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
                    this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

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
            double distance_center = 0.0;
            for(unsigned int i=0; i<GetGeometry().size(); i++)
                distance_center += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            distance_center/=static_cast<double>(GetGeometry().size());

            if(distance_center > 0) //AIR
            {
                Output=0.0;
            }
            else //OTHER MATERIAL
            {
                const unsigned int NumNodes = 4;
                const unsigned int strain_size = 6;

                //create constitutive law parameters:
                const Properties& r_properties = GetProperties();
                ConstitutiveLaw::Parameters Values(GetGeometry(),r_properties,rCurrentProcessInfo);

                //set constitutive law flags:
                Flags& ConstitutiveLawOptions=Values.GetOptions();
                ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
                ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

                Vector Nvec(NumNodes);
                for(unsigned int i=0; i<NumNodes; i++) Nvec[i]=0.25;
                Values.SetShapeFunctionsValues(Nvec);

                Vector strain = CalculateStrain();
                Values.SetStrainVector(strain);

                Vector stress(strain_size);
                Values.SetStressVector(stress); //this is an ouput parameter

    //             Values.SetConstitutiveMatrix(data.C);      //this is an ouput parameter

                //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
                //this is ok under the hypothesis that no history dependent behaviour is employed
                mp_constitutive_law->CalculateMaterialResponseCauchy(Values);

                Output = inner_prod(stress, strain);
            }
        }
        else if(rVariable == EQ_STRAIN_RATE) //compute effective strain rate
        {
            const unsigned int NumNodes = 4;

            const Properties& r_properties = GetProperties();
            ConstitutiveLaw::Parameters Values(GetGeometry(),r_properties,rCurrentProcessInfo);

            Vector Nvec(NumNodes);
            for(unsigned int i=0; i<NumNodes; i++) Nvec[i]=0.25;
            Values.SetShapeFunctionsValues(Nvec);

            Vector strain = CalculateStrain();
            Values.SetStrainVector(strain);

            Output = mp_constitutive_law->CalculateValue(Values,rVariable, Output);
        } 
        else if(rVariable == EFFECTIVE_VISCOSITY) //compute the heat flux per unit volume induced by the shearing
        {
            double distance_center = 0.0;
            for(unsigned int i=0; i<GetGeometry().size(); i++)
                distance_center += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            distance_center/=static_cast<double>(GetGeometry().size());

            if(distance_center > 0) //AIR
            {
                const Properties& r_properties = GetProperties();
                Output = r_properties[DYNAMIC_VISCOSITY];
            }
            else //OTHER MATERIAL
            {
                const unsigned int NumNodes = 4;

                const Properties& r_properties = GetProperties();
                ConstitutiveLaw::Parameters Values(GetGeometry(),r_properties,rCurrentProcessInfo);

                Vector Nvec(NumNodes);
                for(unsigned int i=0; i<NumNodes; i++) Nvec[i]=0.25;
                Values.SetShapeFunctionsValues(Nvec);

                Vector strain = CalculateStrain();
                Values.SetStrainVector(strain);

                Output = mp_constitutive_law->CalculateValue(Values,rVariable, Output);
            }         
        }

        KRATOS_CATCH("")
    }

    template<int MatrixSize, int NumNodes>
    void ComputeElementAsAIR(BoundedMatrix<double,MatrixSize, MatrixSize>& lhs_local,
                             array_1d<double,MatrixSize>& rhs_local,
                             Matrix& rLeftHandSideMatrix,
                             Vector& rRightHandSideVector,
                             const double& Volume,
                             element_data<4,3>& data,
                             BoundedMatrix<double,NumNodes, NumNodes>& Ncontainer,
                             ProcessInfo& rCurrentProcessInfo)
    {
        const double air_density = GetProperties()[DENSITY_AIR];
        for (unsigned int i = 0; i < NumNodes; i++)
            data.rho[i] = air_density;
        const double air_nu = GetProperties()[DYNAMIC_VISCOSITY]; //ATTENTION: not using here the real visosity of air

        const double weight = Volume/static_cast<double>(NumNodes);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<Ncontainer.size1(); igauss++)
        {
             noalias(data.N) = row(Ncontainer, igauss);

             ComputeConstitutiveResponse_AIR(data, air_density, air_nu, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);
            ComputeGaussPointLHSContribution(lhs_local, data);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
            noalias(rLeftHandSideMatrix) += weight*lhs_local;
            noalias(rRightHandSideVector) += weight*rhs_local;
        }

//         //assign AIR_DENSITY to density
//         const double air_density = GetProperties()[DENSITY_AIR];
//         for (unsigned int i = 0; i < NumNodes; i++)
//             data.rho[i] = air_density;
//         const double air_nu = GetProperties()[DYNAMIC_VISCOSITY]; //ATTENTION: not using here the real visosity of air
//
//         for (unsigned int i = 0; i < NumNodes; i++) //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//             data.N[i] = 0.25;
//
//         for(unsigned int igauss = 0; igauss<1; igauss++)
//         {
//             ComputeConstitutiveResponse_AIR(data,air_density, air_nu, rCurrentProcessInfo);
//
//             Stokes3D::ComputeGaussPointRHSContribution(rhs_local, data);
//             Stokes3D::ComputeGaussPointLHSContribution(lhs_local, data);
//
//             //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
//             noalias(rLeftHandSideMatrix) = lhs_local; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//             noalias(rRightHandSideVector) = rhs_local; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//         }
//
//         rLeftHandSideMatrix  *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//         rRightHandSideVector *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!

    }


    template<int MatrixSize, int NumNodes>
    void ComputeElementAsFLUID(BoundedMatrix<double,MatrixSize, MatrixSize>& lhs_local,
                               array_1d<double,MatrixSize>& rhs_local,
                               Matrix& rLeftHandSideMatrix,
                               Vector& rRightHandSideVector,
                               const double& Volume,
                               element_data<4,3>& data,
                               BoundedMatrix<double,NumNodes, NumNodes>& Ncontainer,
                               ProcessInfo& rCurrentProcessInfo)
    {
        const double weight = Volume/static_cast<double>(NumNodes);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
        noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
        for(unsigned int igauss = 0; igauss<Ncontainer.size1(); igauss++)
        {
             noalias(data.N) = row(Ncontainer, igauss);

             ComputeConstitutiveResponse(data, rCurrentProcessInfo);

            ComputeGaussPointRHSContribution(rhs_local, data);
            ComputeGaussPointLHSContribution(lhs_local, data);

            //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
            noalias(rLeftHandSideMatrix) += weight*lhs_local;
            noalias(rRightHandSideVector) += weight*rhs_local;
        }


    }



//         for (unsigned int i = 0; i < NumNodes; i++) //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//             data.N[i] = 0.25;
//
//         for(unsigned int igauss = 0; igauss<1; igauss++) //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//         {
//
//             ComputeConstitutiveResponse(data, rCurrentProcessInfo);
//
//             Stokes3D::ComputeGaussPointRHSContribution(rhs_local, data);
//             Stokes3D::ComputeGaussPointLHSContribution(lhs_local, data);
//
//             //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
//             noalias(rLeftHandSideMatrix) = lhs_local; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//             noalias(rRightHandSideVector) = rhs_local; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//         }
//
//         rLeftHandSideMatrix  *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//         rRightHandSideVector *= Volume; //ATTENTION DELIBERATELY USING ONE SINGLE GAUSS POINT!!
//
//     }
//


    //ATTENTION: here multiple integration points are used. For this reason the methods used must be reimplemented in the current element
    template<int MatrixSize, int NumNodes>
    void ComputeElementAsMIXED(BoundedMatrix<double,MatrixSize, MatrixSize>& lhs_local,
                               array_1d<double,MatrixSize>& rhs_local,
                               Matrix& rLeftHandSideMatrix,
                               Vector& rRightHandSideVector,
                               const double& Volume,
                               element_data<4,3>& data,
                               ProcessInfo& rCurrentProcessInfo,
                               const array_1d<double,NumNodes>& distances
                              )
    {
            //here do the splitting to determine the gauss points
            Matrix Ncontainer;
            Vector volumes;
            Vector signs(6); //ATTENTION: this shall be initialized of size 6
            std::vector< Matrix > DNenr;
            Matrix Nenr;
            unsigned int ndivisions = ComputeSplitting(data,Ncontainer, volumes, DNenr, Nenr, signs, distances);


            if(ndivisions == 1)
            {
                //gauss point position
                BoundedMatrix<double,NumNodes, NumNodes> Ncontainer;
                GetShapeFunctionsOnGauss(Ncontainer);

                //cases exist when the element is like not subdivided due to the characteristics of the provided distance
                //in this cases the element is treated as AIR or FLUID depending on the side
                array_1d<double,NumNodes> Ncenter;
                for(unsigned int i=0; i<NumNodes; i++) Ncenter[i]=0.25;
                const double dgauss = inner_prod(distances, Ncenter);
                if(dgauss > 0)
                {
                    ComputeElementAsAIR<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer,rCurrentProcessInfo);
                }
                else
                {
                    ComputeElementAsFLUID<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer, rCurrentProcessInfo);
                }
            }
            else
            {
                BoundedMatrix<double, MatrixSize, NumNodes > Vtot, V;
                BoundedMatrix<double, NumNodes, MatrixSize > Htot, H;
                BoundedMatrix<double, NumNodes, NumNodes> Kee_tot, Kee;
                array_1d<double, NumNodes> rhs_ee_tot, rhs_ee;
                Vtot.clear();
                Htot.clear();
                Kee_tot.clear();
                rhs_ee_tot.clear();

                //loop on gauss points
                noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
                noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
                for(unsigned int igauss = 0; igauss<signs.size(); igauss++)
                {
                    noalias(data.N) = row(Ncontainer, igauss);

                    const double dgauss = inner_prod(distances, data.N); //compute the distance on the gauss point

                    if(dgauss >= 0) //gauss is AIR
                    {
                        //assign AIR_DENSITY to density
                        const double air_density = GetProperties()[DENSITY_AIR];
                        for (unsigned int i = 0; i < NumNodes; i++)
                            data.rho[i] = air_density;
                        const double air_nu = GetProperties()[DYNAMIC_VISCOSITY];

                        ComputeConstitutiveResponse_AIR(data, air_density, air_nu, rCurrentProcessInfo);
                    }
                    else
                    {
                        for (unsigned int i = 0; i < NumNodes; i++)
                            data.rho[i] = GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
                        ComputeConstitutiveResponse(data, rCurrentProcessInfo);
                    }

                    ComputeGaussPointRHSContribution(rhs_local, data); //ATTENTION: uses implementation within the current element since a general integration rule must be employed
                    ComputeGaussPointLHSContribution(lhs_local, data); //ATTENTION: uses implementation within the current element since a general integration rule must be employed
                    const array_1d<double,4> Nenriched = row(Nenr,igauss);
                    ComputeGaussPointEnrichmentContributions(H,V,Kee,rhs_ee, data, distances, Nenriched, DNenr[igauss]);

                    const double weight = volumes[igauss];
                    noalias(rLeftHandSideMatrix) += weight*lhs_local;
                    noalias(rRightHandSideVector) += weight*rhs_local;

                    noalias(Htot) += weight*H;
                    noalias(Vtot) += weight*V;
                    noalias(Kee_tot) += weight*Kee;
                    noalias(rhs_ee_tot) += weight*rhs_ee;
                }

                 CondenseEnrichment(rLeftHandSideMatrix,rRightHandSideVector,Htot,Vtot,Kee_tot, rhs_ee_tot, volumes, signs, distances);
            }

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
        return "Stokes3DTwoFluid #";
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
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    //this is the symbolic function implementing the element
    virtual void ComputeGaussPointLHSContribution(BoundedMatrix<double,16,16>& lhs, const element_data<4,3>& data);
    virtual void ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data<4,3>& data);
    virtual void ComputeGaussPointEnrichmentContributions(
        BoundedMatrix<double,4,16>& H,
        BoundedMatrix<double,16,4>& V,
        BoundedMatrix<double,4,4>&  Kee,
        array_1d<double,4>& rhs_ee,
        const element_data<4,3>& data,
        const array_1d<double,4>& distances,
        const array_1d<double,4>& Nenr,
        const BoundedMatrix<double,4,4>& DNenr
    );

    ///@}
    ///@name Protected Operators
    ///@{

    Stokes3DTwoFluid() : Stokes3D()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{



    unsigned int ComputeSplitting(
        const element_data<4,3>& data,
        Matrix& Ncontainer, Vector& volumes,
        std::vector< Matrix >& DNenr,
        Matrix& Nenr, Vector& signs,
        const array_1d<double,4>& distances)
    {
        Vector el_distances = distances;
        const unsigned int NumNodes = 4;
        const unsigned int Dim = 3;
        Matrix coords(NumNodes, Dim);

        //fill coordinates
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
            for (unsigned int j = 0; j < Dim; j++)
                coords(i, j) = xyz[j];
        }

        unsigned int ndivisions = EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, data.DN_DX, el_distances, volumes, Ncontainer, signs, DNenr, Nenr);
        return ndivisions;
    }

    void CondenseEnrichment(Matrix& rLeftHandSideMatrix,Vector& rRightHandSideVector,
                            const BoundedMatrix<double,4,16>& Htot,
                            const BoundedMatrix<double,16,4>& Vtot,
                            BoundedMatrix<double,4,4>& Kee_tot,
                            array_1d<double,4>& Renr,
                            const Vector& volumes,
                            const Vector& signs,
                            const array_1d<double,4> distances
                           )
    {
        const double Dim = 3;
        const double min_area_ratio = 1e-6;

        double positive_volume = 0.0;
        double negative_volume = 0.0;
        for (unsigned int igauss = 0; igauss < volumes.size(); igauss++)
        {
            double wGauss = volumes[igauss];

            if(signs[igauss] >= 0) //check positive and negative volume
                positive_volume += wGauss;
            else
                negative_volume += wGauss;
        }
        const double Vol = positive_volume + negative_volume;

        double max_diag = 0.0;
        for(unsigned int k=0; k<Dim+1; k++)
            if(std::abs(Kee_tot(k,k) ) > max_diag) max_diag = std::abs(Kee_tot(k,k) );
        if(max_diag == 0) max_diag = 1.0;

        if(positive_volume/Vol < min_area_ratio)
        {
            for(unsigned int i=0; i<Dim+1; i++)
            {
                if(distances[i] >= 0.0)
                {
                    Kee_tot(i,i) += 1000.0*max_diag;
                }
            }
        }
        if(negative_volume/Vol < min_area_ratio)
        {
            for(unsigned int i=0; i<Dim+1; i++)
            {
                if(distances[i] < 0.0)
                {
                    Kee_tot(i,i) += 1000.0*max_diag;
                }
            }
        }

        //"weakly" impose continuity
        for(unsigned int i=0; i<Dim; i++)
        {
            const double di = std::abs(distances[i]);

            for(unsigned int j=i+1; j<Dim+1; j++)
            {
                const double dj =  std::abs(distances[j]);

                if( distances[i]*distances[j] < 0.0) //cut edge
                {
                    double sum_d = di+dj;
                    double Ni = dj/sum_d;
                    double Nj = di/sum_d;

                    double penalty_coeff = max_diag*0.001; // h/BDFVector[0];
                    Kee_tot(i,i) += penalty_coeff * Ni*Ni;
                    Kee_tot(i,j) -= penalty_coeff * Ni*Nj;
                    Kee_tot(j,i) -= penalty_coeff * Nj*Ni;
                    Kee_tot(j,j) += penalty_coeff * Nj*Nj;

                }
            }
        }

        //add to LHS enrichment contributions
        double det;
        BoundedMatrix<double,4,4> inverse_diag;
        MathUtils<double>::InvertMatrix(Kee_tot, inverse_diag,det);

        const BoundedMatrix<double,4,16> tmp = prod(inverse_diag,Htot);
        noalias(rLeftHandSideMatrix) -= prod(Vtot,tmp);

        const array_1d<double,4> tmp2 = prod(inverse_diag,Renr);
        noalias(rRightHandSideVector) -= prod(Vtot,tmp2);

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
    virtual void ComputeConstitutiveResponse_AIR(element_data<4,3>& data, const double rho, const double nu,  ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int strain_size = 6;

        if(data.C.size1() != strain_size)
            data.C.resize(strain_size,strain_size,false);
        if(data.stress.size() != strain_size)
            data.stress.resize(strain_size,false);

        //compute strain
        Vector strain = CalculateStrain();

        //here we shall call the constitutive law
        data.C.clear();
        data.C(0,0) = 2.0*nu;
        data.C(1,1) = 2.0*nu;
        data.C(2,2) = 2.0*nu;
        data.C(3,3) = nu;
        data.C(4,4) = nu;
        data.C(5,5) = nu;

        const double c2 = nu;
        const double c1 = 2.0*c2;
        data.stress[0] =  c1*strain[0];
        data.stress[1] =  c1*strain[1];
        data.stress[2] =  c1*strain[2];
        data.stress[3] =  c2*strain[3];
        data.stress[4] =  c2*strain[4];
        data.stress[5] =  c2*strain[5];
    }

    //compute strain

    Vector CalculateStrain()
    {
        const unsigned int NumNodes = 4;
        const unsigned int Dim = 3;
        const unsigned int strain_size = 6;

        //struct to pass around the data
        element_data<NumNodes,Dim> data;

        //getting data for the given geometry
        double Volume;
        BoundedMatrix<double,NumNodes,Dim> v, DN;
        array_1d<double,NumNodes> N;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN, N, Volume);

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            for(unsigned int k=0; k<Dim; k++)
                v(i,k)   = vel[k];
        }

        Vector strain(strain_size);
        strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
        strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
        strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
        strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
        strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
        strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
        return strain;
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

#endif // KRATOS_STOKES_ELEMENT_TWOFLUID_3D_INCLUDED  defined


