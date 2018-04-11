//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Daniel Diez
//  Co-authors:      Ruben Zorrilla
//



#if !defined(KRATOS_TWO_FLUID_NAVIER_STOKES)
#define  KRATOS_TWO_FLUID_NAVIER_STOKES

// System includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "utilities/geometry_utilities.h"
#include "includes/cfd_variables.h"
#include "utilities/split_tetrahedra.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "custom_elements/fluid_element.h"



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


template< class TElementData >
class TwoFluidNavierStokes : public FluidElement<TElementData>
{
public:

    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidNavierStokes);

    ///@name Type Definitions
    ///@{
    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

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
    typedef typename FluidElement<TElementData>::ShapeFunctionsType ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesType ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef typename FluidElement<TElementData>::ShapeFunctionDerivativesArrayType ShapeFunctionDerivativesArrayType;

    constexpr static unsigned int Dim = FluidElement<TElementData>::Dim;
    constexpr static unsigned int NumNodes = FluidElement<TElementData>::NumNodes;
    constexpr static unsigned int BlockSize = FluidElement<TElementData>::BlockSize;
    constexpr static unsigned int LocalSize = FluidElement<TElementData>::LocalSize;

    constexpr static unsigned int StrainSize = (Dim - 1) * 3;


    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constuctor.
    /**
    * @param NewId Index number of the new element (optional)
    */
    TwoFluidNavierStokes(IndexType NewId = 0);

    /// Constructor using an array of nodes.
    /**
    * @param NewId Index of the new element
    * @param ThisNodes An array containing the nodes of the new element
    */
    TwoFluidNavierStokes(IndexType NewId, const NodesArrayType& ThisNodes);

    /// Constructor using a geometry object.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    */
    TwoFluidNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry);

    /// Constuctor using geometry and properties.
    /**
    * @param NewId Index of the new element
    * @param pGeometry Pointer to a geometry object
    * @param pProperties Pointer to the element's properties
    */
    TwoFluidNavierStokes(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties);

    /// Destructor.
    virtual ~TwoFluidNavierStokes();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// Create a new element of this type
    /**
    * Returns a pointer to a new TwoFluidNavierStokes element, created using given input.
    * @param NewId the ID of the new element
    * @param ThisNodes the nodes of the new element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(IndexType NewId,
        NodesArrayType const& ThisNodes,
        Properties::Pointer pProperties) const override;

    /// Create a new element of this type using given geometry
    /**
    * Returns a pointer to a new FluidElement element, created using given input.
    * @param NewId the ID of the new element
    * @param pGeom a pointer to the geomerty to be used to create the element
    * @param pProperties the properties assigned to the new element
    * @return a Pointer to the new element
    */
    Element::Pointer Create(IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo &rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;


    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    ///@}


  //  void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  //  {
  //      KRATOS_TRY

  //      const unsigned int NumNodes = 4;
  //      const unsigned int Dim = 3;
  //      const int ndofs = Dim + 1;
  //      const unsigned int MatrixSize = NumNodes*ndofs;



  //      if (rLeftHandSideMatrix.size1() != MatrixSize)
  //          rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

  //      if (rRightHandSideVector.size() != MatrixSize)
  //          rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!


  //      //getting data for the given geometry

  //      double Volume;
  //      GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);

  //      //compute element size
		//data.h = ComputeH<4,3>(data.DN_DX, Volume);
  //      
  //      //gauss point position
  //      bounded_matrix<double,NumNodes, NumNodes> Ncontainer;
  //      GetShapeFunctionsOnGauss(Ncontainer);

  //      //database access to all of the variables needed
  //      const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];
  //      data.bdf0 = BDFVector[0];
  //      data.bdf1 = BDFVector[1];
  //      data.bdf2 = BDFVector[2];
  //      data.dyn_tau_coeff = rCurrentProcessInfo[DYNAMIC_TAU];
  //      data.delta_t = rCurrentProcessInfo[DELTA_TIME]; 

  //      array_1d<double, NumNodes> distances;
  //      for (unsigned int i = 0; i < NumNodes; i++)
  //      {
  //          const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
  //          const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
  //          const array_1d<double,3>& vel_n = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
  //          const array_1d<double,3>& vel_nn = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,2);

  //          for(unsigned int k=0; k<Dim; k++)
  //          {
  //              data.v(i,k)   = vel[k];
  //              data.vn(i,k)  = vel_n[k];
  //              data.vnn(i,k) = vel_nn[k];
  //              data.f(i,k)   = body_force[k];
  //          }

  //          data.p[i] = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
  //          data.rho[i] = GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
  //          data.nu[i] = GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
  //          distances[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
  //      }

  //      //allocate memory needed
  //      bounded_matrix<double,MatrixSize, MatrixSize> lhs_local;
  //      array_1d<double,MatrixSize> rhs_local;

  //      unsigned int npos=0, nneg=0;
  //      for (unsigned int i = 0; i < NumNodes; i++)
  //      {
  //          if(distances[i] > 0)
  //              npos++;
  //          else
  //              nneg++;
  //      }


  //      //here we decide if the element is all FLUID/AIR/MIXED
  //      if(npos == NumNodes) //all AIR 
  //      {
  //          ComputeElementAsAIR<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer, rCurrentProcessInfo,distances);
  //      }
  //      else if (nneg == NumNodes) //all FLUID 
  //      {
  //          ComputeElementAsFLUID<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, Ncontainer, rCurrentProcessInfo, distances);
  //      }
  //      else //element includes both FLUID and AIR 
  //      {
  //          ComputeElementAsMIXED<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data, rCurrentProcessInfo, distances);
  //      }
  //          

  //      KRATOS_CATCH("Error in NavierStokesEnr3D Element Symbolic")
  //  }













    /// Checks the input and that all required Kratos variables have been registered.
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
     * @return 0 if no errors were found.
     */
    //virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    //{
    //    KRATOS_TRY

    //    // Perform basic element checks
    //    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    //    if(ErrorCode != 0) return ErrorCode;

    //    // Check that all required variables have been registered
    //    if(VELOCITY.Key() == 0)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
    //    if(DISTANCE.Key() == 0)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"DISTANCE Key is 0. Check if the application was correctly registered.","");

    //    if(PRESSURE.Key() == 0)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");
    //    if(DENSITY.Key() == 0)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check if the application was correctly registered.","");
    //    if(DYNAMIC_TAU.Key() == 0)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check if the application was correctly registered.","");
    //    if(DELTA_TIME.Key() == 0)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check if the application was correctly registered.","");

    //    // Checks on nodes
    //    
    //    //check Properties
    //    if(GetProperties().Has(DENSITY_AIR) == false)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY_AIR is not set","");
    //    if(GetProperties().Has(CONSTITUTIVE_LAW) == false)
    //        KRATOS_THROW_ERROR(std::invalid_argument,"CONSTITUTIVE_LAW is not set","");
    //    
    //    //check constitutive CONSTITUTIVE_LAW
    //    GetProperties().GetValue(CONSTITUTIVE_LAW)->Check(GetProperties(),GetGeometry(),rCurrentProcessInfo);

    //    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    //    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
    //    {
    //        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
    //            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
    //        if(this->GetGeometry()[i].SolutionStepsDataHas(DISTANCE) == false)
    //            KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data for node ",this->GetGeometry()[i].Id());
    //        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
    //            KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
    //        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
    //            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
    //        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
    //                this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
    //                this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
    //            KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
    //        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
    //            KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
    //    }

    //    return 0;

    //    KRATOS_CATCH("");
    //}

    //virtual void Calculate(const Variable<double>& rVariable,
    //                       double& Output,
    //                       const ProcessInfo& rCurrentProcessInfo)
    //{
    //    KRATOS_TRY

    //    if(rVariable == HEAT_FLUX) //compute the heat flux per unit volume induced by the shearing
    //    {
    //        const unsigned int NumNodes = 4;
    //        const unsigned int Dim = 3;
    //        const unsigned int strain_size = 6;

    //        double distance_center = 0.0;
    //        for(unsigned int i=0; i<GetGeometry().size(); i++)
    //            distance_center += GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    //        distance_center/=static_cast<double>(GetGeometry().size());
    //        
    //        if(distance_center > 0) //AIR 
    //        {
    //            Output=0.0;
    //        }
    //        else //OTHER MATERIAL
    //        {
    //            //struct to pass around the data
    //            element_data<NumNodes,Dim> data;

    //            //getting data for the given geometry
    //            double Volume;
    //            GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, Volume);


    //            for (unsigned int i = 0; i < NumNodes; i++)
    //            {
    //                const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);

    //                for(unsigned int k=0; k<Dim; k++)
    //                {
    //                    data.v(i,k)   = vel[k];
    //                }
    //            }

    //            if (data.stress.size() != strain_size) data.stress.resize(strain_size,false);

    //            //const bounded_matrix<double,NumNodes,Dim>& v = data.v;
    //            //const bounded_matrix<double,NumNodes,Dim>& DN = data.DN_DX;

    //            //compute strain
    //            Vector strain(strain_size);
    //            ComputeStrain(data, strain_size, strain);

    //            //create constitutive law parameters:
    //            ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //            //set constitutive law flags:
    //            Flags& ConstitutiveLawOptions=Values.GetOptions();
    //            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    //            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    //            //this is to pass the shape functions. Unfortunately it is needed to make a copy to a flexible size vector
    //            const Vector Nvec(data.N);
    //            Values.SetShapeFunctionsValues(Nvec);

    //            Values.SetStrainVector(strain); //this is the input parameter
    //            Values.SetStressVector(data.stress); //this is an ouput parameter
    ////             Values.SetConstitutiveMatrix(data.C);      //this is an ouput parameter

    //            //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //            //this is ok under the hypothesis that no history dependent behaviour is employed
    //            mp_constitutive_law->CalculateMaterialResponseCauchy(Values);

    //            Output = inner_prod(data.stress, strain);
    //        }
    //    }

    //    KRATOS_CATCH("")
    //}

    //template<int MatrixSize, int NumNodes>
    //void ComputeElementAsAIR(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
    //                         array_1d<double,MatrixSize>& rhs_local,
    //                         Matrix& rLeftHandSideMatrix,
    //                         Vector& rRightHandSideVector,
    //                         const double& Volume,
    //                         element_data<4,3>& data,
    //                         bounded_matrix<double,NumNodes, NumNodes>& Ncontainer,
    //                         ProcessInfo& rCurrentProcessInfo,
    //                         const array_1d<double,NumNodes>& distances)
    //{
    //    const double air_density = GetProperties()[DENSITY_AIR];
    //    data.tau1_coeff = 1;
    //    for (unsigned int i = 0; i < NumNodes; i++)
    //        data.rho[i] = air_density;
    //    const double air_nu = GetProperties()[DYNAMIC_VISCOSITY]; //ATTENTION: not using here the real visosity of air

    //    const double weight = Volume/static_cast<double>(NumNodes);
    //    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
    //    noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
    //    for(unsigned int igauss = 0; igauss<Ncontainer.size1(); igauss++) 
    //    {
    //         noalias(data.N) = row(Ncontainer, igauss); 
 
    //         ComputeConstitutiveResponse_AIR(data, air_density, air_nu, rCurrentProcessInfo);
	
    //        NavierStokesEnr3D::ComputeGaussPointRHSContribution(rhs_local, data);
    //        NavierStokesEnr3D::ComputeGaussPointLHSContribution(lhs_local, data);

    //        //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    //        noalias(rLeftHandSideMatrix) += weight*lhs_local; 
    //        noalias(rRightHandSideVector) += weight*rhs_local; 
    //    }

    //}


    //template<int MatrixSize, int NumNodes>
    //void ComputeElementAsFLUID(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
    //                           array_1d<double,MatrixSize>& rhs_local,
    //                           Matrix& rLeftHandSideMatrix,
    //                           Vector& rRightHandSideVector,
    //                           const double& Volume,
    //                           element_data<4,3>& data,
    //                           bounded_matrix<double,NumNodes, NumNodes>& Ncontainer,
    //                           ProcessInfo& rCurrentProcessInfo,
    //                           const array_1d<double,NumNodes>& distances)
    //{
    //    const double weight = Volume/static_cast<double>(NumNodes);
    //    data.tau1_coeff = 1;
    //    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
    //    noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
	
    //    for(unsigned int igauss = 0; igauss<Ncontainer.size1(); igauss++) 
    //    {
    //         noalias(data.N) = row(Ncontainer, igauss); 
 
    //         ComputeConstitutiveResponse(data, rCurrentProcessInfo);

    //        NavierStokesEnr3D::ComputeGaussPointRHSContribution(rhs_local, data);
    //        NavierStokesEnr3D::ComputeGaussPointLHSContribution(lhs_local, data);

    //        //here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/NumNodes
    //        noalias(rLeftHandSideMatrix) += weight*lhs_local; 
    //        noalias(rRightHandSideVector) += weight*rhs_local; 
    //    }


    //}
    //
    //    
    //
    ////ATTENTION: here multiple integration points are used. For this reason the methods used must be reimplemented in the current element
    //template<int MatrixSize, int NumNodes>
    //void ComputeElementAsMIXED(bounded_matrix<double,MatrixSize, MatrixSize>& lhs_local,
    //                           array_1d<double,MatrixSize>& rhs_local,
    //                           Matrix& rLeftHandSideMatrix,
    //                           Vector& rRightHandSideVector,
    //                           const double& Volume,
    //                           element_data<4,3>& data,
    //                           ProcessInfo& rCurrentProcessInfo,
    //                           const array_1d<double,NumNodes>& distances
    //                          )
    //{
    //        //here do the splitting to determine the gauss points
    //        Matrix Ncontainer;
    //        Vector volumes;
    //        Vector signs(6); //ATTENTION: this shall be initialized of size 6
    //        std::vector< Matrix > DNenr;
    //        Matrix Nenr;
    //        unsigned int ndivisions = ComputeSplitting(data,Ncontainer, volumes, DNenr, Nenr, signs, distances);
    //        data.tau1_coeff = 1;






    //        if(ndivisions == 1)
    //        {
    //            //gauss point position
    //            bounded_matrix<double,NumNodes, NumNodes> Ncontainer;
    //            GetShapeFunctionsOnGauss(Ncontainer);
    //    
    //            //cases exist when the element is like not subdivided due to the characteristics of the provided distance
    //            //in this cases the element is treated as AIR or FLUID depending on the side
    //            array_1d<double,NumNodes> Ncenter;

    //            for(unsigned int i=0; i<NumNodes; i++) Ncenter[i]=0.25;
    //            const double dgauss = inner_prod(distances, Ncenter);
    //            if(dgauss > 0)
    //            {
    //                ComputeElementAsAIR<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data,Ncontainer,  rCurrentProcessInfo, distances);
    //            }
    //            else
    //            {
    //                ComputeElementAsFLUID<MatrixSize,NumNodes>(lhs_local, rhs_local, rLeftHandSideMatrix, rRightHandSideVector, Volume, data,Ncontainer,  rCurrentProcessInfo, distances);
    //            }
    //        }
    //        else
    //        {
    //            boost::numeric::ublas::bounded_matrix<double, MatrixSize, NumNodes > Vtot, V;
    //            boost::numeric::ublas::bounded_matrix<double, NumNodes, MatrixSize > Htot, H;
    //            boost::numeric::ublas::bounded_matrix<double, NumNodes, NumNodes> Kee_tot, Kee;
    //            array_1d<double, NumNodes> rhs_ee_tot, rhs_ee;
    //            Vtot.clear();
    //            Htot.clear();
    //            Kee_tot.clear();
    //            rhs_ee_tot.clear();

    //            //loop on gauss points
    //            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
    //            noalias(rRightHandSideVector) = ZeroVector(MatrixSize);
    //            for(unsigned int igauss = 0; igauss<signs.size(); igauss++)
    //            {
    //                noalias(data.N) = row(Ncontainer, igauss);

    //                const double dgauss = inner_prod(distances, data.N); //compute the distance on the gauss point

    //                if(dgauss >= 0) //gauss is AIR
    //                {
    //                    //assign AIR_DENSITY to density
    //                    const double air_density = GetProperties()[DENSITY_AIR];
    //                    for (unsigned int i = 0; i < NumNodes; i++)
    //                        data.rho[i] = air_density;
    //                    const double air_nu = GetProperties()[DYNAMIC_VISCOSITY];

    //                    ComputeConstitutiveResponse_AIR(data, air_density, air_nu, rCurrentProcessInfo);
    //                }
    //                else
    //                {
    //                    for (unsigned int i = 0; i < NumNodes; i++)
    //                        data.rho[i] = GetGeometry()[i].FastGetSolutionStepValue(DENSITY);
    //                    ComputeConstitutiveResponse(data, rCurrentProcessInfo);
    //                }

    //                const array_1d<double,4> Nenriched = row(Nenr,igauss);
    //                NavierStokesEnr3D::ComputeGaussPointRHSContribution(rhs_local, data); //ATTENTION: uses implementation within the current element since a general integration rule must be employed
    //                NavierStokesEnr3D::ComputeGaussPointLHSContribution(lhs_local, data); //ATTENTION: uses implementation within the current element since a general integration rule must be employed
    //                NavierStokesEnr3D::ComputeGaussPointEnrichmentContributions(H,V,Kee,rhs_ee, data, distances, Nenriched, DNenr[igauss]);

    //                const double weight = volumes[igauss];
    //                noalias(rLeftHandSideMatrix) += weight*lhs_local;
    //                noalias(rRightHandSideVector) += weight*rhs_local;

    //                noalias(Htot) += weight*H;
    //                noalias(Vtot) += weight*V;
    //                noalias(Kee_tot) += weight*Kee;
    //                noalias(rhs_ee_tot) += weight*rhs_ee;
    //            }

    //             CondenseEnrichment(rLeftHandSideMatrix,rRightHandSideVector,Htot,Vtot,Kee_tot, rhs_ee_tot, volumes, signs, distances);
    //        }

    //    }
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{

protected:

    ///@}
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
    void AddTimeIntegratedSystem(
        TElementData& rData,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    void AddTimeIntegratedLHS(
        TElementData& rData,
        MatrixType& rLHS) override;

    void AddTimeIntegratedRHS(
        TElementData& rData,
        VectorType& rRHS) override;

    void AddBoundaryIntegral(
        TElementData& rData,
        const Vector& rUnitNormal,
        MatrixType& rLHS,
        VectorType& rRHS) override;

    void ComputeGaussPointLHSContribution(
        TElementData& rData,
        MatrixType& rLHS);

    void ComputeGaussPointRHSContribution(
        TElementData& rData,
        VectorType& rRHS);

	void ComputeGaussPointEnrichmentContributions(
		TElementData& rData,
		MatrixType& rV,
		MatrixType& rH,
		MatrixType& rKee,
		VectorType& rRHS_ee);


    //void ComputeGaussPointLHSContribution(bounded_matrix<double,16,16>& lhs, const element_data<4,3>& data);
    //void ComputeGaussPointRHSContribution(array_1d<double,16>& rhs,const element_data<4,3>& data);
    //void ComputeGaussPointEnrichmentContributions(
    //    boost::numeric::ublas::bounded_matrix<double,4,16>& H,
    //    boost::numeric::ublas::bounded_matrix<double,16,4>& V,
    //    boost::numeric::ublas::bounded_matrix<double,4,4>&  Kee,
    //    array_1d<double,4>& rhs_ee,
    //    const element_data<4,3>& data,
    //    const array_1d<double,4>& distances,
    //    const array_1d<double,4>& Nenr,
    //    const boost::numeric::ublas::bounded_matrix<double,4,4>& DNenr
    //);


private:

    //// 3D tetrahedra shape functions values at Gauss points 

    //template< unsigned int TNumNodes, unsigned int TDim>    
    //double ComputeH(boost::numeric::ublas::bounded_matrix<double,TNumNodes, TDim>& DN_DX, const double Volume)
    //{
    //    double h=0.0;
    //             for(unsigned int i=0; i<TNumNodes; i++)
    //    {
    //        double h_inv = 0.0;
    //        for(unsigned int k=0; k<TDim; k++)
    //        {
    //            h_inv += DN_DX(i,k)*DN_DX(i,k);
    //        }
    //        h += 1.0/h_inv;
    //    }
    //    h = sqrt(h)/static_cast<double>(TNumNodes);
    //    return h;
    //}

    //virtual void ComputeConstitutiveResponse(element_data<4,3>& data, ProcessInfo& rCurrentProcessInfo)
    //{
    //    //const unsigned int nnodes = 4;
    //    //const unsigned int dim = 3;
    //    const unsigned int strain_size = 6;
    //    
    //    if(data.C.size1() != strain_size)
    //        data.C.resize(strain_size,strain_size,false);
    //    if(data.stress.size() != strain_size)
    //        data.stress.resize(strain_size,false);
    //   //const double rho = inner_prod(data.rho, data.N);
    //    
    //    //compute strain
    //    Vector strain(strain_size);
    //    ComputeStrain(data, strain_size, strain);      
    //    
    //    //create constitutive law parameters:
    //    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    //    
    //    const Vector Nvec(data.N);
    //    Values.SetShapeFunctionsValues(Nvec);
    //    //set constitutive law flags:
    //    Flags& ConstitutiveLawOptions=Values.GetOptions();
    //    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    //    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    //    
    //    Values.SetStrainVector(strain); //this is the input parameter
    //    Values.SetStressVector(data.stress); //this is an ouput parameter
    //    Values.SetConstitutiveMatrix(data.C);      //this is an ouput parameter

    //    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element. 
    //    //this is ok under the hypothesis that no history dependent behaviour is employed
    //    mp_constitutive_law->CalculateMaterialResponseCauchy(Values);
    //    
    //}

    //void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,4,4>& Ncontainer)
    //{
    //    Ncontainer(0,0) = 0.58541020; Ncontainer(0,1) = 0.13819660; Ncontainer(0,2) = 0.13819660; Ncontainer(0,3) = 0.13819660;
    //    Ncontainer(1,0) = 0.13819660; Ncontainer(1,1) = 0.58541020; Ncontainer(1,2) = 0.13819660; Ncontainer(1,3) = 0.13819660;	
    //    Ncontainer(2,0) = 0.13819660; Ncontainer(2,1) = 0.13819660; Ncontainer(2,2) = 0.58541020; Ncontainer(2,3) = 0.13819660;
    //    Ncontainer(3,0) = 0.13819660; Ncontainer(3,1) = 0.13819660; Ncontainer(3,2) = 0.13819660; Ncontainer(3,3) = 0.58541020;
    //}

    ////2D triangle shape functions values at Gauss points 
    //void GetShapeFunctionsOnGauss(boost::numeric::ublas::bounded_matrix<double,3,3>& Ncontainer)
    //{
    //    const double one_sixt = 1.0/6.0;
    //    const double two_third = 2.0/3.0;
    //    Ncontainer(0,0) = one_sixt; Ncontainer(0,1) = one_sixt; Ncontainer(0,2) = two_third; 
    //    Ncontainer(1,0) = one_sixt; Ncontainer(1,1) = two_third; Ncontainer(1,2) = one_sixt; 
    //    Ncontainer(2,0) = two_third; Ncontainer(2,1) = one_sixt; Ncontainer(2,2) = one_sixt; 
    //}

    //void ComputeStrain(element_data<4,3>& data, const unsigned int& strain_size, Vector& strain)
    //{
    //    const bounded_matrix<double, 4, 3>& v = data.v;
    //    const bounded_matrix<double, 4, 3>& DN = data.DN_DX;
    //    
    //    // Compute strain (B*v)
    //    // 3D strain computation
    //    if (strain_size == 6)
    //    {
    //        strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
    //        strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
    //        strain[2] = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
    //        strain[3] = DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
    //        strain[4] = DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
    //        strain[5] = DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
    //    }
    //    // 2D strain computation
    //    else if (strain_size == 3)
    //    {                
    //        strain[0] = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
    //        strain[1] = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    //        strain[2] = DN(0,1)*v(0,0) + DN(0,0)*v(0,1) + DN(1,1)*v(1,0) + DN(1,0)*v(1,1) + DN(2,1)*v(2,0) + DN(2,0)*v(2,1);
    //    }
    //}
    //
    //void Initialize()
    //{
    //    KRATOS_TRY

    //    mp_constitutive_law = GetProperties()[CONSTITUTIVE_LAW]->Clone();
    //    mp_constitutive_law->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues(), 0 ) );
    //    
    //    KRATOS_CATCH( "" )
    //}

    //ConstitutiveLaw::Pointer mp_constitutive_law = nullptr;


    //unsigned int ComputeSplitting(
    //    const element_data<4,3>& data,
    //    Matrix& Ncontainer, Vector& volumes,
    //    std::vector< Matrix >& DNenr,
    //    Matrix& Nenr,
	   // Vector& signs,
    //    const array_1d<double,4>& distances)
    //{
    //    Vector el_distances = distances;
    //    
    //    const unsigned int NumNodes = 4;
    //    const unsigned int Dim = 3;
    //    Matrix coords(NumNodes, Dim);

    //    //fill coordinates
    //    for (unsigned int i = 0; i < NumNodes; i++)
    //    {
    //        const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
    //        for (unsigned int j = 0; j < Dim; j++)
    //            coords(i, j) = xyz[j];
    //    }

    //    unsigned int ndivisions = EnrichmentUtilitiesDuplicateDofs::CalculateTetrahedraEnrichedShapeFuncions(coords, data.DN_DX, el_distances, volumes, Ncontainer, signs, DNenr, Nenr);
    //    return ndivisions;
    //}

    //void CondenseEnrichment(Matrix& rLeftHandSideMatrix,Vector& rRightHandSideVector,
    //                        const boost::numeric::ublas::bounded_matrix<double,4,16>& Htot,
    //                        const boost::numeric::ublas::bounded_matrix<double,16,4>& Vtot,
    //                        boost::numeric::ublas::bounded_matrix<double,4,4>& Kee_tot,
    //                        array_1d<double,4>& Renr,
    //                        const Vector& volumes,
    //                        const Vector& signs,
    //                        const array_1d<double,4> distances
    //                       )
    //{
    //    const double Dim = 3;
    //    const double min_area_ratio = -1e-6;

    //    double positive_volume = 0.0;
    //    double negative_volume = 0.0;
    //    for (unsigned int igauss = 0; igauss < volumes.size(); igauss++)
    //    {
    //        double wGauss = volumes[igauss];

    //        if(signs[igauss] >= 0) //check positive and negative volume
    //            positive_volume += wGauss;
    //        else
    //            negative_volume += wGauss;
    //    }
    //    const double Vol = positive_volume + negative_volume;
		


    //    double max_diag = 0.0;
    //    for(unsigned int k=0; k<Dim+1; k++)
    //        if(fabs(Kee_tot(k,k) ) > max_diag) max_diag = fabs(Kee_tot(k,k) );
    //    if(max_diag == 0) max_diag = 1.0;
    //    


    //    if(positive_volume/Vol < min_area_ratio)
    //    {
    //        for(unsigned int i=0; i<Dim+1; i++)
    //        {
    //            if(distances[i] >= 0.0)
    //            {
    //                Kee_tot(i,i) += 1000.0*max_diag;
    //            }
    //        }
    //    }
    //    if(negative_volume/Vol < min_area_ratio)
    //    {
    //        for(unsigned int i=0; i<Dim+1; i++)
    //        {
    //            if(distances[i] < 0.0)
    //            {
    //                Kee_tot(i,i) += 1000.0*max_diag;
    //            }
    //        }
    //    }

    //    //"weakly" impose continuity
    //    for(unsigned int i=0; i<Dim; i++)
    //    {
    //        const double di = fabs(distances[i]);

    //        for(unsigned int j=i+1; j<Dim+1; j++)
    //        {
    //            const double dj =  fabs(distances[j]);

    //            if( distances[i]*distances[j] < 0.0) //cut edge
    //            {
    //                double sum_d = di+dj;
    //                double Ni = dj/sum_d;
    //                double Nj = di/sum_d;

    //                double penalty_coeff = max_diag*0.001; // h/BDFVector[0];
    //                Kee_tot(i,i) += penalty_coeff * Ni*Ni;
    //                Kee_tot(i,j) -= penalty_coeff * Ni*Nj;
    //                Kee_tot(j,i) -= penalty_coeff * Nj*Ni;
    //                Kee_tot(j,j) += penalty_coeff * Nj*Nj;

    //            }
    //        }
    //    }

    //    //add to LHS enrichment contributions
    //    boost::numeric::ublas::bounded_matrix<double,4,4> inverse_diag;
    //    bool inversion_successful = InvertMatrix<>(Kee_tot,inverse_diag);

    //    if(!inversion_successful )
    //    {
    //        KRATOS_WATCH(distances)
    //        KRATOS_WATCH(positive_volume/Vol)
    //        KRATOS_WATCH(negative_volume/Vol)
    //        KRATOS_WATCH(Kee_tot)
    //        KRATOS_THROW_ERROR(std::logic_error,"error in the inversion of the enrichment matrix for element ",this->Id());
    //    }

    //    const boost::numeric::ublas::bounded_matrix<double,4,16> tmp = prod(inverse_diag,Htot);
    //    noalias(rLeftHandSideMatrix) -= prod(Vtot,tmp);

    //    const array_1d<double,4> tmp2 = prod(inverse_diag,Renr);
    //    noalias(rRightHandSideVector) -= prod(Vtot,tmp2);

    //}


    //template<class T>
    //bool InvertMatrix(const T& input, T& inverse)
    //{
    //    typedef permutation_matrix<std::size_t> pmatrix;

    //    // create a working copy of the input
    //    T A(input);

    //    // create a permutation matrix for the LU-factorization
    //    pmatrix pm(A.size1());

    //    // perform LU-factorization
    //    int res = lu_factorize(A, pm);
    //    if (res != 0)
    //        return false;

    //    // create identity matrix of "inverse"
    //    inverse.assign(identity_matrix<double> (A.size1()));

    //    // backsubstitute to get the inverse
    //    lu_substitute(A, pm, inverse);

    //    return true;
    //}

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
///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Serialization
///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

///@}

///@name Private Operations
///@{
    //virtual void ComputeConstitutiveResponse_AIR(element_data<4,3>& data, const double rho, const double nu,  ProcessInfo& rCurrentProcessInfo)
    //{
    //    //const unsigned int nnodes = 4;
    //    //const unsigned int dim = 3;
    //    const unsigned int strain_size = 6;

    //    if(data.C.size1() != strain_size)
    //        data.C.resize(strain_size,strain_size,false);
    //    if(data.stress.size() != strain_size)
    //        data.stress.resize(strain_size,false);

    //    //compute strain
    //    Vector strain(strain_size);
    //    ComputeStrain(data, strain_size, strain);

    //    //here we shall call the constitutive law
    //    data.C.clear();
    //    data.C(0,0) = 2.0*nu;
    //    data.C(1,1) = 2.0*nu;
    //    data.C(2,2) = 2.0*nu;
    //    data.C(3,3) = nu;
    //    data.C(4,4) = nu;
    //    data.C(5,5) = nu;

    //    const double c2 = nu;
    //    const double c1 = 2.0*c2;
    //    data.stress[0] =  c1*strain[0];
    //    data.stress[1] =  c1*strain[1];
    //    data.stress[2] =  c1*strain[2];
    //    data.stress[3] =  c2*strain[3];
    //    data.stress[4] =  c2*strain[4];
    //    data.stress[5] =  c2*strain[5];
    //}


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
    TwoFluidNavierStokes& operator=(TwoFluidNavierStokes const& rOther);

    /// Copy constructor.
    TwoFluidNavierStokes(TwoFluidNavierStokes const& rOther);




///@}

}; // Class TwoFluidNavierStokes

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TElementData >
inline std::istream& operator >> (std::istream& rIStream,
    TwoFluidNavierStokes<TElementData>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TElementData >
inline std::ostream& operator <<(std::ostream& rOStream,
    const TwoFluidNavierStokes<TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_TWO_FLUID_NAVIER_STOKES


