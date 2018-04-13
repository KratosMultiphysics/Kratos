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
	void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) override;


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



private:

	unsigned int ComputeSplitting(
		TElementData& data,
		MatrixType& shape_functions,
		ShapeFunctionDerivativesArrayType& shape_derivatives,
		std::vector<MatrixType>& DNenr,
		MatrixType& Nner);

	void CondenseEnrichment(
		TElementData& data,
		Matrix& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		MatrixType& V,
		MatrixType& H,
		MatrixType& K_ee,
		VectorType& rhs_ee);

	void CalculateMaterialPropertiesAtGaussPoint(TElementData& data);

	template<class T>
	bool InvertMatrix(const T& input, T& inverse);


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


