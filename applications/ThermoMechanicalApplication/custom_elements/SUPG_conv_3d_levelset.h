//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Rubio
//


#if !defined(KRATOS_SUPG_CONV_LEVELSET )
#define  KRATOS_SUPG_CONV_LEVELSET


// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "custom_elements/SUPG_conv_3d.h"



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


class SUPGConvLevelSet: public SUPGConv3D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of SUPGConvDiffPhaseChange3DLinearized
    KRATOS_CLASS_POINTER_DEFINITION(SUPGConvLevelSet);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    SUPGConvLevelSet(IndexType NewId, GeometryType::Pointer pGeometry);
    SUPGConvLevelSet(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~SUPGConvLevelSet();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

	void CalculatePenalty(VectorType& penalty);

	void CalculateDistanceGradient(array_1d<double,3>& grad_D);

	//void CalculateRightHandSide(VectorType	ightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    //void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    //void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);



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
        return "SUPGConvLevelSet " ;
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

    ///@}
    ///@name Protected Operators
    ///@{
    //double ComputeTau(const double h, const double k, const array_1d<double,3>& a);

    //void ComputeEffectiveSpeficifHeat(array_1d<double,4>& c,
    //                                  const array_1d<double,4>& temperatures,
    //                                  const double fluid_T,
    //                                  const double solid_T,
    //                                  const double latent_heat
    //                                 );

    //const double ComputeDiscontinuityCapturingDiffusion(
    //    const BoundedMatrix<double, 4, 3 >& DN_DX,
    //    const array_1d<double,3>& gradT,
    //    const double& norm_gradT,
    //    const double& residual,
    //    const double& reference_temperature
    //);

    SUPGConvLevelSet() : SUPGConv3D() {}

    ///@}
    ///@name Protected Operations
    ///@{


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

    // private:

    /// To Calculate tau
    /**
     * @return tau: multiplied by Residual of the momentum equation
    */
//     virtual void CalculateTau(array_1d<double, 2 >& ms_adv_vel, double& tau,const double& K, const double time, const double area, const ProcessInfo& rCurrentProcessInfo);

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SUPGConv3D);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SUPGConv3D);
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


    /// Assignment operator.
    //Fluid2DASGS& operator=(const SUPGConvDiffPhaseChange3D& rOther);

    /// Copy constructor.
    //Fluid2DASGS(const SUPGConvDiffPhaseChange3D& rOther);


    ///@}

}; // Class SUPGConvDiffPhaseChange3D

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
				    SUPGConvDiffPhaseChange3D& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
				    const SUPGConvDiffPhaseChange3D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_SUPG_CONV_DIFF_PHASE_CHANGE_3D_LINEARIZED_INCLUDED  defined
