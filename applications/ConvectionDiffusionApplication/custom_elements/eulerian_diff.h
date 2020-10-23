// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_EULERIAN_DIFFUSION_ELEMENT_INCLUDED )
#define  KRATOS_EULERIAN_DIFFUSION_ELEMENT_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/serializer.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"



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

///formulation described in https://docs.google.com/document/d/13a_zGLj6xORDuLgoOG5LwHI6BwShvfO166opZ815zLY/edit?usp=sharing
template< unsigned int TDim, unsigned int TNumNodes>
class EulerianDiffusionElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EulerianDiffusionElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    EulerianDiffusionElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    EulerianDiffusionElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~EulerianDiffusionElement() {};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<EulerianDiffusionElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }
    
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<EulerianDiffusionElement>(NewId, pGeom, pProperties);
    }


    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false); //false says not to preserve existing storage!!


//         noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);
//         noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

//         //Crank-Nicholson factor
//         const double cr_nk = 0.5;

        const double delta_t = rCurrentProcessInfo[DELTA_TIME];
        const double dt_inv = 1.0 / delta_t;
		const double lumping_factor = 1.00 / double(TNumNodes);

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        //getting data for the given geometry
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
        double Volume;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);

        //here we get all the variables we will need
        array_1d<double,TNumNodes> phi, phi_old, phi_convected;
        //bool active_convection=false;    // to kill some terms in case active_convection=false. For the moment it is inactive.
        //using only one Gauss Point for the material properties and volumetric heat flux
        double conductivity = 0.0;
        double specific_heat = 0.0;
        double density = 0.0;

	//storing locally the flags to avoid repeated check in the nodal loops
        const bool IsDefinedDensityVariable = my_settings->IsDefinedDensityVariable();
        const bool IsDefinedSpecificHeatVariableVariable = my_settings->IsDefinedSpecificHeatVariable();
        const bool IsDefinedDiffusionVariable = my_settings->IsDefinedDiffusionVariable();
        const bool IsDefinedProjectionVariable = my_settings->IsDefinedProjectionVariable();

	//if it is a convection diffusion problem, then the projection variable will exist and therefore we must use it instead of unknownVar(timestep n)
	//that is, to take the convection into account.

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            phi[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
            phi_old[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar,1);
	    if (IsDefinedProjectionVariable)
            	phi_convected[i] = GetGeometry()[i].FastGetSolutionStepValue((my_settings->GetProjectionVariable()),0);
	    else
		phi_convected[i] = phi_old[i];


//             dphi_dt[i] = dt_inv*(phi[i] - phi_old [i];

			if (IsDefinedDensityVariable)
			{
				const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
				density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
			}
			else
				density += 1.0;

			if (IsDefinedSpecificHeatVariableVariable)
			{
				const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();
				specific_heat += GetGeometry()[i].FastGetSolutionStepValue(rSpecificHeatVar);
			}
			else
				specific_heat += 1.0;

			if (IsDefinedDiffusionVariable)
			{
				const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
				conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
			}
			//if not, then the conductivity = 0
        }

        conductivity *= lumping_factor;
        density *= lumping_factor;
        specific_heat *= lumping_factor;
        //heat_flux *= lumping_factor;

        BoundedMatrix<double,TNumNodes, TNumNodes> aux1 = ZeroMatrix(TNumNodes, TNumNodes); //terms multiplying dphi/dt

        BoundedMatrix<double,TNumNodes, TNumNodes> Ncontainer;
        GetShapeFunctionsOnGauss(Ncontainer);
        for(unsigned int igauss=0; igauss<TDim+1; igauss++)
        {
            noalias(N) = row(Ncontainer,igauss);
            noalias(aux1) += outer_prod(N, N);
        }

        //mass matrix
        noalias(rLeftHandSideMatrix)  = (dt_inv*density*specific_heat)*aux1;
        noalias(rRightHandSideVector) = (dt_inv*density*specific_heat)*prod(aux1,phi_convected);
        //adding the diffusion
        noalias(rLeftHandSideMatrix)  += (conductivity *0.5 * prod(DN_DX, trans(DN_DX)))*static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) -= prod((conductivity *0.5 * prod(DN_DX, trans(DN_DX))),phi_convected)*static_cast<double>(TNumNodes) ;



        //take out the dirichlet part to finish computing the residual
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, phi);

        rRightHandSideVector *= Volume/static_cast<double>(TNumNodes);
        rLeftHandSideMatrix *= Volume/static_cast<double>(TNumNodes);

        KRATOS_CATCH("Error in Eulerian ConvDiff Element")
    }

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rRightHandSideVector.size() != TNumNodes) {
            rRightHandSideVector.resize(TNumNodes, false); // False says not to preserve existing storage!!
        }

        ConvectionDiffusionSettings::Pointer p_my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const auto& r_unknown_var = p_my_settings->GetUnknownVariable();

        // Getting data for the given geometry
        double volume;
        array_1d<double, TNumNodes > N;
        BoundedMatrix<double, TNumNodes, TDim > DN_DX;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, volume);

        // Here we get all the variables we will need
        // Using only one Gauss Point for the material properties and volumetric heat flux
        double density = 0.0;
        double conductivity = 0.0;
        double specific_heat = 0.0;
        array_1d<double,TNumNodes> phi, phi_old, phi_convected;

	    // Storing locally the flags to avoid repeated check in the nodal loops
        const bool is_defined_density_variable = p_my_settings->IsDefinedDensityVariable();
        const bool is_defined_specific_heat_variable = p_my_settings->IsDefinedSpecificHeatVariable();
        const bool is_defined_diffusion_variable = p_my_settings->IsDefinedDiffusionVariable();
        const bool is_defined_projection_variable = p_my_settings->IsDefinedProjectionVariable();

        // Get nodal data
        const auto& r_geom = GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; i++) {
            const auto& r_node = r_geom[i];
            phi[i] = r_node.FastGetSolutionStepValue(r_unknown_var);

            // If it is a convection diffusion problem, then the projection variable will exist and 
            // Therefore we must use it instead of UnknownVariable(timestep n), that is, to take the convection into account.
	        if (is_defined_projection_variable) {
                const auto& r_projection_var = p_my_settings->GetProjectionVariable();
            	phi_convected[i] = r_node.FastGetSolutionStepValue(r_projection_var);
            } else {
                phi_old[i] = r_node.FastGetSolutionStepValue(r_unknown_var,1);
		        phi_convected[i] = phi_old[i];
            }

			if (is_defined_density_variable) {
				const auto& r_density_var = p_my_settings->GetDensityVariable();
				density += r_node.FastGetSolutionStepValue(r_density_var);
			} else {
				density += 1.0;
            }

			if (is_defined_specific_heat_variable) {
				const auto& r_specific_heat_var = p_my_settings->GetSpecificHeatVariable();
				specific_heat += r_node.FastGetSolutionStepValue(r_specific_heat_var);
			} else {
				specific_heat += 1.0;
            }

			if (is_defined_diffusion_variable) {
				const auto& r_diffusion_var = p_my_settings->GetDiffusionVariable();
				conductivity += r_node.FastGetSolutionStepValue(r_diffusion_var);
			} // If not, the conductivity is 0
        }

        const double lumping_factor = 1.0 / static_cast<double>(TNumNodes);
        density *= lumping_factor;
        conductivity *= lumping_factor;
        specific_heat *= lumping_factor;

        BoundedMatrix<double, TNumNodes, TNumNodes> aux_NxN = ZeroMatrix(TNumNodes, TNumNodes); // Terms multiplying dphi/dt
        BoundedMatrix<double, TNumNodes, TNumNodes> N_container;
        GetShapeFunctionsOnGauss(N_container);
        for(unsigned int i_gauss=0; i_gauss < TDim+1; ++i_gauss){
            noalias(N) = row(N_container, i_gauss);
            noalias(aux_NxN) += outer_prod(N, N);
        }

        const double dt_inv = 1.0 / rCurrentProcessInfo[DELTA_TIME];

        // Mass matrix
        const double aux_1 = dt_inv * density * specific_heat * volume / static_cast<double>(TNumNodes);
        noalias(rRightHandSideVector) = aux_1 * prod(aux_NxN, phi_convected - phi);
        // Adding the diffusion
        // Note 1: the diffusive term is computed using a Crank-Nicholson scheme
        // Note 2: the gradients already include the corresponding Gauss point weight (no need to divide by TNumNodes)
        const double aux_2 = conductivity * 0.5 * volume;
        noalias(rRightHandSideVector) -= aux_2 * prod(prod(DN_DX, trans(DN_DX)), phi_convected + phi);

        KRATOS_CATCH("Error in Eulerian diffusion element CalculateRightHandSide")
    }

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
        }
        KRATOS_CATCH("")

    }





    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != TNumNodes)
            ElementalDofList.resize(TNumNodes);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

        }
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
        return "LevelSetConvectionElementSimplex #";
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

    EulerianDiffusionElement() : Element()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{

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
    //         ASGS2D() : Element()
    //         {
    //         }

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

#endif // KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED  defined

