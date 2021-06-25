//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla (adapted by Andrea Montanino)
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "custom_elements/compressible_ns_biphase_explicit.h"


namespace Kratos {

unsigned int	conta_new(
    unsigned int i, 
    unsigned int j, 
    unsigned int k, 
    unsigned int l, 
    unsigned int maxi, 
    unsigned int maxj, 
    unsigned int maxk, 
    unsigned int maxl)
{
	
	unsigned int s;
	
	s = l + maxl*k + maxl*maxk*j + maxl*maxk*maxj*i;
	
	return s;
}


void ShockCapturing2d_new(const double mu,
               const double lambda,
               const double c_v,
			   const double ros,
               const double h,
			   array_1d<double,2>& dt_diff,
			   array_1d<double,2>& ds_diff,
               array_1d<double,4>& tau,
               array_1d<double,2>& q,
               double ro_el,
               array_1d<double,10>& gradU, // 5*2
               array_1d<double,5>& Residual,
			   array_1d<double,5>& U,
			   array_1d<double,2>& a_el,
			   double norm_u,
			   double SpeedSound)
{
    const int SpaceDimension  = 2;
	
    double alpha  	  = 0.2;
	double alpha_dc   = 0.0;
	double alpha_de	  = 0.2;
	const double tol  = 1e-32;
    const double tol2 = 1e-32;                               

    unsigned int i;

    double a_ds = 0.0;
	double a_dt = 0.0;
	double v_sc = 0.0;
    double k_sc = 0.0;
    
	array_1d<double,SpaceDimension>  res_m;
	array_1d<double,SpaceDimension>  t_el;
    double dt_ref 	= U(0);
	double ds_ref 	= U(1);
	double res_dt 	= Residual(0);
	double res_ds 	= Residual(1);
	double res_e 	= Residual(4);
    double norm_res_ds;
	double norm_res_dt;
	double norm_res_m;
    double norm_res_e;
    double normgradm = 0.0;

	double normgraddt = sqrt(gradU(0)*gradU(0) + gradU(1)*gradU(1));
	double normgradds = sqrt(gradU(2)*gradU(2) + gradU(3)*gradU(3));
	    
	t_el(0) = -a_el(1);
	t_el(1) =  a_el(0);

    res_m(0) = Residual(2); 
    res_m(1) = Residual(3);

	norm_res_ds = sqrt(res_ds*res_ds);
	norm_res_dt = sqrt(res_dt*res_dt);
    norm_res_m 	= sqrt(res_m(0)*res_m(0) + res_m(1)*res_m(1));
    norm_res_e 	= sqrt(res_e*res_e);

	for (i = 4; i < 8; i++){
        normgradm += gradU(i)*gradU(i);
    }

    normgradm = sqrt(normgradm);

    if (normgradm > tol){         
        v_sc = 0.5*h*alpha*(norm_res_m/normgradm);
    }
    
    
	double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (i = 8; i < 10; i++)      norm_grade += gradU(i)*gradU(i);
    
    norm_grade = sqrt(norm_grade);
    
	if (norm_grade > tol)         k_sc = 0.5*h*alpha_de*(norm_res_e/norm_grade);
	
	if (normgraddt > tol) 		  a_dt = 0.5*h*alpha_dc*norm_res_dt/normgraddt; 
	if (normgradds > tol) 		  a_ds = 0.5*h*alpha_dc*norm_res_ds/normgradds;

	
//	Isotropic DC density

	for (i = 0; i < 2; i++)       dt_diff[i] = a_dt*gradU[i];
	for (i = 0; i < 2; i++)       ds_diff[i] = a_ds*gradU[i+2];


//	Crosswind DC density

/*
	const double ads = a_el[0]*gradU[2] + a_el[1]*gradU[3];
	const double adt = a_el[0]*gradU[0] + a_el[1]*gradU[1];

	const double tds = t_el[0]*gradU[2] + t_el[1]*gradU[3];
	const double tdt = t_el[0]*gradU[0] + t_el[1]*gradU[1];

	const double k_dt_stream = 0.0;
	const double k_ds_stream = 1.0;

	const double k_dt_cross = 0.0;
	const double k_ds_cross = 0.5;

	for (i = 0; i < 2; i++){
		dt_diff[i] = k_dt_stream*a_dt*adt*a_el[i] + k_dt_cross*a_dt*tdt*t_el[i];
		ds_diff[i] = k_ds_stream*a_ds*ads*a_el[i] + k_ds_cross*a_ds*tds*t_el[i];
	}
*/

	for (i = 0; i < 4; i++)       tau[i] *= (1.0 + ro_el*v_sc/mu);

    for (i = 0; i < 2; i++)       q[i] *= (1.0 + ro_el*c_v*k_sc/lambda);


}

void ShockCapturing3d_new(const double mu,
               const double lambda,
               const double c_v,
			   const double ros,
               const double h,
			   array_1d<double,3>& dt_diff,
			   array_1d<double,3>& ds_diff,
               array_1d<double,9>& tau,
               array_1d<double,3>& q,
               double ro_el,
               array_1d<double,18>& gradU, // 6*3
               array_1d<double,6>& Residual,
			   array_1d<double,6>& U,
			   array_1d<double,3>& a_el,
			   double norm_u,
			   double SpeedSound)
{
    const int SpaceDimension  = 3;
	
    double alpha  	  = 0.5;
	double alpha_dc   = 1.0;
	double alpha_de	  = 1.0;
	const double tol  = 1e-32;
    const double tol2 = 1e-32;                               

    unsigned int i;

    double a_dt = 0.0;
	double a_ds = 0.0;
	double v_sc = 0.0;
    double k_sc = 0.0;
    
	array_1d<double,SpaceDimension>  res_m;
	array_1d<double,SpaceDimension>  t_el;
    array_1d<double,SpaceDimension>  v_el;
    double dt_ref 	= U(0);
	double ds_ref 	= U(1);
	double res_dt 	= Residual(0);
	double res_ds 	= Residual(1);
	double res_e 	= Residual(5);
    double norm_res_ds;
	double norm_res_dt;
	double norm_res_m;
    double norm_res_e;
    double normgradm = 0.0;

	double normgraddt = sqrt(gradU(0)*gradU(0) + gradU(1)*gradU(1) + gradU(2)*gradU(2));
	double normgradds = sqrt(gradU(3)*gradU(3) + gradU(4)*gradU(4) + gradU(5)*gradU(5));
	    
	t_el(0) = 1 - a_el(0)*a_el(0);
	t_el(1) =   - a_el(0)*a_el(1);
    t_el(2) =   - a_el(0)*a_el(2);

    v_el(0) =   - a_el(1)*a_el(0);
	v_el(1) = 1 - a_el(1)*a_el(1);
    v_el(2) =   - a_el(1)*a_el(2);

    res_m(0) = Residual(2); 
    res_m(1) = Residual(3);
    res_m(2) = Residual(4);

	norm_res_ds = sqrt(res_ds*res_ds);
	norm_res_dt = sqrt(res_dt*res_dt);
    norm_res_m 	= sqrt(res_m(0)*res_m(0) + res_m(1)*res_m(1) + res_m(2)*res_m(2));
    norm_res_e 	= sqrt(res_e*res_e);

	for (i = 6; i < 15; i++){
        normgradm += gradU(i)*gradU(i);
    }

    normgradm = sqrt(normgradm);

    if (normgradm > tol){         
        v_sc = 0.5*h*alpha*(norm_res_m/normgradm);
    }
    
    
	double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (i = 15; i < 18; i++)      norm_grade += gradU(i)*gradU(i);
    
    norm_grade = sqrt(norm_grade);
    
	if (norm_grade > tol)         k_sc = 0.5*h*alpha_de*(norm_res_e/norm_grade);
	
	if (normgraddt > tol) 		  a_dt = 0.5*h*alpha_dc*norm_res_dt/normgraddt; 
	if (normgradds > tol) 		  a_ds = 0.5*h*alpha_dc*norm_res_ds/normgradds;

	
//	Isotropic DC density

	for (i = 0; i < 3; i++)       dt_diff[i] = a_dt*gradU[i];
	for (i = 0; i < 3; i++)       ds_diff[i] = a_ds*gradU[i + 3];


//	Crosswind DC density

/*
	const double ads = a_el[0]*gradU[2] + a_el[1]*gradU[3];
	const double adt = a_el[0]*gradU[0] + a_el[1]*gradU[1];

	const double tds = t_el[0]*gradU[2] + t_el[1]*gradU[3];
	const double tdt = t_el[0]*gradU[0] + t_el[1]*gradU[1];

	const double k_dt_stream = 0.0;
	const double k_ds_stream = 1.0;

	const double k_dt_cross = 0.0;
	const double k_ds_cross = 0.5;

	for (i = 0; i < 2; i++){
		dt_diff[i] = k_dt_stream*a_dt*adt*a_el[i] + k_dt_cross*a_dt*tdt*t_el[i];
		ds_diff[i] = k_ds_stream*a_ds*ads*a_el[i] + k_ds_cross*a_ds*tds*t_el[i];
	}
*/

	for (i = 0; i < 9; i++)       tau[i] *= (1.0 + ro_el*v_sc/mu);

    for (i = 0; i < 3; i++)       q[i] *= (1.0 + ro_el*c_v*k_sc/lambda);


}

template <>
void CompressibleNSBiphaseExplicit<2,3>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 5;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (rResult.size() != dof_size) {
        rResult.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int den_sol_pos = r_geometry[0].GetDofPosition(DENSITY_SOLID);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY_SOLID, den_sol_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 6;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (rResult.size() != dof_size) {
        rResult.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int den_sol_pos = r_geometry[0].GetDofPosition(DENSITY_SOLID);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY, den_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(DENSITY_SOLID, den_sol_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_X, mom_pos).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Y, mom_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(MOMENTUM_Z, mom_pos + 2).EquationId();
        rResult[local_index++] = r_geometry[i_node].GetDof(TOTAL_ENERGY, enr_pos).EquationId();
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 5;
    constexpr unsigned int dof_size = n_nodes * block_size;

    if (ElementalDofList.size() != dof_size) {
        ElementalDofList.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto& r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int den_sol_pos = r_geometry[0].GetDofPosition(DENSITY_SOLID);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(DENSITY_SOLID, den_sol_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = r_geometry[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::GetDofList(
    DofsVectorType &ElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;
    constexpr unsigned int block_size = 6;
    unsigned int dof_size = n_nodes * block_size;

    if (ElementalDofList.size() != dof_size) {
        ElementalDofList.resize(dof_size);
    }

    unsigned int local_index = 0;
    const auto &r_geometry = GetGeometry();
    const unsigned int den_pos = r_geometry[0].GetDofPosition(DENSITY);
    const unsigned int den_sol_pos = r_geometry[0].GetDofPosition(DENSITY_SOLID);
    const unsigned int mom_pos = r_geometry[0].GetDofPosition(MOMENTUM);
    const unsigned int enr_pos = r_geometry[0].GetDofPosition(TOTAL_ENERGY);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DENSITY, den_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(DENSITY_SOLID, den_sol_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_X, mom_pos);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Y, mom_pos + 1);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(MOMENTUM_Z, mom_pos + 2);
        ElementalDofList[local_index++] = this->GetGeometry()[i_node].pGetDof(TOTAL_ENERGY, enr_pos);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>
int CompressibleNSBiphaseExplicit<TDim, TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Perform basic element checks
    int ErrorCode = Kratos::Element::Check(rCurrentProcessInfo);
    if (ErrorCode != 0) {
        return ErrorCode;
    }

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY)) << "Missing DENSITY variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY_SOLID)) << "Missing DENSITY_SOLID variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(MOMENTUM)) << "Missing MOMENTUM variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE)) << "Missing BODY_FORCE variable on solution step data for node " << this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].SolutionStepsDataHas(HEAT_SOURCE)) << "Missing HEAT_SOURCE variable on solution step data for node " << this->GetGeometry()[i].Id();

        // Activate as soon as we start using the explicit DOF based strategy
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY)) << "Missing DENSITY DOF in node ", this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(DENSITY_SOLID)) << "Missing DENSITY_SOLID DOF in node ", this->GetGeometry()[i].Id();
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_X) || this->GetGeometry()[i].HasDofFor(MOMENTUM_Y)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
        if (TDim == 3) {
            KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(MOMENTUM_Z)) << "Missing MOMENTUM component DOF in node ", this->GetGeometry()[i].Id();
        }
        KRATOS_ERROR_IF_NOT(this->GetGeometry()[i].HasDofFor(TOTAL_ENERGY)) << "Missing TOTAL_ENERGY DOF in node ", this->GetGeometry()[i].Id();
    }

    return 0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes>   // Verify if add DENSOLID_PROJ  
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::Calculate(
    const Variable<double>& rVariable,
    double& Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Lumped projection terms
    if (rVariable == DENSITY_PROJECTION) {
        CalculateDensityProjection(rCurrentProcessInfo);
    } else if (rVariable == TOTAL_ENERGY_PROJECTION) {
        CalculateTotalEnergyProjection(rCurrentProcessInfo);
    } else if (rVariable == VELOCITY_DIVERGENCE) {
        Output = CalculateMidPointVelocityDivergence();
    } else if (rVariable == SOUND_VELOCITY) {
        Output = CalculateMidPointSoundVelocity();
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::Calculate(
    const Variable<array_1d<double, 3 > >& rVariable,
    array_1d<double, 3 > & Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == DENSITY_GRADIENT) {
        Output = CalculateMidPointDensityGradient();
    } else if (rVariable == TEMPERATURE_GRADIENT) {
        Output = CalculateMidPointTemperatureGradient();
    } else if (rVariable == VELOCITY_ROTATIONAL) {
        Output = CalculateMidPointVelocityRotational();
    } else if (rVariable == MOMENTUM_PROJECTION) {
        CalculateMomentumProjection(rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}


template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix & Output,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VELOCITY_GRADIENT) {
        Output = CalculateMidPointVelocityGradient();
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize( r_integration_points.size() );
    }

    if (rVariable == SHOCK_SENSOR) {
        const double sc = this->GetValue(SHOCK_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == SHEAR_SENSOR) {
        const double sc = this->GetValue(SHEAR_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == THERMAL_SENSOR) {
        const double sc = this->GetValue(THERMAL_SENSOR);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = sc;
        }
    } else if (rVariable == ARTIFICIAL_CONDUCTIVITY) {
        const double k_star = this->GetValue(ARTIFICIAL_CONDUCTIVITY);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = k_star;
        }
    } else if (rVariable == ARTIFICIAL_BULK_VISCOSITY) {
        const double beta_star = this->GetValue(ARTIFICIAL_BULK_VISCOSITY);
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = beta_star;
        }
    } else if (rVariable == VELOCITY_DIVERGENCE) {
        const double div_v = CalculateMidPointVelocityDivergence();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = div_v;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double,3>>& rVariable,
    std::vector<array_1d<double,3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints();
    if (rOutput.size() != r_integration_points.size()) {
        rOutput.resize( r_integration_points.size() );
    }

    if (rVariable == DENSITY_GRADIENT) {
        const auto& rho_grad = CalculateMidPointDensityGradient();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = rho_grad;
        }
    } else if (rVariable == TEMPERATURE_GRADIENT) {
        const auto& rho_grad = CalculateMidPointTemperatureGradient();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = rho_grad;
        }
    } else if (rVariable == VELOCITY_ROTATIONAL) {
        const auto rot_v = CalculateMidPointVelocityRotational();
        for (unsigned int i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss) {
            rOutput[i_gauss] = rot_v;
        }
    } else {
        KRATOS_ERROR << "Variable not implemented." << std::endl;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{

    // Getting data for the given geometry
    const auto& r_geometry = GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, rData.DN_DX, rData.N, rData.volume);

    // Compute element size
    rData.h = ElementSizeCalculator<TDim, TNumNodes>::GradientsElementSize(rData.DN_DX);

    // Database access to all of the variables needed
    Properties &r_properties = this->GetProperties();
    rData.mu = r_properties.GetValue(DYNAMIC_VISCOSITY);
    rData.lambda = r_properties.GetValue(CONDUCTIVITY);
    rData.c_v = r_properties.GetValue(SPECIFIC_HEAT); // TODO: WE SHOULD SPECIFY WHICH ONE --> CREATE SPECIFIC_HEAT_CONSTANT_VOLUME
    rData.gamma = r_properties.GetValue(HEAT_CAPACITY_RATIO);
    rData.c_s = r_properties.GetValue(SOLID_MATERIAL_SPECIFIC_HEAT);
    rData.ros = r_properties.GetValue(SOLID_MATERIAL_DENSITY);

    rData.time_step = rCurrentProcessInfo[DELTA_TIME];

    rData.UseOSS = rCurrentProcessInfo[OSS_SWITCH];
    rData.ShockCapturing = rCurrentProcessInfo[SHOCK_CAPTURING_SWITCH];

    // Magnitudes to calculate the time derivatives
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    const double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    const double aux_theta = theta > 0 ? 1.0 / (theta * time_step) : 0.0;

    // Get nodal values
    if (rData.UseOSS) {
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            const auto& r_node = r_geometry[i];
            // Vector data
            const array_1d<double,3>& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double,3>& r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
            const array_1d<double,3>& r_momentum_projection = r_node.GetValue(MOMENTUM_PROJECTION);
            const array_1d<double,3> mom_inc = r_momentum - r_momentum_old;
            const auto& r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);
            for (unsigned int k = 0; k < TDim; ++k) {
                rData.U(i, k + 2) = r_momentum[k];
                rData.dUdt(i, k + 2) = aux_theta * mom_inc[k];
                rData.ResProj(i, k + 2) = r_momentum_projection[k];
                rData.f_ext(i, k) = r_body_force[k];
            }
            // Density data
            const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
            const double rho_inc = r_rho - r_rho_old;
            rData.U(i, 0) = r_rho;
            rData.dUdt(i, 0) = aux_theta * rho_inc;
            rData.ResProj(i, 0) = r_node.GetValue(DENSITY_PROJECTION);
            // Total energy data
            const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            const double tot_ener_inc = r_tot_ener - r_tot_ener_old;
            rData.U(i, TDim + 2) = r_tot_ener;
            rData.dUdt(i, TDim + 2) = aux_theta * tot_ener_inc;
            rData.ResProj(i, TDim + 2) = r_node.GetValue(TOTAL_ENERGY_PROJECTION);
            // Source data
            rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
            rData.m_ext(i) = r_node.FastGetSolutionStepValue(MASS_SOURCE);
            // Shock capturing data
            rData.mu_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
            rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
            rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
        }
    } else {
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            const auto& r_node = r_geometry[i];
            // Vector data
            const array_1d<double,3>& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
            const array_1d<double,3>& r_momentum_old = r_node.FastGetSolutionStepValue(MOMENTUM, 1);
            const array_1d<double,3> mom_inc = r_momentum - r_momentum_old;
            const auto& r_body_force = r_node.FastGetSolutionStepValue(BODY_FORCE);
            for (unsigned int k = 0; k < TDim; ++k) {
                rData.U(i, k + 2) = r_momentum[k];
                rData.dUdt(i, k + 2) = aux_theta * mom_inc[k];
                rData.f_ext(i, k) = r_body_force[k];
            }
            // Density data
            const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = r_node.FastGetSolutionStepValue(DENSITY, 1);
            rData.U(i, 0) = r_rho;
            rData.dUdt(i, 0) = aux_theta * (r_rho - r_rho_old);
            // Density solid data
            const double& r_rho_sol = r_node.FastGetSolutionStepValue(DENSITY_SOLID);
            const double& r_rho_sol_old = r_node.FastGetSolutionStepValue(DENSITY_SOLID, 1);
            rData.U(i, 1) = r_rho_sol;
            rData.dUdt(i, 1) = aux_theta * (r_rho_sol - r_rho_sol_old);
            // Total energy data
            const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_ener_old = r_node.FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            rData.U(i, TDim + 2) = r_tot_ener;
            rData.dUdt(i, TDim + 2) = aux_theta * (r_tot_ener - r_tot_ener_old);
            // Source data
            rData.r_ext(i) = r_node.FastGetSolutionStepValue(HEAT_SOURCE);
            rData.m_ext(i) = r_node.FastGetSolutionStepValue(MASS_SOURCE);
            // Shock capturing data
            rData.mu_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY);
            rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
            rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double,3> CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateMidPointDensityGradient() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    array_1d<double,3> midpoint_grad_rho = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_grad_rho[d1] += node_dNdX(d1) * r_rho;
        }
    }

    return midpoint_grad_rho;
}

// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
// TODO: IN HERE WE HAVE LINEARIZED THE DERIVATIVE... CHECK IF WE SHOULD PROPERLY COMPUTE IT
template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double,3> CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateMidPointTemperatureGradient() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    const double c_v = GetProperties()[SPECIFIC_HEAT];
    array_1d<double,3> midpoint_grad_temp = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
        const array_1d<double, 3> vel = r_mom / r_rho;
        const double temp = (r_tot_ener / r_rho - 0.5 * inner_prod(vel, vel)) / c_v;
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_grad_temp[d1] += node_dNdX(d1) * temp;
        }
    }

    return midpoint_grad_temp;
}

template <unsigned int TDim, unsigned int TNumNodes>
double CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateMidPointSoundVelocity() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_tot_ener = 0.0;
    array_1d<double,TDim> midpoint_mom = ZeroVector(TDim);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
        midpoint_rho += r_rho;
        midpoint_tot_ener += r_tot_ener;
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_mom[d1] += r_mom(d1);
        }
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;
    midpoint_tot_ener /= n_nodes;

    // Calculate midpoint speed of sound
    const auto& r_prop = GetProperties();
    const double c_v = r_prop.GetValue(SPECIFIC_HEAT);
    const double gamma = r_prop.GetValue(HEAT_CAPACITY_RATIO);
    const double temp = (midpoint_tot_ener / midpoint_rho - inner_prod(midpoint_mom, midpoint_mom) / (2 * std::pow(midpoint_rho, 2))) / c_v;
    double midpoint_c = std::sqrt(gamma * (gamma - 1.0) * c_v * temp);
    return midpoint_c;
}

template <unsigned int TDim, unsigned int TNumNodes>
double CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateMidPointVelocityDivergence() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_div_mom = 0.0;
    array_1d<double,TDim> midpoint_mom = ZeroVector(TDim);
    array_1d<double,TDim> midpoint_grad_rho = ZeroVector(TDim);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        for (unsigned int d1 = 0; d1 < TDim; ++d1) {
            midpoint_mom[d1] += r_mom(d1);
            midpoint_div_mom += node_dNdX(d1) * r_mom(d1);
            midpoint_grad_rho[d1] += node_dNdX(d1) * r_rho;
        }
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity divergence
    // Note that the formulation is written in conservative variables. Hence we do div(mom/rho).
    double midpoint_div_v = (midpoint_rho * midpoint_div_mom - inner_prod(midpoint_mom, midpoint_grad_rho)) / std::pow(midpoint_rho, 2);
    return midpoint_div_v;
}

template <>
array_1d<double,3> CompressibleNSBiphaseExplicit<2,3>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / std::pow(midpoint_rho, 2);
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / std::pow(midpoint_rho, 2);
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = 0.0;
    midpoint_rot_v[1] = 0.0;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
array_1d<double,3> CompressibleNSBiphaseExplicit<3,4>::CalculateMidPointVelocityRotational() const
{
    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmx_dz = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dz = 0.0;
    double midpoint_dmz_dx = 0.0;
    double midpoint_dmz_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    double midpoint_rho_dz = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmx_dz += r_mom[0] * node_dNdX[2];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dz += r_mom[1] * node_dNdX[2];
        midpoint_dmz_dx += r_mom[2] * node_dNdX[0];
        midpoint_dmz_dy += r_mom[2] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
        midpoint_rho_dz += r_rho * node_dNdX[2];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity rotational
    // Note that the formulation is written in conservative variables. Hence we do rot(mom/rho).
    const double rho_pow = std::pow(midpoint_rho, 2);
    const double dvz_dy = (midpoint_dmz_dy * midpoint_rho - midpoint_mom[2] * midpoint_rho_dy) / rho_pow;
    const double dvy_dz = (midpoint_dmy_dz * midpoint_rho - midpoint_mom[1] * midpoint_rho_dz) / rho_pow;
    const double dvx_dz = (midpoint_dmx_dz * midpoint_rho - midpoint_mom[0] * midpoint_rho_dz) / rho_pow;
    const double dvz_dx = (midpoint_dmz_dx * midpoint_rho - midpoint_mom[2] * midpoint_rho_dx) / rho_pow;
    const double dvy_dx = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx) / rho_pow;
    const double dvx_dy = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy) / rho_pow;
    array_1d<double,3> midpoint_rot_v;
    midpoint_rot_v[0] = dvz_dy - dvy_dz;
    midpoint_rot_v[1] = dvx_dz - dvz_dx;
    midpoint_rot_v[2] = dvy_dx - dvx_dy;
    return midpoint_rot_v;
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNSBiphaseExplicit<2, 3>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v = ZeroMatrix(3, 3);
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
BoundedMatrix<double, 3, 3> CompressibleNSBiphaseExplicit<3, 4>::CalculateMidPointVelocityGradient() const
{
    KRATOS_TRY

    // Get geometry data
    const auto& r_geom = GetGeometry();
    const unsigned int n_nodes = r_geom.PointsNumber();
    Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
    r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
    const auto& r_dNdX = dNdX_container[0];

    // Calculate midpoint magnitudes
    double midpoint_rho = 0.0;
    double midpoint_dmx_dx = 0.0;
    double midpoint_dmx_dy = 0.0;
    double midpoint_dmx_dz = 0.0;
    double midpoint_dmy_dx = 0.0;
    double midpoint_dmy_dy = 0.0;
    double midpoint_dmy_dz = 0.0;
    double midpoint_dmz_dx = 0.0;
    double midpoint_dmz_dy = 0.0;
    double midpoint_dmz_dz = 0.0;
    double midpoint_rho_dx = 0.0;
    double midpoint_rho_dy = 0.0;
    double midpoint_rho_dz = 0.0;
    array_1d<double,3> midpoint_mom = ZeroVector(3);
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto& r_node = r_geom[i_node];
        const auto node_dNdX = row(r_dNdX, i_node);
        const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        midpoint_rho += r_rho;
        midpoint_mom += r_mom;
        midpoint_dmx_dx += r_mom[0] * node_dNdX[0];
        midpoint_dmx_dy += r_mom[0] * node_dNdX[1];
        midpoint_dmx_dz += r_mom[0] * node_dNdX[2];
        midpoint_dmy_dx += r_mom[1] * node_dNdX[0];
        midpoint_dmy_dy += r_mom[1] * node_dNdX[1];
        midpoint_dmy_dz += r_mom[1] * node_dNdX[2];
        midpoint_dmz_dx += r_mom[2] * node_dNdX[0];
        midpoint_dmz_dy += r_mom[2] * node_dNdX[1];
        midpoint_dmz_dz += r_mom[2] * node_dNdX[2];
        midpoint_rho_dx += r_rho * node_dNdX[0];
        midpoint_rho_dy += r_rho * node_dNdX[1];
        midpoint_rho_dz += r_rho * node_dNdX[2];
    }
    midpoint_rho /= n_nodes;
    midpoint_mom /= n_nodes;

    // Calculate velocity gradient
    // Note that the formulation is written in conservative variables. Hence we do grad(mom/rho).
    BoundedMatrix<double, 3, 3> midpoint_grad_v;
    midpoint_grad_v(0,0) = (midpoint_dmx_dx * midpoint_rho - midpoint_mom[0] * midpoint_rho_dx);
    midpoint_grad_v(0,1) = (midpoint_dmx_dy * midpoint_rho - midpoint_mom[0] * midpoint_rho_dy);
    midpoint_grad_v(0,2) = (midpoint_dmx_dz * midpoint_rho - midpoint_mom[0] * midpoint_rho_dz);
    midpoint_grad_v(1,0) = (midpoint_dmy_dx * midpoint_rho - midpoint_mom[1] * midpoint_rho_dx);
    midpoint_grad_v(1,1) = (midpoint_dmy_dy * midpoint_rho - midpoint_mom[1] * midpoint_rho_dy);
    midpoint_grad_v(1,2) = (midpoint_dmy_dz * midpoint_rho - midpoint_mom[1] * midpoint_rho_dz);
    midpoint_grad_v(2,0) = (midpoint_dmz_dx * midpoint_rho - midpoint_mom[2] * midpoint_rho_dx);
    midpoint_grad_v(2,1) = (midpoint_dmz_dy * midpoint_rho - midpoint_mom[2] * midpoint_rho_dy);
    midpoint_grad_v(2,2) = (midpoint_dmz_dz * midpoint_rho - midpoint_mom[2] * midpoint_rho_dz);
    midpoint_grad_v /= std::pow(midpoint_rho, 2);

    return midpoint_grad_v;

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);

    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    // Calculate shock capturing values
    BoundedVector<double, 6> mom_proj;

    const double cmom_proj0 =             -0.25*dUdt_1_1;
const double cmom_proj1 =             0.166666666666667*U_0_0;
const double cmom_proj2 =             0.166666666666667*U_2_0;
const double cmom_proj3 =             0.666666666666667*U_1_0 + cmom_proj1 + cmom_proj2;
const double cmom_proj4 =             0.166666666666667*f_ext(0,0);
const double cmom_proj5 =             0.166666666666667*f_ext(2,0);
const double cmom_proj6 =             cmom_proj3*(cmom_proj4 + cmom_proj5 + 0.666666666666667*f_ext(1,0));
const double cmom_proj7 =             0.166666666666667*cmom_proj6;
const double cmom_proj8 =             0.166666666666667*U_0_1;
const double cmom_proj9 =             0.166666666666667*U_2_1;
const double cmom_proj10 =             0.666666666666667*U_1_1 + cmom_proj8 + cmom_proj9;
const double cmom_proj11 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2;
const double cmom_proj12 =             1.0/cmom_proj3;
const double cmom_proj13 =             cmom_proj10*cmom_proj11*cmom_proj12;
const double cmom_proj14 =             -0.166666666666667*cmom_proj13;
const double cmom_proj15 =             0.166666666666667*U_0_2;
const double cmom_proj16 =             0.166666666666667*U_2_2;
const double cmom_proj17 =             0.666666666666667*U_1_2 + cmom_proj15 + cmom_proj16;
const double cmom_proj18 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double cmom_proj19 =             cmom_proj12*cmom_proj17*cmom_proj18;
const double cmom_proj20 =             -0.166666666666667*cmom_proj19;
const double cmom_proj21 =             gamma - 1;
const double cmom_proj22 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double cmom_proj23 =             cmom_proj12*cmom_proj17*cmom_proj21*cmom_proj22;
const double cmom_proj24 =             0.166666666666667*cmom_proj23;
const double cmom_proj25 =             1.0*gamma - 3.0;
const double cmom_proj26 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1;
const double cmom_proj27 =             cmom_proj10*cmom_proj12*cmom_proj25*cmom_proj26;
const double cmom_proj28 =             0.166666666666667*cmom_proj27;
const double cmom_proj29 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double cmom_proj30 =             pow(cmom_proj3, -2);
const double cmom_proj31 =             cmom_proj10*cmom_proj17*cmom_proj29*cmom_proj30;
const double cmom_proj32 =             0.166666666666667*cmom_proj31;
const double cmom_proj33 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double cmom_proj34 =             pow(cmom_proj10, 2);
const double cmom_proj35 =             0.5*gamma - 0.5;
const double cmom_proj36 =             pow(cmom_proj17, 2);
const double cmom_proj37 =             cmom_proj35*(cmom_proj34 + cmom_proj36);
const double cmom_proj38 =             cmom_proj30*cmom_proj33*(-cmom_proj34 + cmom_proj37);
const double cmom_proj39 =             -0.166666666666667*cmom_proj38;
const double cmom_proj40 =             -0.25*dUdt_2_1;
const double cmom_proj41 =             cmom_proj21*(DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3);
const double cmom_proj42 =             -1.0*cmom_proj41;
const double cmom_proj43 =             0.166666666666667*U_1_0;
const double cmom_proj44 =             0.666666666666667*U_2_0 + cmom_proj1 + cmom_proj43;
const double cmom_proj45 =             0.166666666666667*f_ext(1,0);
const double cmom_proj46 =             cmom_proj44*(cmom_proj4 + cmom_proj45 + 0.666666666666667*f_ext(2,0));
const double cmom_proj47 =             0.166666666666667*cmom_proj46;
const double cmom_proj48 =             0.666666666666667*U_0_0 + cmom_proj2 + cmom_proj43;
const double cmom_proj49 =             cmom_proj48*(cmom_proj45 + cmom_proj5 + 0.666666666666667*f_ext(0,0));
const double cmom_proj50 =             0.166666666666667*U_1_1;
const double cmom_proj51 =             0.666666666666667*U_2_1 + cmom_proj50 + cmom_proj8;
const double cmom_proj52 =             1.0/cmom_proj44;
const double cmom_proj53 =             cmom_proj11*cmom_proj51*cmom_proj52;
const double cmom_proj54 =             -0.166666666666667*cmom_proj53;
const double cmom_proj55 =             0.166666666666667*U_1_2;
const double cmom_proj56 =             0.666666666666667*U_2_2 + cmom_proj15 + cmom_proj55;
const double cmom_proj57 =             cmom_proj18*cmom_proj52*cmom_proj56;
const double cmom_proj58 =             -0.166666666666667*cmom_proj57;
const double cmom_proj59 =             0.666666666666667*U_0_1 + cmom_proj50 + cmom_proj9;
const double cmom_proj60 =             1.0/cmom_proj48;
const double cmom_proj61 =             cmom_proj11*cmom_proj59*cmom_proj60;
const double cmom_proj62 =             0.666666666666667*U_0_2 + cmom_proj16 + cmom_proj55;
const double cmom_proj63 =             cmom_proj18*cmom_proj60*cmom_proj62;
const double cmom_proj64 =             cmom_proj21*cmom_proj22*cmom_proj52*cmom_proj56;
const double cmom_proj65 =             0.166666666666667*cmom_proj64;
const double cmom_proj66 =             cmom_proj21*cmom_proj22*cmom_proj60*cmom_proj62;
const double cmom_proj67 =             cmom_proj25*cmom_proj26*cmom_proj51*cmom_proj52;
const double cmom_proj68 =             0.166666666666667*cmom_proj67;
const double cmom_proj69 =             cmom_proj25*cmom_proj26*cmom_proj59*cmom_proj60;
const double cmom_proj70 =             pow(cmom_proj44, -2);
const double cmom_proj71 =             cmom_proj29*cmom_proj51*cmom_proj56*cmom_proj70;
const double cmom_proj72 =             0.166666666666667*cmom_proj71;
const double cmom_proj73 =             pow(cmom_proj48, -2);
const double cmom_proj74 =             cmom_proj29*cmom_proj59*cmom_proj62*cmom_proj73;
const double cmom_proj75 =             pow(cmom_proj51, 2);
const double cmom_proj76 =             pow(cmom_proj56, 2);
const double cmom_proj77 =             cmom_proj35*(cmom_proj75 + cmom_proj76);
const double cmom_proj78 =             cmom_proj33*cmom_proj70*(-cmom_proj75 + cmom_proj77);
const double cmom_proj79 =             -0.166666666666667*cmom_proj78;
const double cmom_proj80 =             pow(cmom_proj59, 2);
const double cmom_proj81 =             pow(cmom_proj62, 2);
const double cmom_proj82 =             cmom_proj35*(cmom_proj80 + cmom_proj81);
const double cmom_proj83 =             cmom_proj33*cmom_proj73*(-cmom_proj80 + cmom_proj82);
const double cmom_proj84 =             -0.25*dUdt_1_2;
const double cmom_proj85 =             0.166666666666667*f_ext(0,1);
const double cmom_proj86 =             0.166666666666667*f_ext(2,1);
const double cmom_proj87 =             cmom_proj3*(cmom_proj85 + cmom_proj86 + 0.666666666666667*f_ext(1,1));
const double cmom_proj88 =             0.166666666666667*cmom_proj87;
const double cmom_proj89 =             cmom_proj10*cmom_proj12*cmom_proj22;
const double cmom_proj90 =             -0.166666666666667*cmom_proj89;
const double cmom_proj91 =             cmom_proj12*cmom_proj17*cmom_proj26;
const double cmom_proj92 =             -0.166666666666667*cmom_proj91;
const double cmom_proj93 =             cmom_proj10*cmom_proj12*cmom_proj18*cmom_proj21;
const double cmom_proj94 =             0.166666666666667*cmom_proj93;
const double cmom_proj95 =             cmom_proj11*cmom_proj12*cmom_proj17*cmom_proj25;
const double cmom_proj96 =             0.166666666666667*cmom_proj95;
const double cmom_proj97 =             cmom_proj10*cmom_proj17*cmom_proj30*cmom_proj33;
const double cmom_proj98 =             0.166666666666667*cmom_proj97;
const double cmom_proj99 =             cmom_proj29*cmom_proj30*(-cmom_proj36 + cmom_proj37);
const double cmom_proj100 =             -0.166666666666667*cmom_proj99;
const double cmom_proj101 =             -0.25*dUdt_2_2;
const double cmom_proj102 =             cmom_proj21*(DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3);
const double cmom_proj103 =             -1.0*cmom_proj102;
const double cmom_proj104 =             0.166666666666667*f_ext(1,1);
const double cmom_proj105 =             cmom_proj44*(cmom_proj104 + cmom_proj85 + 0.666666666666667*f_ext(2,1));
const double cmom_proj106 =             0.166666666666667*cmom_proj105;
const double cmom_proj107 =             cmom_proj48*(cmom_proj104 + cmom_proj86 + 0.666666666666667*f_ext(0,1));
const double cmom_proj108 =             cmom_proj22*cmom_proj51*cmom_proj52;
const double cmom_proj109 =             -0.166666666666667*cmom_proj108;
const double cmom_proj110 =             cmom_proj26*cmom_proj52*cmom_proj56;
const double cmom_proj111 =             -0.166666666666667*cmom_proj110;
const double cmom_proj112 =             cmom_proj22*cmom_proj59*cmom_proj60;
const double cmom_proj113 =             cmom_proj26*cmom_proj60*cmom_proj62;
const double cmom_proj114 =             cmom_proj18*cmom_proj21*cmom_proj51*cmom_proj52;
const double cmom_proj115 =             0.166666666666667*cmom_proj114;
const double cmom_proj116 =             cmom_proj18*cmom_proj21*cmom_proj59*cmom_proj60;
const double cmom_proj117 =             cmom_proj11*cmom_proj25*cmom_proj52*cmom_proj56;
const double cmom_proj118 =             0.166666666666667*cmom_proj117;
const double cmom_proj119 =             cmom_proj11*cmom_proj25*cmom_proj60*cmom_proj62;
const double cmom_proj120 =             cmom_proj33*cmom_proj51*cmom_proj56*cmom_proj70;
const double cmom_proj121 =             0.166666666666667*cmom_proj120;
const double cmom_proj122 =             cmom_proj33*cmom_proj59*cmom_proj62*cmom_proj73;
const double cmom_proj123 =             cmom_proj29*cmom_proj70*(-cmom_proj76 + cmom_proj77);
const double cmom_proj124 =             -0.166666666666667*cmom_proj123;
const double cmom_proj125 =             cmom_proj29*cmom_proj73*(-cmom_proj81 + cmom_proj82);
const double cmom_proj126 =             0.166666666666667*cmom_proj49 - 0.166666666666667*cmom_proj61 - 0.166666666666667*cmom_proj63 + 0.166666666666667*cmom_proj66 + 0.166666666666667*cmom_proj69 + 0.166666666666667*cmom_proj74 - 0.166666666666667*cmom_proj83 - 0.25*dUdt_0_1;
const double cmom_proj127 =             0.166666666666667*cmom_proj107 - 0.166666666666667*cmom_proj112 - 0.166666666666667*cmom_proj113 + 0.166666666666667*cmom_proj116 + 0.166666666666667*cmom_proj119 + 0.166666666666667*cmom_proj122 - 0.166666666666667*cmom_proj125 - 0.25*dUdt_0_2;
            mom_proj[0]=cmom_proj0 + cmom_proj14 + cmom_proj20 + cmom_proj24 + cmom_proj28 + cmom_proj32 + cmom_proj39 + cmom_proj40 + cmom_proj42 + cmom_proj47 + 0.666666666666667*cmom_proj49 + cmom_proj54 + cmom_proj58 - 0.666666666666667*cmom_proj61 - 0.666666666666667*cmom_proj63 + cmom_proj65 + 0.666666666666667*cmom_proj66 + cmom_proj68 + 0.666666666666667*cmom_proj69 + cmom_proj7 + cmom_proj72 + 0.666666666666667*cmom_proj74 + cmom_proj79 - 0.666666666666667*cmom_proj83 - 0.5*dUdt_0_1;
            mom_proj[1]=cmom_proj100 + cmom_proj101 + cmom_proj103 + cmom_proj106 + 0.666666666666667*cmom_proj107 + cmom_proj109 + cmom_proj111 - 0.666666666666667*cmom_proj112 - 0.666666666666667*cmom_proj113 + cmom_proj115 + 0.666666666666667*cmom_proj116 + cmom_proj118 + 0.666666666666667*cmom_proj119 + cmom_proj121 + 0.666666666666667*cmom_proj122 + cmom_proj124 - 0.666666666666667*cmom_proj125 + cmom_proj84 + cmom_proj88 + cmom_proj90 + cmom_proj92 + cmom_proj94 + cmom_proj96 + cmom_proj98 - 0.5*dUdt_0_2;
            mom_proj[2]=cmom_proj126 - 0.666666666666667*cmom_proj13 - 0.666666666666667*cmom_proj19 + 0.666666666666667*cmom_proj23 + 0.666666666666667*cmom_proj27 + 0.666666666666667*cmom_proj31 - 0.666666666666667*cmom_proj38 + cmom_proj40 + cmom_proj42 + cmom_proj47 + cmom_proj54 + cmom_proj58 + 0.666666666666667*cmom_proj6 + cmom_proj65 + cmom_proj68 + cmom_proj72 + cmom_proj79 - 0.5*dUdt_1_1;
            mom_proj[3]=cmom_proj101 + cmom_proj103 + cmom_proj106 + cmom_proj109 + cmom_proj111 + cmom_proj115 + cmom_proj118 + cmom_proj121 + cmom_proj124 + cmom_proj127 + 0.666666666666667*cmom_proj87 - 0.666666666666667*cmom_proj89 - 0.666666666666667*cmom_proj91 + 0.666666666666667*cmom_proj93 + 0.666666666666667*cmom_proj95 + 0.666666666666667*cmom_proj97 - 0.666666666666667*cmom_proj99 - 0.5*dUdt_1_2;
            mom_proj[4]=cmom_proj0 + cmom_proj126 + cmom_proj14 + cmom_proj20 + cmom_proj24 + cmom_proj28 + cmom_proj32 + cmom_proj39 - 1.0*cmom_proj41 + 0.666666666666667*cmom_proj46 - 0.666666666666667*cmom_proj53 - 0.666666666666667*cmom_proj57 + 0.666666666666667*cmom_proj64 + 0.666666666666667*cmom_proj67 + cmom_proj7 + 0.666666666666667*cmom_proj71 - 0.666666666666667*cmom_proj78 - 0.5*dUdt_2_1;
            mom_proj[5]=cmom_proj100 - 1.0*cmom_proj102 + 0.666666666666667*cmom_proj105 - 0.666666666666667*cmom_proj108 - 0.666666666666667*cmom_proj110 + 0.666666666666667*cmom_proj114 + 0.666666666666667*cmom_proj117 + 0.666666666666667*cmom_proj120 - 0.666666666666667*cmom_proj123 + cmom_proj127 + cmom_proj84 + cmom_proj88 + cmom_proj90 + cmom_proj92 + cmom_proj94 + cmom_proj96 + cmom_proj98 - 0.5*dUdt_2_2;

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    mom_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * dim;
        auto& r_mom_proj = r_geometry[i_node].GetValue(MOMENTUM_PROJECTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom_proj[d] += mom_proj[aux + d];
        }
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::CalculateMomentumProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_0_4 = data.U(0, 4);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_1_4 = data.U(1, 4);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_2_4 = data.U(2, 4);
    const double &U_3_0 = data.U(3, 0);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);
    const double &U_3_4 = data.U(3, 4);

    const double &dUdt_0_1 = data.dUdt(0, 1);
    const double &dUdt_0_2 = data.dUdt(0, 2);
    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_1 = data.dUdt(1, 1);
    const double &dUdt_1_2 = data.dUdt(1, 2);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_1 = data.dUdt(2, 1);
    const double &dUdt_2_2 = data.dUdt(2, 2);
    const double &dUdt_2_3 = data.dUdt(2, 3);
    const double &dUdt_3_1 = data.dUdt(3, 1);
    const double &dUdt_3_2 = data.dUdt(3, 2);
    const double &dUdt_3_3 = data.dUdt(3, 3);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    // Calculate shock capturing values
    BoundedVector<double, 12> mom_proj;

    const double cmom_proj0 =             0.1381966*U_0_0;
const double cmom_proj1 =             0.1381966*U_2_0;
const double cmom_proj2 =             0.1381966*U_3_0;
const double cmom_proj3 =             cmom_proj1 + cmom_proj2;
const double cmom_proj4 =             0.5854102*U_1_0 + cmom_proj0 + cmom_proj3;
const double cmom_proj5 =             0.1381966*f_ext(0,0);
const double cmom_proj6 =             0.1381966*f_ext(2,0);
const double cmom_proj7 =             0.1381966*f_ext(3,0);
const double cmom_proj8 =             cmom_proj6 + cmom_proj7;
const double cmom_proj9 =             cmom_proj4*(cmom_proj5 + cmom_proj8 + 0.5854102*f_ext(1,0));
const double cmom_proj10 =             0.1381966*cmom_proj9;
const double cmom_proj11 =             0.1381966*U_0_1;
const double cmom_proj12 =             0.1381966*U_2_1;
const double cmom_proj13 =             0.1381966*U_3_1;
const double cmom_proj14 =             cmom_proj12 + cmom_proj13;
const double cmom_proj15 =             0.5854102*U_1_1 + cmom_proj11 + cmom_proj14;
const double cmom_proj16 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2 + DN_DX_3_1*U_3_2;
const double cmom_proj17 =             1.0/cmom_proj4;
const double cmom_proj18 =             cmom_proj15*cmom_proj16*cmom_proj17;
const double cmom_proj19 =             -0.1381966*cmom_proj18;
const double cmom_proj20 =             DN_DX_0_2*U_0_3 + DN_DX_1_2*U_1_3 + DN_DX_2_2*U_2_3 + DN_DX_3_2*U_3_3;
const double cmom_proj21 =             cmom_proj15*cmom_proj17*cmom_proj20;
const double cmom_proj22 =             -0.1381966*cmom_proj21;
const double cmom_proj23 =             0.1381966*U_0_2;
const double cmom_proj24 =             0.1381966*U_2_2;
const double cmom_proj25 =             0.1381966*U_3_2;
const double cmom_proj26 =             cmom_proj24 + cmom_proj25;
const double cmom_proj27 =             0.5854102*U_1_2 + cmom_proj23 + cmom_proj26;
const double cmom_proj28 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1 + DN_DX_3_1*U_3_1;
const double cmom_proj29 =             cmom_proj17*cmom_proj27*cmom_proj28;
const double cmom_proj30 =             -0.1381966*cmom_proj29;
const double cmom_proj31 =             0.1381966*U_0_3;
const double cmom_proj32 =             0.1381966*U_2_3;
const double cmom_proj33 =             0.1381966*U_3_3;
const double cmom_proj34 =             cmom_proj32 + cmom_proj33;
const double cmom_proj35 =             0.5854102*U_1_3 + cmom_proj31 + cmom_proj34;
const double cmom_proj36 =             DN_DX_0_2*U_0_1 + DN_DX_1_2*U_1_1 + DN_DX_2_2*U_2_1 + DN_DX_3_2*U_3_1;
const double cmom_proj37 =             cmom_proj17*cmom_proj35*cmom_proj36;
const double cmom_proj38 =             -0.1381966*cmom_proj37;
const double cmom_proj39 =             gamma - 1;
const double cmom_proj40 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2 + DN_DX_3_0*U_3_2;
const double cmom_proj41 =             cmom_proj17*cmom_proj27*cmom_proj39*cmom_proj40;
const double cmom_proj42 =             0.1381966*cmom_proj41;
const double cmom_proj43 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3 + DN_DX_3_0*U_3_3;
const double cmom_proj44 =             cmom_proj17*cmom_proj35*cmom_proj39*cmom_proj43;
const double cmom_proj45 =             0.1381966*cmom_proj44;
const double cmom_proj46 =             1.0*gamma;
const double cmom_proj47 =             cmom_proj46 - 3.0;
const double cmom_proj48 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1 + DN_DX_3_0*U_3_1;
const double cmom_proj49 =             cmom_proj15*cmom_proj17*cmom_proj47*cmom_proj48;
const double cmom_proj50 =             0.1381966*cmom_proj49;
const double cmom_proj51 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0 + DN_DX_3_1*U_3_0;
const double cmom_proj52 =             pow(cmom_proj4, -2);
const double cmom_proj53 =             cmom_proj15*cmom_proj27*cmom_proj51*cmom_proj52;
const double cmom_proj54 =             0.1381966*cmom_proj53;
const double cmom_proj55 =             DN_DX_0_2*U_0_0 + DN_DX_1_2*U_1_0 + DN_DX_2_2*U_2_0 + DN_DX_3_2*U_3_0;
const double cmom_proj56 =             cmom_proj15*cmom_proj35*cmom_proj52*cmom_proj55;
const double cmom_proj57 =             0.1381966*cmom_proj56;
const double cmom_proj58 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0 + DN_DX_3_0*U_3_0;
const double cmom_proj59 =             pow(cmom_proj15, 2);
const double cmom_proj60 =             0.5*gamma - 0.5;
const double cmom_proj61 =             pow(cmom_proj27, 2);
const double cmom_proj62 =             pow(cmom_proj35, 2);
const double cmom_proj63 =             cmom_proj60*(cmom_proj59 + cmom_proj61 + cmom_proj62);
const double cmom_proj64 =             cmom_proj52*cmom_proj58*(-cmom_proj59 + cmom_proj63);
const double cmom_proj65 =             -0.1381966*cmom_proj64;
const double cmom_proj66 =             -0.19999999899376*dUdt_2_1;
const double cmom_proj67 =             -0.19999999899376*dUdt_3_1;
const double cmom_proj68 =             cmom_proj46 - 1.0;
const double cmom_proj69 =             -cmom_proj68*(DN_DX_0_0*U_0_4 + DN_DX_1_0*U_1_4 + DN_DX_2_0*U_2_4 + DN_DX_3_0*U_3_4);
const double cmom_proj70 =             0.1381966*U_1_0;
const double cmom_proj71 =             cmom_proj0 + cmom_proj70;
const double cmom_proj72 =             0.5854102*U_3_0 + cmom_proj1 + cmom_proj71;
const double cmom_proj73 =             0.1381966*f_ext(1,0);
const double cmom_proj74 =             cmom_proj5 + cmom_proj73;
const double cmom_proj75 =             cmom_proj72*(cmom_proj6 + cmom_proj74 + 0.5854102*f_ext(3,0));
const double cmom_proj76 =             0.1381966*cmom_proj75;
const double cmom_proj77 =             0.5854102*U_2_0 + cmom_proj2 + cmom_proj71;
const double cmom_proj78 =             cmom_proj77*(cmom_proj7 + cmom_proj74 + 0.5854102*f_ext(2,0));
const double cmom_proj79 =             0.1381966*cmom_proj78;
const double cmom_proj80 =             0.5854102*U_0_0 + cmom_proj3 + cmom_proj70;
const double cmom_proj81 =             cmom_proj80*(cmom_proj73 + cmom_proj8 + 0.5854102*f_ext(0,0));
const double cmom_proj82 =             0.1381966*U_1_1;
const double cmom_proj83 =             cmom_proj11 + cmom_proj82;
const double cmom_proj84 =             0.5854102*U_3_1 + cmom_proj12 + cmom_proj83;
const double cmom_proj85 =             1.0/cmom_proj72;
const double cmom_proj86 =             cmom_proj16*cmom_proj84*cmom_proj85;
const double cmom_proj87 =             -0.1381966*cmom_proj86;
const double cmom_proj88 =             cmom_proj20*cmom_proj84*cmom_proj85;
const double cmom_proj89 =             -0.1381966*cmom_proj88;
const double cmom_proj90 =             0.1381966*U_1_2;
const double cmom_proj91 =             cmom_proj23 + cmom_proj90;
const double cmom_proj92 =             0.5854102*U_3_2 + cmom_proj24 + cmom_proj91;
const double cmom_proj93 =             cmom_proj28*cmom_proj85*cmom_proj92;
const double cmom_proj94 =             -0.1381966*cmom_proj93;
const double cmom_proj95 =             0.1381966*U_1_3;
const double cmom_proj96 =             cmom_proj31 + cmom_proj95;
const double cmom_proj97 =             0.5854102*U_3_3 + cmom_proj32 + cmom_proj96;
const double cmom_proj98 =             cmom_proj36*cmom_proj85*cmom_proj97;
const double cmom_proj99 =             -0.1381966*cmom_proj98;
const double cmom_proj100 =             0.5854102*U_2_1 + cmom_proj13 + cmom_proj83;
const double cmom_proj101 =             1.0/cmom_proj77;
const double cmom_proj102 =             cmom_proj100*cmom_proj101*cmom_proj16;
const double cmom_proj103 =             -0.1381966*cmom_proj102;
const double cmom_proj104 =             cmom_proj100*cmom_proj101*cmom_proj20;
const double cmom_proj105 =             -0.1381966*cmom_proj104;
const double cmom_proj106 =             0.5854102*U_2_2 + cmom_proj25 + cmom_proj91;
const double cmom_proj107 =             cmom_proj101*cmom_proj106*cmom_proj28;
const double cmom_proj108 =             -0.1381966*cmom_proj107;
const double cmom_proj109 =             0.5854102*U_2_3 + cmom_proj33 + cmom_proj96;
const double cmom_proj110 =             cmom_proj101*cmom_proj109*cmom_proj36;
const double cmom_proj111 =             -0.1381966*cmom_proj110;
const double cmom_proj112 =             0.5854102*U_0_1 + cmom_proj14 + cmom_proj82;
const double cmom_proj113 =             1.0/cmom_proj80;
const double cmom_proj114 =             cmom_proj112*cmom_proj113*cmom_proj16;
const double cmom_proj115 =             cmom_proj112*cmom_proj113*cmom_proj20;
const double cmom_proj116 =             0.5854102*U_0_2 + cmom_proj26 + cmom_proj90;
const double cmom_proj117 =             cmom_proj113*cmom_proj116*cmom_proj28;
const double cmom_proj118 =             0.5854102*U_0_3 + cmom_proj34 + cmom_proj95;
const double cmom_proj119 =             cmom_proj113*cmom_proj118*cmom_proj36;
const double cmom_proj120 =             cmom_proj39*cmom_proj40*cmom_proj85*cmom_proj92;
const double cmom_proj121 =             0.1381966*cmom_proj120;
const double cmom_proj122 =             cmom_proj39*cmom_proj43*cmom_proj85*cmom_proj97;
const double cmom_proj123 =             0.1381966*cmom_proj122;
const double cmom_proj124 =             cmom_proj101*cmom_proj106*cmom_proj39*cmom_proj40;
const double cmom_proj125 =             0.1381966*cmom_proj124;
const double cmom_proj126 =             cmom_proj101*cmom_proj109*cmom_proj39*cmom_proj43;
const double cmom_proj127 =             0.1381966*cmom_proj126;
const double cmom_proj128 =             cmom_proj113*cmom_proj116*cmom_proj39*cmom_proj40;
const double cmom_proj129 =             cmom_proj113*cmom_proj118*cmom_proj39*cmom_proj43;
const double cmom_proj130 =             cmom_proj47*cmom_proj48*cmom_proj84*cmom_proj85;
const double cmom_proj131 =             0.1381966*cmom_proj130;
const double cmom_proj132 =             cmom_proj100*cmom_proj101*cmom_proj47*cmom_proj48;
const double cmom_proj133 =             0.1381966*cmom_proj132;
const double cmom_proj134 =             cmom_proj112*cmom_proj113*cmom_proj47*cmom_proj48;
const double cmom_proj135 =             pow(cmom_proj72, -2);
const double cmom_proj136 =             cmom_proj135*cmom_proj51*cmom_proj84*cmom_proj92;
const double cmom_proj137 =             0.1381966*cmom_proj136;
const double cmom_proj138 =             cmom_proj135*cmom_proj55*cmom_proj84*cmom_proj97;
const double cmom_proj139 =             0.1381966*cmom_proj138;
const double cmom_proj140 =             pow(cmom_proj77, -2);
const double cmom_proj141 =             cmom_proj100*cmom_proj106*cmom_proj140*cmom_proj51;
const double cmom_proj142 =             0.1381966*cmom_proj141;
const double cmom_proj143 =             cmom_proj100*cmom_proj109*cmom_proj140*cmom_proj55;
const double cmom_proj144 =             0.1381966*cmom_proj143;
const double cmom_proj145 =             pow(cmom_proj80, -2);
const double cmom_proj146 =             cmom_proj112*cmom_proj116*cmom_proj145*cmom_proj51;
const double cmom_proj147 =             cmom_proj112*cmom_proj118*cmom_proj145*cmom_proj55;
const double cmom_proj148 =             pow(cmom_proj84, 2);
const double cmom_proj149 =             pow(cmom_proj92, 2);
const double cmom_proj150 =             pow(cmom_proj97, 2);
const double cmom_proj151 =             cmom_proj60*(cmom_proj148 + cmom_proj149 + cmom_proj150);
const double cmom_proj152 =             cmom_proj135*cmom_proj58*(-cmom_proj148 + cmom_proj151);
const double cmom_proj153 =             -0.1381966*cmom_proj152;
const double cmom_proj154 =             pow(cmom_proj100, 2);
const double cmom_proj155 =             pow(cmom_proj106, 2);
const double cmom_proj156 =             pow(cmom_proj109, 2);
const double cmom_proj157 =             cmom_proj60*(cmom_proj154 + cmom_proj155 + cmom_proj156);
const double cmom_proj158 =             cmom_proj140*cmom_proj58*(-cmom_proj154 + cmom_proj157);
const double cmom_proj159 =             -0.1381966*cmom_proj158;
const double cmom_proj160 =             pow(cmom_proj112, 2);
const double cmom_proj161 =             pow(cmom_proj116, 2);
const double cmom_proj162 =             pow(cmom_proj118, 2);
const double cmom_proj163 =             cmom_proj60*(cmom_proj160 + cmom_proj161 + cmom_proj162);
const double cmom_proj164 =             cmom_proj145*cmom_proj58*(-cmom_proj160 + cmom_proj163);
const double cmom_proj165 =             0.1381966*f_ext(0,1);
const double cmom_proj166 =             0.1381966*f_ext(2,1);
const double cmom_proj167 =             0.1381966*f_ext(3,1);
const double cmom_proj168 =             cmom_proj166 + cmom_proj167;
const double cmom_proj169 =             cmom_proj4*(cmom_proj165 + cmom_proj168 + 0.5854102*f_ext(1,1));
const double cmom_proj170 =             0.1381966*cmom_proj169;
const double cmom_proj171 =             cmom_proj15*cmom_proj17*cmom_proj40;
const double cmom_proj172 =             -0.1381966*cmom_proj171;
const double cmom_proj173 =             cmom_proj17*cmom_proj27*cmom_proj48;
const double cmom_proj174 =             -0.1381966*cmom_proj173;
const double cmom_proj175 =             cmom_proj17*cmom_proj20*cmom_proj27;
const double cmom_proj176 =             -0.1381966*cmom_proj175;
const double cmom_proj177 =             DN_DX_0_2*U_0_2 + DN_DX_1_2*U_1_2 + DN_DX_2_2*U_2_2 + DN_DX_3_2*U_3_2;
const double cmom_proj178 =             cmom_proj17*cmom_proj177*cmom_proj35;
const double cmom_proj179 =             -0.1381966*cmom_proj178;
const double cmom_proj180 =             cmom_proj15*cmom_proj17*cmom_proj28*cmom_proj39;
const double cmom_proj181 =             0.1381966*cmom_proj180;
const double cmom_proj182 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3 + DN_DX_3_1*U_3_3;
const double cmom_proj183 =             cmom_proj17*cmom_proj182*cmom_proj35*cmom_proj39;
const double cmom_proj184 =             0.1381966*cmom_proj183;
const double cmom_proj185 =             cmom_proj16*cmom_proj17*cmom_proj27*cmom_proj47;
const double cmom_proj186 =             0.1381966*cmom_proj185;
const double cmom_proj187 =             cmom_proj15*cmom_proj27*cmom_proj52*cmom_proj58;
const double cmom_proj188 =             0.1381966*cmom_proj187;
const double cmom_proj189 =             cmom_proj27*cmom_proj35*cmom_proj52*cmom_proj55;
const double cmom_proj190 =             0.1381966*cmom_proj189;
const double cmom_proj191 =             cmom_proj51*cmom_proj52*(-cmom_proj61 + cmom_proj63);
const double cmom_proj192 =             -0.1381966*cmom_proj191;
const double cmom_proj193 =             -0.19999999899376*dUdt_2_2;
const double cmom_proj194 =             -0.19999999899376*dUdt_3_2;
const double cmom_proj195 =             -cmom_proj68*(DN_DX_0_1*U_0_4 + DN_DX_1_1*U_1_4 + DN_DX_2_1*U_2_4 + DN_DX_3_1*U_3_4);
const double cmom_proj196 =             0.1381966*f_ext(1,1);
const double cmom_proj197 =             cmom_proj165 + cmom_proj196;
const double cmom_proj198 =             cmom_proj72*(cmom_proj166 + cmom_proj197 + 0.5854102*f_ext(3,1));
const double cmom_proj199 =             0.1381966*cmom_proj198;
const double cmom_proj200 =             cmom_proj77*(cmom_proj167 + cmom_proj197 + 0.5854102*f_ext(2,1));
const double cmom_proj201 =             0.1381966*cmom_proj200;
const double cmom_proj202 =             cmom_proj80*(cmom_proj168 + cmom_proj196 + 0.5854102*f_ext(0,1));
const double cmom_proj203 =             cmom_proj40*cmom_proj84*cmom_proj85;
const double cmom_proj204 =             -0.1381966*cmom_proj203;
const double cmom_proj205 =             cmom_proj48*cmom_proj85*cmom_proj92;
const double cmom_proj206 =             -0.1381966*cmom_proj205;
const double cmom_proj207 =             cmom_proj20*cmom_proj85*cmom_proj92;
const double cmom_proj208 =             -0.1381966*cmom_proj207;
const double cmom_proj209 =             cmom_proj177*cmom_proj85*cmom_proj97;
const double cmom_proj210 =             -0.1381966*cmom_proj209;
const double cmom_proj211 =             cmom_proj100*cmom_proj101*cmom_proj40;
const double cmom_proj212 =             -0.1381966*cmom_proj211;
const double cmom_proj213 =             cmom_proj101*cmom_proj106*cmom_proj48;
const double cmom_proj214 =             -0.1381966*cmom_proj213;
const double cmom_proj215 =             cmom_proj101*cmom_proj106*cmom_proj20;
const double cmom_proj216 =             -0.1381966*cmom_proj215;
const double cmom_proj217 =             cmom_proj101*cmom_proj109*cmom_proj177;
const double cmom_proj218 =             -0.1381966*cmom_proj217;
const double cmom_proj219 =             cmom_proj112*cmom_proj113*cmom_proj40;
const double cmom_proj220 =             cmom_proj113*cmom_proj116*cmom_proj48;
const double cmom_proj221 =             cmom_proj113*cmom_proj116*cmom_proj20;
const double cmom_proj222 =             cmom_proj113*cmom_proj118*cmom_proj177;
const double cmom_proj223 =             cmom_proj28*cmom_proj39*cmom_proj84*cmom_proj85;
const double cmom_proj224 =             0.1381966*cmom_proj223;
const double cmom_proj225 =             cmom_proj182*cmom_proj39*cmom_proj85*cmom_proj97;
const double cmom_proj226 =             0.1381966*cmom_proj225;
const double cmom_proj227 =             cmom_proj100*cmom_proj101*cmom_proj28*cmom_proj39;
const double cmom_proj228 =             0.1381966*cmom_proj227;
const double cmom_proj229 =             cmom_proj101*cmom_proj109*cmom_proj182*cmom_proj39;
const double cmom_proj230 =             0.1381966*cmom_proj229;
const double cmom_proj231 =             cmom_proj112*cmom_proj113*cmom_proj28*cmom_proj39;
const double cmom_proj232 =             cmom_proj113*cmom_proj118*cmom_proj182*cmom_proj39;
const double cmom_proj233 =             cmom_proj16*cmom_proj47*cmom_proj85*cmom_proj92;
const double cmom_proj234 =             0.1381966*cmom_proj233;
const double cmom_proj235 =             cmom_proj101*cmom_proj106*cmom_proj16*cmom_proj47;
const double cmom_proj236 =             0.1381966*cmom_proj235;
const double cmom_proj237 =             cmom_proj113*cmom_proj116*cmom_proj16*cmom_proj47;
const double cmom_proj238 =             cmom_proj135*cmom_proj58*cmom_proj84*cmom_proj92;
const double cmom_proj239 =             0.1381966*cmom_proj238;
const double cmom_proj240 =             cmom_proj135*cmom_proj55*cmom_proj92*cmom_proj97;
const double cmom_proj241 =             0.1381966*cmom_proj240;
const double cmom_proj242 =             cmom_proj100*cmom_proj106*cmom_proj140*cmom_proj58;
const double cmom_proj243 =             0.1381966*cmom_proj242;
const double cmom_proj244 =             cmom_proj106*cmom_proj109*cmom_proj140*cmom_proj55;
const double cmom_proj245 =             0.1381966*cmom_proj244;
const double cmom_proj246 =             cmom_proj112*cmom_proj116*cmom_proj145*cmom_proj58;
const double cmom_proj247 =             cmom_proj116*cmom_proj118*cmom_proj145*cmom_proj55;
const double cmom_proj248 =             cmom_proj135*cmom_proj51*(-cmom_proj149 + cmom_proj151);
const double cmom_proj249 =             -0.1381966*cmom_proj248;
const double cmom_proj250 =             cmom_proj140*cmom_proj51*(-cmom_proj155 + cmom_proj157);
const double cmom_proj251 =             -0.1381966*cmom_proj250;
const double cmom_proj252 =             cmom_proj145*cmom_proj51*(-cmom_proj161 + cmom_proj163);
const double cmom_proj253 =             0.1381966*f_ext(0,2);
const double cmom_proj254 =             0.1381966*f_ext(2,2);
const double cmom_proj255 =             0.1381966*f_ext(3,2);
const double cmom_proj256 =             cmom_proj254 + cmom_proj255;
const double cmom_proj257 =             cmom_proj4*(cmom_proj253 + cmom_proj256 + 0.5854102*f_ext(1,2));
const double cmom_proj258 =             0.1381966*cmom_proj257;
const double cmom_proj259 =             cmom_proj15*cmom_proj17*cmom_proj43;
const double cmom_proj260 =             -0.1381966*cmom_proj259;
const double cmom_proj261 =             cmom_proj17*cmom_proj182*cmom_proj27;
const double cmom_proj262 =             -0.1381966*cmom_proj261;
const double cmom_proj263 =             cmom_proj17*cmom_proj35*cmom_proj48;
const double cmom_proj264 =             -0.1381966*cmom_proj263;
const double cmom_proj265 =             cmom_proj16*cmom_proj17*cmom_proj35;
const double cmom_proj266 =             -0.1381966*cmom_proj265;
const double cmom_proj267 =             cmom_proj15*cmom_proj17*cmom_proj36*cmom_proj39;
const double cmom_proj268 =             0.1381966*cmom_proj267;
const double cmom_proj269 =             cmom_proj17*cmom_proj177*cmom_proj27*cmom_proj39;
const double cmom_proj270 =             0.1381966*cmom_proj269;
const double cmom_proj271 =             cmom_proj17*cmom_proj20*cmom_proj35*cmom_proj47;
const double cmom_proj272 =             0.1381966*cmom_proj271;
const double cmom_proj273 =             cmom_proj15*cmom_proj35*cmom_proj52*cmom_proj58;
const double cmom_proj274 =             0.1381966*cmom_proj273;
const double cmom_proj275 =             cmom_proj27*cmom_proj35*cmom_proj51*cmom_proj52;
const double cmom_proj276 =             0.1381966*cmom_proj275;
const double cmom_proj277 =             cmom_proj52*cmom_proj55*(-cmom_proj62 + cmom_proj63);
const double cmom_proj278 =             -0.1381966*cmom_proj277;
const double cmom_proj279 =             -0.19999999899376*dUdt_2_3;
const double cmom_proj280 =             -0.19999999899376*dUdt_3_3;
const double cmom_proj281 =             -cmom_proj68*(DN_DX_0_2*U_0_4 + DN_DX_1_2*U_1_4 + DN_DX_2_2*U_2_4 + DN_DX_3_2*U_3_4);
const double cmom_proj282 =             0.1381966*f_ext(1,2);
const double cmom_proj283 =             cmom_proj253 + cmom_proj282;
const double cmom_proj284 =             cmom_proj72*(cmom_proj254 + cmom_proj283 + 0.5854102*f_ext(3,2));
const double cmom_proj285 =             0.1381966*cmom_proj284;
const double cmom_proj286 =             cmom_proj77*(cmom_proj255 + cmom_proj283 + 0.5854102*f_ext(2,2));
const double cmom_proj287 =             0.1381966*cmom_proj286;
const double cmom_proj288 =             cmom_proj80*(cmom_proj256 + cmom_proj282 + 0.5854102*f_ext(0,2));
const double cmom_proj289 =             cmom_proj43*cmom_proj84*cmom_proj85;
const double cmom_proj290 =             -0.1381966*cmom_proj289;
const double cmom_proj291 =             cmom_proj182*cmom_proj85*cmom_proj92;
const double cmom_proj292 =             -0.1381966*cmom_proj291;
const double cmom_proj293 =             cmom_proj48*cmom_proj85*cmom_proj97;
const double cmom_proj294 =             -0.1381966*cmom_proj293;
const double cmom_proj295 =             cmom_proj16*cmom_proj85*cmom_proj97;
const double cmom_proj296 =             -0.1381966*cmom_proj295;
const double cmom_proj297 =             cmom_proj100*cmom_proj101*cmom_proj43;
const double cmom_proj298 =             -0.1381966*cmom_proj297;
const double cmom_proj299 =             cmom_proj101*cmom_proj106*cmom_proj182;
const double cmom_proj300 =             -0.1381966*cmom_proj299;
const double cmom_proj301 =             cmom_proj101*cmom_proj109*cmom_proj48;
const double cmom_proj302 =             -0.1381966*cmom_proj301;
const double cmom_proj303 =             cmom_proj101*cmom_proj109*cmom_proj16;
const double cmom_proj304 =             -0.1381966*cmom_proj303;
const double cmom_proj305 =             cmom_proj112*cmom_proj113*cmom_proj43;
const double cmom_proj306 =             cmom_proj113*cmom_proj116*cmom_proj182;
const double cmom_proj307 =             cmom_proj113*cmom_proj118*cmom_proj48;
const double cmom_proj308 =             cmom_proj113*cmom_proj118*cmom_proj16;
const double cmom_proj309 =             cmom_proj36*cmom_proj39*cmom_proj84*cmom_proj85;
const double cmom_proj310 =             0.1381966*cmom_proj309;
const double cmom_proj311 =             cmom_proj177*cmom_proj39*cmom_proj85*cmom_proj92;
const double cmom_proj312 =             0.1381966*cmom_proj311;
const double cmom_proj313 =             cmom_proj100*cmom_proj101*cmom_proj36*cmom_proj39;
const double cmom_proj314 =             0.1381966*cmom_proj313;
const double cmom_proj315 =             cmom_proj101*cmom_proj106*cmom_proj177*cmom_proj39;
const double cmom_proj316 =             0.1381966*cmom_proj315;
const double cmom_proj317 =             cmom_proj112*cmom_proj113*cmom_proj36*cmom_proj39;
const double cmom_proj318 =             cmom_proj113*cmom_proj116*cmom_proj177*cmom_proj39;
const double cmom_proj319 =             cmom_proj20*cmom_proj47*cmom_proj85*cmom_proj97;
const double cmom_proj320 =             0.1381966*cmom_proj319;
const double cmom_proj321 =             cmom_proj101*cmom_proj109*cmom_proj20*cmom_proj47;
const double cmom_proj322 =             0.1381966*cmom_proj321;
const double cmom_proj323 =             cmom_proj113*cmom_proj118*cmom_proj20*cmom_proj47;
const double cmom_proj324 =             cmom_proj135*cmom_proj58*cmom_proj84*cmom_proj97;
const double cmom_proj325 =             0.1381966*cmom_proj324;
const double cmom_proj326 =             cmom_proj135*cmom_proj51*cmom_proj92*cmom_proj97;
const double cmom_proj327 =             0.1381966*cmom_proj326;
const double cmom_proj328 =             cmom_proj100*cmom_proj109*cmom_proj140*cmom_proj58;
const double cmom_proj329 =             0.1381966*cmom_proj328;
const double cmom_proj330 =             cmom_proj106*cmom_proj109*cmom_proj140*cmom_proj51;
const double cmom_proj331 =             0.1381966*cmom_proj330;
const double cmom_proj332 =             cmom_proj112*cmom_proj118*cmom_proj145*cmom_proj58;
const double cmom_proj333 =             cmom_proj116*cmom_proj118*cmom_proj145*cmom_proj51;
const double cmom_proj334 =             cmom_proj135*cmom_proj55*(-cmom_proj150 + cmom_proj151);
const double cmom_proj335 =             -0.1381966*cmom_proj334;
const double cmom_proj336 =             cmom_proj140*cmom_proj55*(-cmom_proj156 + cmom_proj157);
const double cmom_proj337 =             -0.1381966*cmom_proj336;
const double cmom_proj338 =             cmom_proj145*cmom_proj55*(-cmom_proj162 + cmom_proj163);
const double cmom_proj339 =             0.1381966*cmom_proj81;
const double cmom_proj340 =             -0.1381966*cmom_proj114;
const double cmom_proj341 =             -0.1381966*cmom_proj115;
const double cmom_proj342 =             -0.1381966*cmom_proj117;
const double cmom_proj343 =             -0.1381966*cmom_proj119;
const double cmom_proj344 =             0.1381966*cmom_proj128;
const double cmom_proj345 =             0.1381966*cmom_proj129;
const double cmom_proj346 =             0.1381966*cmom_proj134;
const double cmom_proj347 =             0.1381966*cmom_proj146;
const double cmom_proj348 =             0.1381966*cmom_proj147;
const double cmom_proj349 =             -0.1381966*cmom_proj164;
const double cmom_proj350 =             0.1381966*cmom_proj202;
const double cmom_proj351 =             -0.1381966*cmom_proj219;
const double cmom_proj352 =             -0.1381966*cmom_proj220;
const double cmom_proj353 =             -0.1381966*cmom_proj221;
const double cmom_proj354 =             -0.1381966*cmom_proj222;
const double cmom_proj355 =             0.1381966*cmom_proj231;
const double cmom_proj356 =             0.1381966*cmom_proj232;
const double cmom_proj357 =             0.1381966*cmom_proj237;
const double cmom_proj358 =             0.1381966*cmom_proj246;
const double cmom_proj359 =             0.1381966*cmom_proj247;
const double cmom_proj360 =             -0.1381966*cmom_proj252;
const double cmom_proj361 =             0.1381966*cmom_proj288;
const double cmom_proj362 =             -0.1381966*cmom_proj305;
const double cmom_proj363 =             -0.1381966*cmom_proj306;
const double cmom_proj364 =             -0.1381966*cmom_proj307;
const double cmom_proj365 =             -0.1381966*cmom_proj308;
const double cmom_proj366 =             0.1381966*cmom_proj317;
const double cmom_proj367 =             0.1381966*cmom_proj318;
const double cmom_proj368 =             0.1381966*cmom_proj323;
const double cmom_proj369 =             0.1381966*cmom_proj332;
const double cmom_proj370 =             0.1381966*cmom_proj333;
const double cmom_proj371 =             -0.1381966*cmom_proj338;
const double cmom_proj372 =             cmom_proj10 + cmom_proj19 + cmom_proj22 + cmom_proj30 + cmom_proj339 + cmom_proj340 + cmom_proj341 + cmom_proj342 + cmom_proj343 + cmom_proj344 + cmom_proj345 + cmom_proj346 + cmom_proj347 + cmom_proj348 + cmom_proj349 + cmom_proj38 + cmom_proj42 + cmom_proj45 + cmom_proj50 + cmom_proj54 + cmom_proj57 + cmom_proj65 + cmom_proj69 - 0.19999999899376*dUdt_0_1 - 0.19999999899376*dUdt_1_1;
const double cmom_proj373 =             cmom_proj170 + cmom_proj172 + cmom_proj174 + cmom_proj176 + cmom_proj179 + cmom_proj181 + cmom_proj184 + cmom_proj186 + cmom_proj188 + cmom_proj190 + cmom_proj192 + cmom_proj195 + cmom_proj350 + cmom_proj351 + cmom_proj352 + cmom_proj353 + cmom_proj354 + cmom_proj355 + cmom_proj356 + cmom_proj357 + cmom_proj358 + cmom_proj359 + cmom_proj360 - 0.19999999899376*dUdt_0_2 - 0.19999999899376*dUdt_1_2;
const double cmom_proj374 =             cmom_proj258 + cmom_proj260 + cmom_proj262 + cmom_proj264 + cmom_proj266 + cmom_proj268 + cmom_proj270 + cmom_proj272 + cmom_proj274 + cmom_proj276 + cmom_proj278 + cmom_proj281 + cmom_proj361 + cmom_proj362 + cmom_proj363 + cmom_proj364 + cmom_proj365 + cmom_proj366 + cmom_proj367 + cmom_proj368 + cmom_proj369 + cmom_proj370 + cmom_proj371 - 0.19999999899376*dUdt_0_3 - 0.19999999899376*dUdt_1_3;
            mom_proj[0]=cmom_proj10 + cmom_proj103 + cmom_proj105 + cmom_proj108 + cmom_proj111 - 0.5854102*cmom_proj114 - 0.5854102*cmom_proj115 - 0.5854102*cmom_proj117 - 0.5854102*cmom_proj119 + cmom_proj121 + cmom_proj123 + cmom_proj125 + cmom_proj127 + 0.5854102*cmom_proj128 + 0.5854102*cmom_proj129 + cmom_proj131 + cmom_proj133 + 0.5854102*cmom_proj134 + cmom_proj137 + cmom_proj139 + cmom_proj142 + cmom_proj144 + 0.5854102*cmom_proj146 + 0.5854102*cmom_proj147 + cmom_proj153 + cmom_proj159 - 0.5854102*cmom_proj164 + cmom_proj19 + cmom_proj22 + cmom_proj30 + cmom_proj38 + cmom_proj42 + cmom_proj45 + cmom_proj50 + cmom_proj54 + cmom_proj57 + cmom_proj65 + cmom_proj66 + cmom_proj67 + cmom_proj69 + cmom_proj76 + cmom_proj79 + 0.5854102*cmom_proj81 + cmom_proj87 + cmom_proj89 + cmom_proj94 + cmom_proj99 - 0.40000000301872*dUdt_0_1 - 0.19999999899376*dUdt_1_1;
            mom_proj[1]=cmom_proj170 + cmom_proj172 + cmom_proj174 + cmom_proj176 + cmom_proj179 + cmom_proj181 + cmom_proj184 + cmom_proj186 + cmom_proj188 + cmom_proj190 + cmom_proj192 + cmom_proj193 + cmom_proj194 + cmom_proj195 + cmom_proj199 + cmom_proj201 + 0.5854102*cmom_proj202 + cmom_proj204 + cmom_proj206 + cmom_proj208 + cmom_proj210 + cmom_proj212 + cmom_proj214 + cmom_proj216 + cmom_proj218 - 0.5854102*cmom_proj219 - 0.5854102*cmom_proj220 - 0.5854102*cmom_proj221 - 0.5854102*cmom_proj222 + cmom_proj224 + cmom_proj226 + cmom_proj228 + cmom_proj230 + 0.5854102*cmom_proj231 + 0.5854102*cmom_proj232 + cmom_proj234 + cmom_proj236 + 0.5854102*cmom_proj237 + cmom_proj239 + cmom_proj241 + cmom_proj243 + cmom_proj245 + 0.5854102*cmom_proj246 + 0.5854102*cmom_proj247 + cmom_proj249 + cmom_proj251 - 0.5854102*cmom_proj252 - 0.40000000301872*dUdt_0_2 - 0.19999999899376*dUdt_1_2;
            mom_proj[2]=cmom_proj258 + cmom_proj260 + cmom_proj262 + cmom_proj264 + cmom_proj266 + cmom_proj268 + cmom_proj270 + cmom_proj272 + cmom_proj274 + cmom_proj276 + cmom_proj278 + cmom_proj279 + cmom_proj280 + cmom_proj281 + cmom_proj285 + cmom_proj287 + 0.5854102*cmom_proj288 + cmom_proj290 + cmom_proj292 + cmom_proj294 + cmom_proj296 + cmom_proj298 + cmom_proj300 + cmom_proj302 + cmom_proj304 - 0.5854102*cmom_proj305 - 0.5854102*cmom_proj306 - 0.5854102*cmom_proj307 - 0.5854102*cmom_proj308 + cmom_proj310 + cmom_proj312 + cmom_proj314 + cmom_proj316 + 0.5854102*cmom_proj317 + 0.5854102*cmom_proj318 + cmom_proj320 + cmom_proj322 + 0.5854102*cmom_proj323 + cmom_proj325 + cmom_proj327 + cmom_proj329 + cmom_proj331 + 0.5854102*cmom_proj332 + 0.5854102*cmom_proj333 + cmom_proj335 + cmom_proj337 - 0.5854102*cmom_proj338 - 0.40000000301872*dUdt_0_3 - 0.19999999899376*dUdt_1_3;
            mom_proj[3]=cmom_proj103 + cmom_proj105 + cmom_proj108 + cmom_proj111 + cmom_proj121 + cmom_proj123 + cmom_proj125 + cmom_proj127 + cmom_proj131 + cmom_proj133 + cmom_proj137 + cmom_proj139 + cmom_proj142 + cmom_proj144 + cmom_proj153 + cmom_proj159 - 0.5854102*cmom_proj18 - 0.5854102*cmom_proj21 - 0.5854102*cmom_proj29 + cmom_proj339 + cmom_proj340 + cmom_proj341 + cmom_proj342 + cmom_proj343 + cmom_proj344 + cmom_proj345 + cmom_proj346 + cmom_proj347 + cmom_proj348 + cmom_proj349 - 0.5854102*cmom_proj37 + 0.5854102*cmom_proj41 + 0.5854102*cmom_proj44 + 0.5854102*cmom_proj49 + 0.5854102*cmom_proj53 + 0.5854102*cmom_proj56 - 0.5854102*cmom_proj64 + cmom_proj66 + cmom_proj67 + cmom_proj69 + cmom_proj76 + cmom_proj79 + cmom_proj87 + cmom_proj89 + 0.5854102*cmom_proj9 + cmom_proj94 + cmom_proj99 - 0.19999999899376*dUdt_0_1 - 0.40000000301872*dUdt_1_1;
            mom_proj[4]=0.5854102*cmom_proj169 - 0.5854102*cmom_proj171 - 0.5854102*cmom_proj173 - 0.5854102*cmom_proj175 - 0.5854102*cmom_proj178 + 0.5854102*cmom_proj180 + 0.5854102*cmom_proj183 + 0.5854102*cmom_proj185 + 0.5854102*cmom_proj187 + 0.5854102*cmom_proj189 - 0.5854102*cmom_proj191 + cmom_proj193 + cmom_proj194 + cmom_proj195 + cmom_proj199 + cmom_proj201 + cmom_proj204 + cmom_proj206 + cmom_proj208 + cmom_proj210 + cmom_proj212 + cmom_proj214 + cmom_proj216 + cmom_proj218 + cmom_proj224 + cmom_proj226 + cmom_proj228 + cmom_proj230 + cmom_proj234 + cmom_proj236 + cmom_proj239 + cmom_proj241 + cmom_proj243 + cmom_proj245 + cmom_proj249 + cmom_proj251 + cmom_proj350 + cmom_proj351 + cmom_proj352 + cmom_proj353 + cmom_proj354 + cmom_proj355 + cmom_proj356 + cmom_proj357 + cmom_proj358 + cmom_proj359 + cmom_proj360 - 0.19999999899376*dUdt_0_2 - 0.40000000301872*dUdt_1_2;
            mom_proj[5]=0.5854102*cmom_proj257 - 0.5854102*cmom_proj259 - 0.5854102*cmom_proj261 - 0.5854102*cmom_proj263 - 0.5854102*cmom_proj265 + 0.5854102*cmom_proj267 + 0.5854102*cmom_proj269 + 0.5854102*cmom_proj271 + 0.5854102*cmom_proj273 + 0.5854102*cmom_proj275 - 0.5854102*cmom_proj277 + cmom_proj279 + cmom_proj280 + cmom_proj281 + cmom_proj285 + cmom_proj287 + cmom_proj290 + cmom_proj292 + cmom_proj294 + cmom_proj296 + cmom_proj298 + cmom_proj300 + cmom_proj302 + cmom_proj304 + cmom_proj310 + cmom_proj312 + cmom_proj314 + cmom_proj316 + cmom_proj320 + cmom_proj322 + cmom_proj325 + cmom_proj327 + cmom_proj329 + cmom_proj331 + cmom_proj335 + cmom_proj337 + cmom_proj361 + cmom_proj362 + cmom_proj363 + cmom_proj364 + cmom_proj365 + cmom_proj366 + cmom_proj367 + cmom_proj368 + cmom_proj369 + cmom_proj370 + cmom_proj371 - 0.19999999899376*dUdt_0_3 - 0.40000000301872*dUdt_1_3;
            mom_proj[6]=-0.5854102*cmom_proj102 - 0.5854102*cmom_proj104 - 0.5854102*cmom_proj107 - 0.5854102*cmom_proj110 + cmom_proj121 + cmom_proj123 + 0.5854102*cmom_proj124 + 0.5854102*cmom_proj126 + cmom_proj131 + 0.5854102*cmom_proj132 + cmom_proj137 + cmom_proj139 + 0.5854102*cmom_proj141 + 0.5854102*cmom_proj143 + cmom_proj153 - 0.5854102*cmom_proj158 + cmom_proj372 + cmom_proj67 + cmom_proj76 + 0.5854102*cmom_proj78 + cmom_proj87 + cmom_proj89 + cmom_proj94 + cmom_proj99 - 0.40000000301872*dUdt_2_1;
            mom_proj[7]=cmom_proj194 + cmom_proj199 + 0.5854102*cmom_proj200 + cmom_proj204 + cmom_proj206 + cmom_proj208 + cmom_proj210 - 0.5854102*cmom_proj211 - 0.5854102*cmom_proj213 - 0.5854102*cmom_proj215 - 0.5854102*cmom_proj217 + cmom_proj224 + cmom_proj226 + 0.5854102*cmom_proj227 + 0.5854102*cmom_proj229 + cmom_proj234 + 0.5854102*cmom_proj235 + cmom_proj239 + cmom_proj241 + 0.5854102*cmom_proj242 + 0.5854102*cmom_proj244 + cmom_proj249 - 0.5854102*cmom_proj250 + cmom_proj373 - 0.40000000301872*dUdt_2_2;
            mom_proj[8]=cmom_proj280 + cmom_proj285 + 0.5854102*cmom_proj286 + cmom_proj290 + cmom_proj292 + cmom_proj294 + cmom_proj296 - 0.5854102*cmom_proj297 - 0.5854102*cmom_proj299 - 0.5854102*cmom_proj301 - 0.5854102*cmom_proj303 + cmom_proj310 + cmom_proj312 + 0.5854102*cmom_proj313 + 0.5854102*cmom_proj315 + cmom_proj320 + 0.5854102*cmom_proj321 + cmom_proj325 + cmom_proj327 + 0.5854102*cmom_proj328 + 0.5854102*cmom_proj330 + cmom_proj335 - 0.5854102*cmom_proj336 + cmom_proj374 - 0.40000000301872*dUdt_2_3;
            mom_proj[9]=cmom_proj103 + cmom_proj105 + cmom_proj108 + cmom_proj111 + 0.5854102*cmom_proj120 + 0.5854102*cmom_proj122 + cmom_proj125 + cmom_proj127 + 0.5854102*cmom_proj130 + cmom_proj133 + 0.5854102*cmom_proj136 + 0.5854102*cmom_proj138 + cmom_proj142 + cmom_proj144 - 0.5854102*cmom_proj152 + cmom_proj159 + cmom_proj372 + cmom_proj66 + 0.5854102*cmom_proj75 + cmom_proj79 - 0.5854102*cmom_proj86 - 0.5854102*cmom_proj88 - 0.5854102*cmom_proj93 - 0.5854102*cmom_proj98 - 0.40000000301872*dUdt_3_1;
            mom_proj[10]=cmom_proj193 + 0.5854102*cmom_proj198 + cmom_proj201 - 0.5854102*cmom_proj203 - 0.5854102*cmom_proj205 - 0.5854102*cmom_proj207 - 0.5854102*cmom_proj209 + cmom_proj212 + cmom_proj214 + cmom_proj216 + cmom_proj218 + 0.5854102*cmom_proj223 + 0.5854102*cmom_proj225 + cmom_proj228 + cmom_proj230 + 0.5854102*cmom_proj233 + cmom_proj236 + 0.5854102*cmom_proj238 + 0.5854102*cmom_proj240 + cmom_proj243 + cmom_proj245 - 0.5854102*cmom_proj248 + cmom_proj251 + cmom_proj373 - 0.40000000301872*dUdt_3_2;
            mom_proj[11]=cmom_proj279 + 0.5854102*cmom_proj284 + cmom_proj287 - 0.5854102*cmom_proj289 - 0.5854102*cmom_proj291 - 0.5854102*cmom_proj293 - 0.5854102*cmom_proj295 + cmom_proj298 + cmom_proj300 + cmom_proj302 + cmom_proj304 + 0.5854102*cmom_proj309 + 0.5854102*cmom_proj311 + cmom_proj314 + cmom_proj316 + 0.5854102*cmom_proj319 + cmom_proj322 + 0.5854102*cmom_proj324 + 0.5854102*cmom_proj326 + cmom_proj329 + cmom_proj331 - 0.5854102*cmom_proj334 + cmom_proj337 + cmom_proj374 - 0.40000000301872*dUdt_3_3;

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    mom_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * dim;
        auto& r_mom_proj = r_geometry[i_node].GetValue(MOMENTUM_PROJECTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom_proj[d] += mom_proj[aux + d];
        }
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &m_ext = data.m_ext;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_2_0 = data.dUdt(2, 0);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    // Calculate shock capturing values
    BoundedVector<double, 3> rho_proj;

    const double crho_proj0 =             -0.25*dUdt_1_0;
const double crho_proj1 =             0.25*m_ext[1];
const double crho_proj2 =             DN_DX_0_0*U_0_1;
const double crho_proj3 =             DN_DX_0_1*U_0_2;
const double crho_proj4 =             DN_DX_1_0*U_1_1;
const double crho_proj5 =             DN_DX_1_1*U_1_2;
const double crho_proj6 =             DN_DX_2_0*U_2_1;
const double crho_proj7 =             DN_DX_2_1*U_2_2;
const double crho_proj8 =             -1.0*crho_proj2 - 1.0*crho_proj3 - 1.0*crho_proj4 - 1.0*crho_proj5 - 1.0*crho_proj6 - 1.0*crho_proj7 - 0.25*dUdt_2_0 + 0.25*m_ext[2];
const double crho_proj9 =             -0.25*dUdt_0_0;
const double crho_proj10 =             0.25*m_ext[0];
            rho_proj[0]=crho_proj0 + crho_proj1 + crho_proj8 - 0.5*dUdt_0_0 + 0.5*m_ext[0];
            rho_proj[1]=crho_proj10 + crho_proj8 + crho_proj9 - 0.5*dUdt_1_0 + 0.5*m_ext[1];
            rho_proj[2]=crho_proj0 + crho_proj1 + crho_proj10 - 1.0*crho_proj2 - 1.0*crho_proj3 - 1.0*crho_proj4 - 1.0*crho_proj5 - 1.0*crho_proj6 - 1.0*crho_proj7 + crho_proj9 - 0.5*dUdt_2_0 + 0.5*m_ext[2];

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rho_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(DENSITY_PROJECTION) += rho_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::CalculateDensityProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &m_ext = data.m_ext;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);

    const double &dUdt_0_0 = data.dUdt(0, 0);
    const double &dUdt_1_0 = data.dUdt(1, 0);
    const double &dUdt_2_0 = data.dUdt(2, 0);
    const double &dUdt_3_0 = data.dUdt(3, 0);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    // Calculate shock capturing values
    BoundedVector<double, 4> rho_proj;

    const double crho_proj0 =             -0.19999999899376*dUdt_2_0;
const double crho_proj1 =             -0.19999999899376*dUdt_3_0;
const double crho_proj2 =             0.19999999899376*m_ext[2];
const double crho_proj3 =             0.19999999899376*m_ext[3];
const double crho_proj4 =             -1.0*DN_DX_0_0*U_0_1;
const double crho_proj5 =             -1.0*DN_DX_0_1*U_0_2;
const double crho_proj6 =             -1.0*DN_DX_0_2*U_0_3;
const double crho_proj7 =             -1.0*DN_DX_1_0*U_1_1;
const double crho_proj8 =             -1.0*DN_DX_1_1*U_1_2;
const double crho_proj9 =             -1.0*DN_DX_1_2*U_1_3;
const double crho_proj10 =             -1.0*DN_DX_2_0*U_2_1;
const double crho_proj11 =             -1.0*DN_DX_2_1*U_2_2;
const double crho_proj12 =             -1.0*DN_DX_2_2*U_2_3;
const double crho_proj13 =             -1.0*DN_DX_3_0*U_3_1;
const double crho_proj14 =             -1.0*DN_DX_3_1*U_3_2;
const double crho_proj15 =             -1.0*DN_DX_3_2*U_3_3;
const double crho_proj16 =             crho_proj0 + crho_proj1 + crho_proj10 + crho_proj11 + crho_proj12 + crho_proj13 + crho_proj14 + crho_proj15 + crho_proj2 + crho_proj3 + crho_proj4 + crho_proj5 + crho_proj6 + crho_proj7 + crho_proj8 + crho_proj9;
const double crho_proj17 =             crho_proj10 + crho_proj11 + crho_proj12 + crho_proj13 + crho_proj14 + crho_proj15 + crho_proj4 + crho_proj5 + crho_proj6 + crho_proj7 + crho_proj8 + crho_proj9 - 0.19999999899376*dUdt_0_0 - 0.19999999899376*dUdt_1_0 + 0.19999999899376*m_ext[0] + 0.19999999899376*m_ext[1];
            rho_proj[0]=crho_proj16 - 0.40000000301872*dUdt_0_0 - 0.19999999899376*dUdt_1_0 + 0.40000000301872*m_ext[0] + 0.19999999899376*m_ext[1];
            rho_proj[1]=crho_proj16 - 0.19999999899376*dUdt_0_0 - 0.40000000301872*dUdt_1_0 + 0.19999999899376*m_ext[0] + 0.40000000301872*m_ext[1];
            rho_proj[2]=crho_proj1 + crho_proj17 + crho_proj3 - 0.40000000301872*dUdt_2_0 + 0.40000000301872*m_ext[2];
            rho_proj[3]=crho_proj0 + crho_proj17 + crho_proj2 - 0.40000000301872*dUdt_3_0 + 0.40000000301872*m_ext[3];

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rho_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(DENSITY_PROJECTION) += rho_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 2> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);

    const double &dUdt_0_3 = data.dUdt(0, 3);
    const double &dUdt_1_3 = data.dUdt(1, 3);
    const double &dUdt_2_3 = data.dUdt(2, 3);

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    // Calculate shock capturing values
    BoundedVector<double, 3> tot_ener_proj;

    const double ctot_ener_proj0 =             -0.25*dUdt_1_3;
const double ctot_ener_proj1 =             0.166666666666667*U_0_0;
const double ctot_ener_proj2 =             0.166666666666667*U_2_0;
const double ctot_ener_proj3 =             0.666666666666667*U_1_0 + ctot_ener_proj1 + ctot_ener_proj2;
const double ctot_ener_proj4 =             0.166666666666667*r_ext[0];
const double ctot_ener_proj5 =             0.166666666666667*r_ext[2];
const double ctot_ener_proj6 =             ctot_ener_proj3*(ctot_ener_proj4 + ctot_ener_proj5 + 0.666666666666667*r_ext[1]);
const double ctot_ener_proj7 =             0.166666666666667*ctot_ener_proj6;
const double ctot_ener_proj8 =             0.166666666666667*U_0_1;
const double ctot_ener_proj9 =             0.166666666666667*U_2_1;
const double ctot_ener_proj10 =             0.666666666666667*U_1_1 + ctot_ener_proj8 + ctot_ener_proj9;
const double ctot_ener_proj11 =             0.166666666666667*f_ext(0,0);
const double ctot_ener_proj12 =             0.166666666666667*f_ext(2,0);
const double ctot_ener_proj13 =             ctot_ener_proj10*(ctot_ener_proj11 + ctot_ener_proj12 + 0.666666666666667*f_ext(1,0));
const double ctot_ener_proj14 =             0.166666666666667*ctot_ener_proj13;
const double ctot_ener_proj15 =             0.166666666666667*U_0_2;
const double ctot_ener_proj16 =             0.166666666666667*U_2_2;
const double ctot_ener_proj17 =             0.666666666666667*U_1_2 + ctot_ener_proj15 + ctot_ener_proj16;
const double ctot_ener_proj18 =             0.166666666666667*f_ext(0,1);
const double ctot_ener_proj19 =             0.166666666666667*f_ext(2,1);
const double ctot_ener_proj20 =             ctot_ener_proj17*(ctot_ener_proj18 + ctot_ener_proj19 + 0.666666666666667*f_ext(1,1));
const double ctot_ener_proj21 =             0.166666666666667*ctot_ener_proj20;
const double ctot_ener_proj22 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3;
const double ctot_ener_proj23 =             1.0/ctot_ener_proj3;
const double ctot_ener_proj24 =             ctot_ener_proj10*ctot_ener_proj22*ctot_ener_proj23*gamma;
const double ctot_ener_proj25 =             -0.166666666666667*ctot_ener_proj24;
const double ctot_ener_proj26 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3;
const double ctot_ener_proj27 =             ctot_ener_proj17*ctot_ener_proj23*ctot_ener_proj26*gamma;
const double ctot_ener_proj28 =             -0.166666666666667*ctot_ener_proj27;
const double ctot_ener_proj29 =             gamma - 1;
const double ctot_ener_proj30 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2;
const double ctot_ener_proj31 =             pow(ctot_ener_proj3, -2);
const double ctot_ener_proj32 =             ctot_ener_proj10*ctot_ener_proj17*ctot_ener_proj29*ctot_ener_proj30*ctot_ener_proj31;
const double ctot_ener_proj33 =             0.166666666666667*ctot_ener_proj32;
const double ctot_ener_proj34 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1;
const double ctot_ener_proj35 =             ctot_ener_proj10*ctot_ener_proj17*ctot_ener_proj29*ctot_ener_proj31*ctot_ener_proj34;
const double ctot_ener_proj36 =             0.166666666666667*ctot_ener_proj35;
const double ctot_ener_proj37 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1;
const double ctot_ener_proj38 =             pow(ctot_ener_proj10, 2);
const double ctot_ener_proj39 =             1.0*ctot_ener_proj23*ctot_ener_proj29;
const double ctot_ener_proj40 =             0.166666666666667*U_0_3;
const double ctot_ener_proj41 =             -ctot_ener_proj40;
const double ctot_ener_proj42 =             0.666666666666667*U_1_3;
const double ctot_ener_proj43 =             0.166666666666667*U_2_3;
const double ctot_ener_proj44 =             -ctot_ener_proj43;
const double ctot_ener_proj45 =             pow(ctot_ener_proj17, 2);
const double ctot_ener_proj46 =             -ctot_ener_proj29*(-ctot_ener_proj23*(0.5*ctot_ener_proj38 + 0.5*ctot_ener_proj45) + ctot_ener_proj40 + ctot_ener_proj42 + ctot_ener_proj43) + ctot_ener_proj41 - ctot_ener_proj42 + ctot_ener_proj44;
const double ctot_ener_proj47 =             ctot_ener_proj23*ctot_ener_proj37*(ctot_ener_proj38*ctot_ener_proj39 + ctot_ener_proj46);
const double ctot_ener_proj48 =             0.166666666666667*ctot_ener_proj47;
const double ctot_ener_proj49 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2;
const double ctot_ener_proj50 =             ctot_ener_proj23*ctot_ener_proj49*(ctot_ener_proj39*ctot_ener_proj45 + ctot_ener_proj46);
const double ctot_ener_proj51 =             0.166666666666667*ctot_ener_proj50;
const double ctot_ener_proj52 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0;
const double ctot_ener_proj53 =             0.5*gamma - 0.5;
const double ctot_ener_proj54 =             ctot_ener_proj23*ctot_ener_proj53*(ctot_ener_proj38 + ctot_ener_proj45) + ctot_ener_proj46;
const double ctot_ener_proj55 =             ctot_ener_proj10*ctot_ener_proj31*ctot_ener_proj52*ctot_ener_proj54;
const double ctot_ener_proj56 =             -0.166666666666667*ctot_ener_proj55;
const double ctot_ener_proj57 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0;
const double ctot_ener_proj58 =             ctot_ener_proj17*ctot_ener_proj31*ctot_ener_proj54*ctot_ener_proj57;
const double ctot_ener_proj59 =             -0.166666666666667*ctot_ener_proj58;
const double ctot_ener_proj60 =             -0.25*dUdt_2_3;
const double ctot_ener_proj61 =             0.166666666666667*U_1_0;
const double ctot_ener_proj62 =             0.666666666666667*U_2_0 + ctot_ener_proj1 + ctot_ener_proj61;
const double ctot_ener_proj63 =             0.166666666666667*r_ext[1];
const double ctot_ener_proj64 =             ctot_ener_proj62*(ctot_ener_proj4 + ctot_ener_proj63 + 0.666666666666667*r_ext[2]);
const double ctot_ener_proj65 =             0.166666666666667*ctot_ener_proj64;
const double ctot_ener_proj66 =             0.666666666666667*U_0_0 + ctot_ener_proj2 + ctot_ener_proj61;
const double ctot_ener_proj67 =             ctot_ener_proj66*(ctot_ener_proj5 + ctot_ener_proj63 + 0.666666666666667*r_ext[0]);
const double ctot_ener_proj68 =             0.166666666666667*U_1_1;
const double ctot_ener_proj69 =             0.666666666666667*U_2_1 + ctot_ener_proj68 + ctot_ener_proj8;
const double ctot_ener_proj70 =             0.166666666666667*f_ext(1,0);
const double ctot_ener_proj71 =             ctot_ener_proj69*(ctot_ener_proj11 + ctot_ener_proj70 + 0.666666666666667*f_ext(2,0));
const double ctot_ener_proj72 =             0.166666666666667*ctot_ener_proj71;
const double ctot_ener_proj73 =             0.666666666666667*U_0_1 + ctot_ener_proj68 + ctot_ener_proj9;
const double ctot_ener_proj74 =             ctot_ener_proj73*(ctot_ener_proj12 + ctot_ener_proj70 + 0.666666666666667*f_ext(0,0));
const double ctot_ener_proj75 =             0.166666666666667*U_1_2;
const double ctot_ener_proj76 =             0.666666666666667*U_2_2 + ctot_ener_proj15 + ctot_ener_proj75;
const double ctot_ener_proj77 =             0.166666666666667*f_ext(1,1);
const double ctot_ener_proj78 =             ctot_ener_proj76*(ctot_ener_proj18 + ctot_ener_proj77 + 0.666666666666667*f_ext(2,1));
const double ctot_ener_proj79 =             0.166666666666667*ctot_ener_proj78;
const double ctot_ener_proj80 =             0.666666666666667*U_0_2 + ctot_ener_proj16 + ctot_ener_proj75;
const double ctot_ener_proj81 =             ctot_ener_proj80*(ctot_ener_proj19 + ctot_ener_proj77 + 0.666666666666667*f_ext(0,1));
const double ctot_ener_proj82 =             1.0/ctot_ener_proj62;
const double ctot_ener_proj83 =             ctot_ener_proj22*ctot_ener_proj69*ctot_ener_proj82*gamma;
const double ctot_ener_proj84 =             -0.166666666666667*ctot_ener_proj83;
const double ctot_ener_proj85 =             ctot_ener_proj26*ctot_ener_proj76*ctot_ener_proj82*gamma;
const double ctot_ener_proj86 =             -0.166666666666667*ctot_ener_proj85;
const double ctot_ener_proj87 =             1.0/ctot_ener_proj66;
const double ctot_ener_proj88 =             ctot_ener_proj22*ctot_ener_proj73*ctot_ener_proj87*gamma;
const double ctot_ener_proj89 =             ctot_ener_proj26*ctot_ener_proj80*ctot_ener_proj87*gamma;
const double ctot_ener_proj90 =             pow(ctot_ener_proj62, -2);
const double ctot_ener_proj91 =             ctot_ener_proj29*ctot_ener_proj30*ctot_ener_proj69*ctot_ener_proj76*ctot_ener_proj90;
const double ctot_ener_proj92 =             0.166666666666667*ctot_ener_proj91;
const double ctot_ener_proj93 =             ctot_ener_proj29*ctot_ener_proj34*ctot_ener_proj69*ctot_ener_proj76*ctot_ener_proj90;
const double ctot_ener_proj94 =             0.166666666666667*ctot_ener_proj93;
const double ctot_ener_proj95 =             pow(ctot_ener_proj66, -2);
const double ctot_ener_proj96 =             ctot_ener_proj29*ctot_ener_proj30*ctot_ener_proj73*ctot_ener_proj80*ctot_ener_proj95;
const double ctot_ener_proj97 =             ctot_ener_proj29*ctot_ener_proj34*ctot_ener_proj73*ctot_ener_proj80*ctot_ener_proj95;
const double ctot_ener_proj98 =             pow(ctot_ener_proj69, 2);
const double ctot_ener_proj99 =             1.0*ctot_ener_proj29*ctot_ener_proj82;
const double ctot_ener_proj100 =             0.166666666666667*U_1_3;
const double ctot_ener_proj101 =             -ctot_ener_proj100;
const double ctot_ener_proj102 =             0.666666666666667*U_2_3;
const double ctot_ener_proj103 =             pow(ctot_ener_proj76, 2);
const double ctot_ener_proj104 =             ctot_ener_proj101 - ctot_ener_proj102 - ctot_ener_proj29*(ctot_ener_proj100 + ctot_ener_proj102 + ctot_ener_proj40 - ctot_ener_proj82*(0.5*ctot_ener_proj103 + 0.5*ctot_ener_proj98)) + ctot_ener_proj41;
const double ctot_ener_proj105 =             ctot_ener_proj37*ctot_ener_proj82*(ctot_ener_proj104 + ctot_ener_proj98*ctot_ener_proj99);
const double ctot_ener_proj106 =             0.166666666666667*ctot_ener_proj105;
const double ctot_ener_proj107 =             ctot_ener_proj49*ctot_ener_proj82*(ctot_ener_proj103*ctot_ener_proj99 + ctot_ener_proj104);
const double ctot_ener_proj108 =             0.166666666666667*ctot_ener_proj107;
const double ctot_ener_proj109 =             pow(ctot_ener_proj73, 2);
const double ctot_ener_proj110 =             1.0*ctot_ener_proj29*ctot_ener_proj87;
const double ctot_ener_proj111 =             0.666666666666667*U_0_3;
const double ctot_ener_proj112 =             pow(ctot_ener_proj80, 2);
const double ctot_ener_proj113 =             ctot_ener_proj101 - ctot_ener_proj111 - ctot_ener_proj29*(ctot_ener_proj100 + ctot_ener_proj111 + ctot_ener_proj43 - ctot_ener_proj87*(0.5*ctot_ener_proj109 + 0.5*ctot_ener_proj112)) + ctot_ener_proj44;
const double ctot_ener_proj114 =             ctot_ener_proj37*ctot_ener_proj87*(ctot_ener_proj109*ctot_ener_proj110 + ctot_ener_proj113);
const double ctot_ener_proj115 =             ctot_ener_proj49*ctot_ener_proj87*(ctot_ener_proj110*ctot_ener_proj112 + ctot_ener_proj113);
const double ctot_ener_proj116 =             ctot_ener_proj104 + ctot_ener_proj53*ctot_ener_proj82*(ctot_ener_proj103 + ctot_ener_proj98);
const double ctot_ener_proj117 =             ctot_ener_proj116*ctot_ener_proj52*ctot_ener_proj69*ctot_ener_proj90;
const double ctot_ener_proj118 =             -0.166666666666667*ctot_ener_proj117;
const double ctot_ener_proj119 =             ctot_ener_proj116*ctot_ener_proj57*ctot_ener_proj76*ctot_ener_proj90;
const double ctot_ener_proj120 =             -0.166666666666667*ctot_ener_proj119;
const double ctot_ener_proj121 =             ctot_ener_proj113 + ctot_ener_proj53*ctot_ener_proj87*(ctot_ener_proj109 + ctot_ener_proj112);
const double ctot_ener_proj122 =             ctot_ener_proj121*ctot_ener_proj52*ctot_ener_proj73*ctot_ener_proj95;
const double ctot_ener_proj123 =             ctot_ener_proj121*ctot_ener_proj57*ctot_ener_proj80*ctot_ener_proj95;
const double ctot_ener_proj124 =             0.166666666666667*ctot_ener_proj114 + 0.166666666666667*ctot_ener_proj115 - 0.166666666666667*ctot_ener_proj122 - 0.166666666666667*ctot_ener_proj123 + 0.166666666666667*ctot_ener_proj67 + 0.166666666666667*ctot_ener_proj74 + 0.166666666666667*ctot_ener_proj81 - 0.166666666666667*ctot_ener_proj88 - 0.166666666666667*ctot_ener_proj89 + 0.166666666666667*ctot_ener_proj96 + 0.166666666666667*ctot_ener_proj97 - 0.25*dUdt_0_3;
            tot_ener_proj[0]=ctot_ener_proj0 + ctot_ener_proj106 + ctot_ener_proj108 + 0.666666666666667*ctot_ener_proj114 + 0.666666666666667*ctot_ener_proj115 + ctot_ener_proj118 + ctot_ener_proj120 - 0.666666666666667*ctot_ener_proj122 - 0.666666666666667*ctot_ener_proj123 + ctot_ener_proj14 + ctot_ener_proj21 + ctot_ener_proj25 + ctot_ener_proj28 + ctot_ener_proj33 + ctot_ener_proj36 + ctot_ener_proj48 + ctot_ener_proj51 + ctot_ener_proj56 + ctot_ener_proj59 + ctot_ener_proj60 + ctot_ener_proj65 + 0.666666666666667*ctot_ener_proj67 + ctot_ener_proj7 + ctot_ener_proj72 + 0.666666666666667*ctot_ener_proj74 + ctot_ener_proj79 + 0.666666666666667*ctot_ener_proj81 + ctot_ener_proj84 + ctot_ener_proj86 - 0.666666666666667*ctot_ener_proj88 - 0.666666666666667*ctot_ener_proj89 + ctot_ener_proj92 + ctot_ener_proj94 + 0.666666666666667*ctot_ener_proj96 + 0.666666666666667*ctot_ener_proj97 - 0.5*dUdt_0_3;
            tot_ener_proj[1]=ctot_ener_proj106 + ctot_ener_proj108 + ctot_ener_proj118 + ctot_ener_proj120 + ctot_ener_proj124 + 0.666666666666667*ctot_ener_proj13 + 0.666666666666667*ctot_ener_proj20 - 0.666666666666667*ctot_ener_proj24 - 0.666666666666667*ctot_ener_proj27 + 0.666666666666667*ctot_ener_proj32 + 0.666666666666667*ctot_ener_proj35 + 0.666666666666667*ctot_ener_proj47 + 0.666666666666667*ctot_ener_proj50 - 0.666666666666667*ctot_ener_proj55 - 0.666666666666667*ctot_ener_proj58 + 0.666666666666667*ctot_ener_proj6 + ctot_ener_proj60 + ctot_ener_proj65 + ctot_ener_proj72 + ctot_ener_proj79 + ctot_ener_proj84 + ctot_ener_proj86 + ctot_ener_proj92 + ctot_ener_proj94 - 0.5*dUdt_1_3;
            tot_ener_proj[2]=ctot_ener_proj0 + 0.666666666666667*ctot_ener_proj105 + 0.666666666666667*ctot_ener_proj107 - 0.666666666666667*ctot_ener_proj117 - 0.666666666666667*ctot_ener_proj119 + ctot_ener_proj124 + ctot_ener_proj14 + ctot_ener_proj21 + ctot_ener_proj25 + ctot_ener_proj28 + ctot_ener_proj33 + ctot_ener_proj36 + ctot_ener_proj48 + ctot_ener_proj51 + ctot_ener_proj56 + ctot_ener_proj59 + 0.666666666666667*ctot_ener_proj64 + ctot_ener_proj7 + 0.666666666666667*ctot_ener_proj71 + 0.666666666666667*ctot_ener_proj78 - 0.666666666666667*ctot_ener_proj83 - 0.666666666666667*ctot_ener_proj85 + 0.666666666666667*ctot_ener_proj91 + 0.666666666666667*ctot_ener_proj93 - 0.5*dUdt_2_3;

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    tot_ener_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION) += tot_ener_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::CalculateTotalEnergyProjection(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const array_1d<double, n_nodes> &r_ext = data.r_ext;
    const BoundedMatrix<double, n_nodes, 3> &f_ext = data.f_ext;
    const double gamma = data.gamma;

    // Solution vector values and time derivatives from nodal data
    // This is intentionally done in this way to limit the matrix acceses
    // The notation U_i_j DOF j value in node i
    const double &U_0_0 = data.U(0, 0);
    const double &U_0_1 = data.U(0, 1);
    const double &U_0_2 = data.U(0, 2);
    const double &U_0_3 = data.U(0, 3);
    const double &U_0_4 = data.U(0, 4);
    const double &U_1_0 = data.U(1, 0);
    const double &U_1_1 = data.U(1, 1);
    const double &U_1_2 = data.U(1, 2);
    const double &U_1_3 = data.U(1, 3);
    const double &U_1_4 = data.U(1, 4);
    const double &U_2_0 = data.U(2, 0);
    const double &U_2_1 = data.U(2, 1);
    const double &U_2_2 = data.U(2, 2);
    const double &U_2_3 = data.U(2, 3);
    const double &U_2_4 = data.U(2, 4);
    const double &U_3_0 = data.U(3, 0);
    const double &U_3_1 = data.U(3, 1);
    const double &U_3_2 = data.U(3, 2);
    const double &U_3_3 = data.U(3, 3);
    const double &U_3_4 = data.U(3, 4);

    const double &dUdt_0_4 = data.dUdt(0, 4);
    const double &dUdt_1_4 = data.dUdt(1, 4);
    const double &dUdt_2_4 = data.dUdt(2, 4);
    const double &dUdt_3_4 = data.dUdt(3, 4);

    // Hardcoded shape functions gradients for linear tetrahedra element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    // Calculate shock capturing values
    BoundedVector<double, 4> tot_ener_proj;

    const double ctot_ener_proj0 =             0.1381966*U_0_0;
const double ctot_ener_proj1 =             0.1381966*U_2_0;
const double ctot_ener_proj2 =             0.1381966*U_3_0;
const double ctot_ener_proj3 =             ctot_ener_proj1 + ctot_ener_proj2;
const double ctot_ener_proj4 =             0.5854102*U_1_0 + ctot_ener_proj0 + ctot_ener_proj3;
const double ctot_ener_proj5 =             0.1381966*r_ext[0];
const double ctot_ener_proj6 =             0.1381966*r_ext[2];
const double ctot_ener_proj7 =             0.1381966*r_ext[3];
const double ctot_ener_proj8 =             ctot_ener_proj6 + ctot_ener_proj7;
const double ctot_ener_proj9 =             ctot_ener_proj4*(ctot_ener_proj5 + ctot_ener_proj8 + 0.5854102*r_ext[1]);
const double ctot_ener_proj10 =             0.1381966*ctot_ener_proj9;
const double ctot_ener_proj11 =             0.1381966*U_0_1;
const double ctot_ener_proj12 =             0.1381966*U_2_1;
const double ctot_ener_proj13 =             0.1381966*U_3_1;
const double ctot_ener_proj14 =             ctot_ener_proj12 + ctot_ener_proj13;
const double ctot_ener_proj15 =             0.5854102*U_1_1 + ctot_ener_proj11 + ctot_ener_proj14;
const double ctot_ener_proj16 =             0.1381966*f_ext(0,0);
const double ctot_ener_proj17 =             0.1381966*f_ext(2,0);
const double ctot_ener_proj18 =             0.1381966*f_ext(3,0);
const double ctot_ener_proj19 =             ctot_ener_proj17 + ctot_ener_proj18;
const double ctot_ener_proj20 =             ctot_ener_proj15*(ctot_ener_proj16 + ctot_ener_proj19 + 0.5854102*f_ext(1,0));
const double ctot_ener_proj21 =             0.1381966*ctot_ener_proj20;
const double ctot_ener_proj22 =             0.1381966*U_0_2;
const double ctot_ener_proj23 =             0.1381966*U_2_2;
const double ctot_ener_proj24 =             0.1381966*U_3_2;
const double ctot_ener_proj25 =             ctot_ener_proj23 + ctot_ener_proj24;
const double ctot_ener_proj26 =             0.5854102*U_1_2 + ctot_ener_proj22 + ctot_ener_proj25;
const double ctot_ener_proj27 =             0.1381966*f_ext(0,1);
const double ctot_ener_proj28 =             0.1381966*f_ext(2,1);
const double ctot_ener_proj29 =             0.1381966*f_ext(3,1);
const double ctot_ener_proj30 =             ctot_ener_proj28 + ctot_ener_proj29;
const double ctot_ener_proj31 =             ctot_ener_proj26*(ctot_ener_proj27 + ctot_ener_proj30 + 0.5854102*f_ext(1,1));
const double ctot_ener_proj32 =             0.1381966*ctot_ener_proj31;
const double ctot_ener_proj33 =             0.1381966*U_0_3;
const double ctot_ener_proj34 =             0.1381966*U_2_3;
const double ctot_ener_proj35 =             0.1381966*U_3_3;
const double ctot_ener_proj36 =             ctot_ener_proj34 + ctot_ener_proj35;
const double ctot_ener_proj37 =             0.5854102*U_1_3 + ctot_ener_proj33 + ctot_ener_proj36;
const double ctot_ener_proj38 =             0.1381966*f_ext(0,2);
const double ctot_ener_proj39 =             0.1381966*f_ext(2,2);
const double ctot_ener_proj40 =             0.1381966*f_ext(3,2);
const double ctot_ener_proj41 =             ctot_ener_proj39 + ctot_ener_proj40;
const double ctot_ener_proj42 =             ctot_ener_proj37*(ctot_ener_proj38 + ctot_ener_proj41 + 0.5854102*f_ext(1,2));
const double ctot_ener_proj43 =             0.1381966*ctot_ener_proj42;
const double ctot_ener_proj44 =             DN_DX_0_0*U_0_4 + DN_DX_1_0*U_1_4 + DN_DX_2_0*U_2_4 + DN_DX_3_0*U_3_4;
const double ctot_ener_proj45 =             1.0/ctot_ener_proj4;
const double ctot_ener_proj46 =             ctot_ener_proj15*ctot_ener_proj44*ctot_ener_proj45*gamma;
const double ctot_ener_proj47 =             -0.1381966*ctot_ener_proj46;
const double ctot_ener_proj48 =             DN_DX_0_1*U_0_4 + DN_DX_1_1*U_1_4 + DN_DX_2_1*U_2_4 + DN_DX_3_1*U_3_4;
const double ctot_ener_proj49 =             ctot_ener_proj26*ctot_ener_proj45*ctot_ener_proj48*gamma;
const double ctot_ener_proj50 =             -0.1381966*ctot_ener_proj49;
const double ctot_ener_proj51 =             DN_DX_0_2*U_0_4 + DN_DX_1_2*U_1_4 + DN_DX_2_2*U_2_4 + DN_DX_3_2*U_3_4;
const double ctot_ener_proj52 =             ctot_ener_proj37*ctot_ener_proj45*ctot_ener_proj51*gamma;
const double ctot_ener_proj53 =             -0.1381966*ctot_ener_proj52;
const double ctot_ener_proj54 =             gamma - 1;
const double ctot_ener_proj55 =             DN_DX_0_0*U_0_2 + DN_DX_1_0*U_1_2 + DN_DX_2_0*U_2_2 + DN_DX_3_0*U_3_2;
const double ctot_ener_proj56 =             pow(ctot_ener_proj4, -2);
const double ctot_ener_proj57 =             ctot_ener_proj15*ctot_ener_proj26*ctot_ener_proj54*ctot_ener_proj55*ctot_ener_proj56;
const double ctot_ener_proj58 =             0.1381966*ctot_ener_proj57;
const double ctot_ener_proj59 =             DN_DX_0_1*U_0_1 + DN_DX_1_1*U_1_1 + DN_DX_2_1*U_2_1 + DN_DX_3_1*U_3_1;
const double ctot_ener_proj60 =             ctot_ener_proj15*ctot_ener_proj26*ctot_ener_proj54*ctot_ener_proj56*ctot_ener_proj59;
const double ctot_ener_proj61 =             0.1381966*ctot_ener_proj60;
const double ctot_ener_proj62 =             DN_DX_0_0*U_0_3 + DN_DX_1_0*U_1_3 + DN_DX_2_0*U_2_3 + DN_DX_3_0*U_3_3;
const double ctot_ener_proj63 =             ctot_ener_proj15*ctot_ener_proj37*ctot_ener_proj54*ctot_ener_proj56*ctot_ener_proj62;
const double ctot_ener_proj64 =             0.1381966*ctot_ener_proj63;
const double ctot_ener_proj65 =             DN_DX_0_2*U_0_1 + DN_DX_1_2*U_1_1 + DN_DX_2_2*U_2_1 + DN_DX_3_2*U_3_1;
const double ctot_ener_proj66 =             ctot_ener_proj15*ctot_ener_proj37*ctot_ener_proj54*ctot_ener_proj56*ctot_ener_proj65;
const double ctot_ener_proj67 =             0.1381966*ctot_ener_proj66;
const double ctot_ener_proj68 =             DN_DX_0_1*U_0_3 + DN_DX_1_1*U_1_3 + DN_DX_2_1*U_2_3 + DN_DX_3_1*U_3_3;
const double ctot_ener_proj69 =             ctot_ener_proj26*ctot_ener_proj37*ctot_ener_proj54*ctot_ener_proj56*ctot_ener_proj68;
const double ctot_ener_proj70 =             0.1381966*ctot_ener_proj69;
const double ctot_ener_proj71 =             DN_DX_0_2*U_0_2 + DN_DX_1_2*U_1_2 + DN_DX_2_2*U_2_2 + DN_DX_3_2*U_3_2;
const double ctot_ener_proj72 =             ctot_ener_proj26*ctot_ener_proj37*ctot_ener_proj54*ctot_ener_proj56*ctot_ener_proj71;
const double ctot_ener_proj73 =             0.1381966*ctot_ener_proj72;
const double ctot_ener_proj74 =             DN_DX_0_0*U_0_1 + DN_DX_1_0*U_1_1 + DN_DX_2_0*U_2_1 + DN_DX_3_0*U_3_1;
const double ctot_ener_proj75 =             pow(ctot_ener_proj15, 2);
const double ctot_ener_proj76 =             1.0*ctot_ener_proj45*ctot_ener_proj54;
const double ctot_ener_proj77 =             0.1381966*U_0_4;
const double ctot_ener_proj78 =             -ctot_ener_proj77;
const double ctot_ener_proj79 =             0.5854102*U_1_4;
const double ctot_ener_proj80 =             0.1381966*U_2_4;
const double ctot_ener_proj81 =             -ctot_ener_proj80;
const double ctot_ener_proj82 =             0.1381966*U_3_4;
const double ctot_ener_proj83 =             -ctot_ener_proj82;
const double ctot_ener_proj84 =             ctot_ener_proj80 + ctot_ener_proj82;
const double ctot_ener_proj85 =             pow(ctot_ener_proj26, 2);
const double ctot_ener_proj86 =             pow(ctot_ener_proj37, 2);
const double ctot_ener_proj87 =             -ctot_ener_proj54*(-ctot_ener_proj45*(0.5*ctot_ener_proj75 + 0.5*ctot_ener_proj85 + 0.5*ctot_ener_proj86) + ctot_ener_proj77 + ctot_ener_proj79 + ctot_ener_proj84) + ctot_ener_proj78 - ctot_ener_proj79 + ctot_ener_proj81 + ctot_ener_proj83;
const double ctot_ener_proj88 =             ctot_ener_proj45*ctot_ener_proj74*(ctot_ener_proj75*ctot_ener_proj76 + ctot_ener_proj87);
const double ctot_ener_proj89 =             0.1381966*ctot_ener_proj88;
const double ctot_ener_proj90 =             DN_DX_0_1*U_0_2 + DN_DX_1_1*U_1_2 + DN_DX_2_1*U_2_2 + DN_DX_3_1*U_3_2;
const double ctot_ener_proj91 =             ctot_ener_proj45*ctot_ener_proj90*(ctot_ener_proj76*ctot_ener_proj85 + ctot_ener_proj87);
const double ctot_ener_proj92 =             0.1381966*ctot_ener_proj91;
const double ctot_ener_proj93 =             DN_DX_0_2*U_0_3 + DN_DX_1_2*U_1_3 + DN_DX_2_2*U_2_3 + DN_DX_3_2*U_3_3;
const double ctot_ener_proj94 =             ctot_ener_proj45*ctot_ener_proj93*(ctot_ener_proj76*ctot_ener_proj86 + ctot_ener_proj87);
const double ctot_ener_proj95 =             0.1381966*ctot_ener_proj94;
const double ctot_ener_proj96 =             DN_DX_0_0*U_0_0 + DN_DX_1_0*U_1_0 + DN_DX_2_0*U_2_0 + DN_DX_3_0*U_3_0;
const double ctot_ener_proj97 =             0.5*gamma - 0.5;
const double ctot_ener_proj98 =             ctot_ener_proj45*ctot_ener_proj97*(ctot_ener_proj75 + ctot_ener_proj85 + ctot_ener_proj86) + ctot_ener_proj87;
const double ctot_ener_proj99 =             ctot_ener_proj15*ctot_ener_proj56*ctot_ener_proj96*ctot_ener_proj98;
const double ctot_ener_proj100 =             -0.1381966*ctot_ener_proj99;
const double ctot_ener_proj101 =             DN_DX_0_1*U_0_0 + DN_DX_1_1*U_1_0 + DN_DX_2_1*U_2_0 + DN_DX_3_1*U_3_0;
const double ctot_ener_proj102 =             ctot_ener_proj101*ctot_ener_proj26*ctot_ener_proj56*ctot_ener_proj98;
const double ctot_ener_proj103 =             -0.1381966*ctot_ener_proj102;
const double ctot_ener_proj104 =             DN_DX_0_2*U_0_0 + DN_DX_1_2*U_1_0 + DN_DX_2_2*U_2_0 + DN_DX_3_2*U_3_0;
const double ctot_ener_proj105 =             ctot_ener_proj104*ctot_ener_proj37*ctot_ener_proj56*ctot_ener_proj98;
const double ctot_ener_proj106 =             -0.1381966*ctot_ener_proj105;
const double ctot_ener_proj107 =             -0.19999999899376*dUdt_2_4;
const double ctot_ener_proj108 =             -0.19999999899376*dUdt_3_4;
const double ctot_ener_proj109 =             0.1381966*U_1_0;
const double ctot_ener_proj110 =             ctot_ener_proj0 + ctot_ener_proj109;
const double ctot_ener_proj111 =             0.5854102*U_3_0 + ctot_ener_proj1 + ctot_ener_proj110;
const double ctot_ener_proj112 =             0.1381966*r_ext[1];
const double ctot_ener_proj113 =             ctot_ener_proj112 + ctot_ener_proj5;
const double ctot_ener_proj114 =             ctot_ener_proj111*(ctot_ener_proj113 + ctot_ener_proj6 + 0.5854102*r_ext[3]);
const double ctot_ener_proj115 =             0.1381966*ctot_ener_proj114;
const double ctot_ener_proj116 =             0.5854102*U_2_0 + ctot_ener_proj110 + ctot_ener_proj2;
const double ctot_ener_proj117 =             ctot_ener_proj116*(ctot_ener_proj113 + ctot_ener_proj7 + 0.5854102*r_ext[2]);
const double ctot_ener_proj118 =             0.1381966*ctot_ener_proj117;
const double ctot_ener_proj119 =             0.5854102*U_0_0 + ctot_ener_proj109 + ctot_ener_proj3;
const double ctot_ener_proj120 =             ctot_ener_proj119*(ctot_ener_proj112 + ctot_ener_proj8 + 0.5854102*r_ext[0]);
const double ctot_ener_proj121 =             0.1381966*U_1_1;
const double ctot_ener_proj122 =             ctot_ener_proj11 + ctot_ener_proj121;
const double ctot_ener_proj123 =             0.5854102*U_3_1 + ctot_ener_proj12 + ctot_ener_proj122;
const double ctot_ener_proj124 =             0.1381966*f_ext(1,0);
const double ctot_ener_proj125 =             ctot_ener_proj124 + ctot_ener_proj16;
const double ctot_ener_proj126 =             ctot_ener_proj123*(ctot_ener_proj125 + ctot_ener_proj17 + 0.5854102*f_ext(3,0));
const double ctot_ener_proj127 =             0.1381966*ctot_ener_proj126;
const double ctot_ener_proj128 =             0.5854102*U_2_1 + ctot_ener_proj122 + ctot_ener_proj13;
const double ctot_ener_proj129 =             ctot_ener_proj128*(ctot_ener_proj125 + ctot_ener_proj18 + 0.5854102*f_ext(2,0));
const double ctot_ener_proj130 =             0.1381966*ctot_ener_proj129;
const double ctot_ener_proj131 =             0.5854102*U_0_1 + ctot_ener_proj121 + ctot_ener_proj14;
const double ctot_ener_proj132 =             ctot_ener_proj131*(ctot_ener_proj124 + ctot_ener_proj19 + 0.5854102*f_ext(0,0));
const double ctot_ener_proj133 =             0.1381966*U_1_2;
const double ctot_ener_proj134 =             ctot_ener_proj133 + ctot_ener_proj22;
const double ctot_ener_proj135 =             0.5854102*U_3_2 + ctot_ener_proj134 + ctot_ener_proj23;
const double ctot_ener_proj136 =             0.1381966*f_ext(1,1);
const double ctot_ener_proj137 =             ctot_ener_proj136 + ctot_ener_proj27;
const double ctot_ener_proj138 =             ctot_ener_proj135*(ctot_ener_proj137 + ctot_ener_proj28 + 0.5854102*f_ext(3,1));
const double ctot_ener_proj139 =             0.1381966*ctot_ener_proj138;
const double ctot_ener_proj140 =             0.5854102*U_2_2 + ctot_ener_proj134 + ctot_ener_proj24;
const double ctot_ener_proj141 =             ctot_ener_proj140*(ctot_ener_proj137 + ctot_ener_proj29 + 0.5854102*f_ext(2,1));
const double ctot_ener_proj142 =             0.1381966*ctot_ener_proj141;
const double ctot_ener_proj143 =             0.5854102*U_0_2 + ctot_ener_proj133 + ctot_ener_proj25;
const double ctot_ener_proj144 =             ctot_ener_proj143*(ctot_ener_proj136 + ctot_ener_proj30 + 0.5854102*f_ext(0,1));
const double ctot_ener_proj145 =             0.1381966*U_1_3;
const double ctot_ener_proj146 =             ctot_ener_proj145 + ctot_ener_proj33;
const double ctot_ener_proj147 =             0.5854102*U_3_3 + ctot_ener_proj146 + ctot_ener_proj34;
const double ctot_ener_proj148 =             0.1381966*f_ext(1,2);
const double ctot_ener_proj149 =             ctot_ener_proj148 + ctot_ener_proj38;
const double ctot_ener_proj150 =             ctot_ener_proj147*(ctot_ener_proj149 + ctot_ener_proj39 + 0.5854102*f_ext(3,2));
const double ctot_ener_proj151 =             0.1381966*ctot_ener_proj150;
const double ctot_ener_proj152 =             0.5854102*U_2_3 + ctot_ener_proj146 + ctot_ener_proj35;
const double ctot_ener_proj153 =             ctot_ener_proj152*(ctot_ener_proj149 + ctot_ener_proj40 + 0.5854102*f_ext(2,2));
const double ctot_ener_proj154 =             0.1381966*ctot_ener_proj153;
const double ctot_ener_proj155 =             0.5854102*U_0_3 + ctot_ener_proj145 + ctot_ener_proj36;
const double ctot_ener_proj156 =             ctot_ener_proj155*(ctot_ener_proj148 + ctot_ener_proj41 + 0.5854102*f_ext(0,2));
const double ctot_ener_proj157 =             1.0/ctot_ener_proj111;
const double ctot_ener_proj158 =             ctot_ener_proj123*ctot_ener_proj157*ctot_ener_proj44*gamma;
const double ctot_ener_proj159 =             -0.1381966*ctot_ener_proj158;
const double ctot_ener_proj160 =             ctot_ener_proj135*ctot_ener_proj157*ctot_ener_proj48*gamma;
const double ctot_ener_proj161 =             -0.1381966*ctot_ener_proj160;
const double ctot_ener_proj162 =             ctot_ener_proj147*ctot_ener_proj157*ctot_ener_proj51*gamma;
const double ctot_ener_proj163 =             -0.1381966*ctot_ener_proj162;
const double ctot_ener_proj164 =             1.0/ctot_ener_proj116;
const double ctot_ener_proj165 =             ctot_ener_proj128*ctot_ener_proj164*ctot_ener_proj44*gamma;
const double ctot_ener_proj166 =             -0.1381966*ctot_ener_proj165;
const double ctot_ener_proj167 =             ctot_ener_proj140*ctot_ener_proj164*ctot_ener_proj48*gamma;
const double ctot_ener_proj168 =             -0.1381966*ctot_ener_proj167;
const double ctot_ener_proj169 =             ctot_ener_proj152*ctot_ener_proj164*ctot_ener_proj51*gamma;
const double ctot_ener_proj170 =             -0.1381966*ctot_ener_proj169;
const double ctot_ener_proj171 =             1.0/ctot_ener_proj119;
const double ctot_ener_proj172 =             ctot_ener_proj131*ctot_ener_proj171*ctot_ener_proj44*gamma;
const double ctot_ener_proj173 =             ctot_ener_proj143*ctot_ener_proj171*ctot_ener_proj48*gamma;
const double ctot_ener_proj174 =             ctot_ener_proj155*ctot_ener_proj171*ctot_ener_proj51*gamma;
const double ctot_ener_proj175 =             pow(ctot_ener_proj111, -2);
const double ctot_ener_proj176 =             ctot_ener_proj123*ctot_ener_proj135*ctot_ener_proj175*ctot_ener_proj54*ctot_ener_proj55;
const double ctot_ener_proj177 =             0.1381966*ctot_ener_proj176;
const double ctot_ener_proj178 =             ctot_ener_proj123*ctot_ener_proj135*ctot_ener_proj175*ctot_ener_proj54*ctot_ener_proj59;
const double ctot_ener_proj179 =             0.1381966*ctot_ener_proj178;
const double ctot_ener_proj180 =             ctot_ener_proj123*ctot_ener_proj147*ctot_ener_proj175*ctot_ener_proj54*ctot_ener_proj62;
const double ctot_ener_proj181 =             0.1381966*ctot_ener_proj180;
const double ctot_ener_proj182 =             ctot_ener_proj123*ctot_ener_proj147*ctot_ener_proj175*ctot_ener_proj54*ctot_ener_proj65;
const double ctot_ener_proj183 =             0.1381966*ctot_ener_proj182;
const double ctot_ener_proj184 =             ctot_ener_proj135*ctot_ener_proj147*ctot_ener_proj175*ctot_ener_proj54*ctot_ener_proj68;
const double ctot_ener_proj185 =             0.1381966*ctot_ener_proj184;
const double ctot_ener_proj186 =             ctot_ener_proj135*ctot_ener_proj147*ctot_ener_proj175*ctot_ener_proj54*ctot_ener_proj71;
const double ctot_ener_proj187 =             0.1381966*ctot_ener_proj186;
const double ctot_ener_proj188 =             pow(ctot_ener_proj116, -2);
const double ctot_ener_proj189 =             ctot_ener_proj128*ctot_ener_proj140*ctot_ener_proj188*ctot_ener_proj54*ctot_ener_proj55;
const double ctot_ener_proj190 =             0.1381966*ctot_ener_proj189;
const double ctot_ener_proj191 =             ctot_ener_proj128*ctot_ener_proj140*ctot_ener_proj188*ctot_ener_proj54*ctot_ener_proj59;
const double ctot_ener_proj192 =             0.1381966*ctot_ener_proj191;
const double ctot_ener_proj193 =             ctot_ener_proj128*ctot_ener_proj152*ctot_ener_proj188*ctot_ener_proj54*ctot_ener_proj62;
const double ctot_ener_proj194 =             0.1381966*ctot_ener_proj193;
const double ctot_ener_proj195 =             ctot_ener_proj128*ctot_ener_proj152*ctot_ener_proj188*ctot_ener_proj54*ctot_ener_proj65;
const double ctot_ener_proj196 =             0.1381966*ctot_ener_proj195;
const double ctot_ener_proj197 =             ctot_ener_proj140*ctot_ener_proj152*ctot_ener_proj188*ctot_ener_proj54*ctot_ener_proj68;
const double ctot_ener_proj198 =             0.1381966*ctot_ener_proj197;
const double ctot_ener_proj199 =             ctot_ener_proj140*ctot_ener_proj152*ctot_ener_proj188*ctot_ener_proj54*ctot_ener_proj71;
const double ctot_ener_proj200 =             0.1381966*ctot_ener_proj199;
const double ctot_ener_proj201 =             pow(ctot_ener_proj119, -2);
const double ctot_ener_proj202 =             ctot_ener_proj131*ctot_ener_proj143*ctot_ener_proj201*ctot_ener_proj54*ctot_ener_proj55;
const double ctot_ener_proj203 =             ctot_ener_proj131*ctot_ener_proj143*ctot_ener_proj201*ctot_ener_proj54*ctot_ener_proj59;
const double ctot_ener_proj204 =             ctot_ener_proj131*ctot_ener_proj155*ctot_ener_proj201*ctot_ener_proj54*ctot_ener_proj62;
const double ctot_ener_proj205 =             ctot_ener_proj131*ctot_ener_proj155*ctot_ener_proj201*ctot_ener_proj54*ctot_ener_proj65;
const double ctot_ener_proj206 =             ctot_ener_proj143*ctot_ener_proj155*ctot_ener_proj201*ctot_ener_proj54*ctot_ener_proj68;
const double ctot_ener_proj207 =             ctot_ener_proj143*ctot_ener_proj155*ctot_ener_proj201*ctot_ener_proj54*ctot_ener_proj71;
const double ctot_ener_proj208 =             pow(ctot_ener_proj123, 2);
const double ctot_ener_proj209 =             1.0*ctot_ener_proj157*ctot_ener_proj54;
const double ctot_ener_proj210 =             0.1381966*U_1_4;
const double ctot_ener_proj211 =             -ctot_ener_proj210;
const double ctot_ener_proj212 =             0.5854102*U_3_4;
const double ctot_ener_proj213 =             ctot_ener_proj210 + ctot_ener_proj77;
const double ctot_ener_proj214 =             pow(ctot_ener_proj135, 2);
const double ctot_ener_proj215 =             pow(ctot_ener_proj147, 2);
const double ctot_ener_proj216 =             ctot_ener_proj211 - ctot_ener_proj212 - ctot_ener_proj54*(-ctot_ener_proj157*(0.5*ctot_ener_proj208 + 0.5*ctot_ener_proj214 + 0.5*ctot_ener_proj215) + ctot_ener_proj212 + ctot_ener_proj213 + ctot_ener_proj80) + ctot_ener_proj78 + ctot_ener_proj81;
const double ctot_ener_proj217 =             ctot_ener_proj157*ctot_ener_proj74*(ctot_ener_proj208*ctot_ener_proj209 + ctot_ener_proj216);
const double ctot_ener_proj218 =             0.1381966*ctot_ener_proj217;
const double ctot_ener_proj219 =             ctot_ener_proj157*ctot_ener_proj90*(ctot_ener_proj209*ctot_ener_proj214 + ctot_ener_proj216);
const double ctot_ener_proj220 =             0.1381966*ctot_ener_proj219;
const double ctot_ener_proj221 =             ctot_ener_proj157*ctot_ener_proj93*(ctot_ener_proj209*ctot_ener_proj215 + ctot_ener_proj216);
const double ctot_ener_proj222 =             0.1381966*ctot_ener_proj221;
const double ctot_ener_proj223 =             pow(ctot_ener_proj128, 2);
const double ctot_ener_proj224 =             1.0*ctot_ener_proj164*ctot_ener_proj54;
const double ctot_ener_proj225 =             0.5854102*U_2_4;
const double ctot_ener_proj226 =             pow(ctot_ener_proj140, 2);
const double ctot_ener_proj227 =             pow(ctot_ener_proj152, 2);
const double ctot_ener_proj228 =             ctot_ener_proj211 - ctot_ener_proj225 - ctot_ener_proj54*(-ctot_ener_proj164*(0.5*ctot_ener_proj223 + 0.5*ctot_ener_proj226 + 0.5*ctot_ener_proj227) + ctot_ener_proj213 + ctot_ener_proj225 + ctot_ener_proj82) + ctot_ener_proj78 + ctot_ener_proj83;
const double ctot_ener_proj229 =             ctot_ener_proj164*ctot_ener_proj74*(ctot_ener_proj223*ctot_ener_proj224 + ctot_ener_proj228);
const double ctot_ener_proj230 =             0.1381966*ctot_ener_proj229;
const double ctot_ener_proj231 =             ctot_ener_proj164*ctot_ener_proj90*(ctot_ener_proj224*ctot_ener_proj226 + ctot_ener_proj228);
const double ctot_ener_proj232 =             0.1381966*ctot_ener_proj231;
const double ctot_ener_proj233 =             ctot_ener_proj164*ctot_ener_proj93*(ctot_ener_proj224*ctot_ener_proj227 + ctot_ener_proj228);
const double ctot_ener_proj234 =             0.1381966*ctot_ener_proj233;
const double ctot_ener_proj235 =             pow(ctot_ener_proj131, 2);
const double ctot_ener_proj236 =             1.0*ctot_ener_proj171*ctot_ener_proj54;
const double ctot_ener_proj237 =             0.5854102*U_0_4;
const double ctot_ener_proj238 =             pow(ctot_ener_proj143, 2);
const double ctot_ener_proj239 =             pow(ctot_ener_proj155, 2);
const double ctot_ener_proj240 =             ctot_ener_proj211 - ctot_ener_proj237 - ctot_ener_proj54*(-ctot_ener_proj171*(0.5*ctot_ener_proj235 + 0.5*ctot_ener_proj238 + 0.5*ctot_ener_proj239) + ctot_ener_proj210 + ctot_ener_proj237 + ctot_ener_proj84) + ctot_ener_proj81 + ctot_ener_proj83;
const double ctot_ener_proj241 =             ctot_ener_proj171*ctot_ener_proj74*(ctot_ener_proj235*ctot_ener_proj236 + ctot_ener_proj240);
const double ctot_ener_proj242 =             ctot_ener_proj171*ctot_ener_proj90*(ctot_ener_proj236*ctot_ener_proj238 + ctot_ener_proj240);
const double ctot_ener_proj243 =             ctot_ener_proj171*ctot_ener_proj93*(ctot_ener_proj236*ctot_ener_proj239 + ctot_ener_proj240);
const double ctot_ener_proj244 =             ctot_ener_proj157*ctot_ener_proj97*(ctot_ener_proj208 + ctot_ener_proj214 + ctot_ener_proj215) + ctot_ener_proj216;
const double ctot_ener_proj245 =             ctot_ener_proj123*ctot_ener_proj175*ctot_ener_proj244*ctot_ener_proj96;
const double ctot_ener_proj246 =             -0.1381966*ctot_ener_proj245;
const double ctot_ener_proj247 =             ctot_ener_proj101*ctot_ener_proj135*ctot_ener_proj175*ctot_ener_proj244;
const double ctot_ener_proj248 =             -0.1381966*ctot_ener_proj247;
const double ctot_ener_proj249 =             ctot_ener_proj104*ctot_ener_proj147*ctot_ener_proj175*ctot_ener_proj244;
const double ctot_ener_proj250 =             -0.1381966*ctot_ener_proj249;
const double ctot_ener_proj251 =             ctot_ener_proj164*ctot_ener_proj97*(ctot_ener_proj223 + ctot_ener_proj226 + ctot_ener_proj227) + ctot_ener_proj228;
const double ctot_ener_proj252 =             ctot_ener_proj128*ctot_ener_proj188*ctot_ener_proj251*ctot_ener_proj96;
const double ctot_ener_proj253 =             -0.1381966*ctot_ener_proj252;
const double ctot_ener_proj254 =             ctot_ener_proj101*ctot_ener_proj140*ctot_ener_proj188*ctot_ener_proj251;
const double ctot_ener_proj255 =             -0.1381966*ctot_ener_proj254;
const double ctot_ener_proj256 =             ctot_ener_proj104*ctot_ener_proj152*ctot_ener_proj188*ctot_ener_proj251;
const double ctot_ener_proj257 =             -0.1381966*ctot_ener_proj256;
const double ctot_ener_proj258 =             ctot_ener_proj171*ctot_ener_proj97*(ctot_ener_proj235 + ctot_ener_proj238 + ctot_ener_proj239) + ctot_ener_proj240;
const double ctot_ener_proj259 =             ctot_ener_proj131*ctot_ener_proj201*ctot_ener_proj258*ctot_ener_proj96;
const double ctot_ener_proj260 =             ctot_ener_proj101*ctot_ener_proj143*ctot_ener_proj201*ctot_ener_proj258;
const double ctot_ener_proj261 =             ctot_ener_proj104*ctot_ener_proj155*ctot_ener_proj201*ctot_ener_proj258;
const double ctot_ener_proj262 =             0.1381966*ctot_ener_proj120;
const double ctot_ener_proj263 =             0.1381966*ctot_ener_proj132;
const double ctot_ener_proj264 =             0.1381966*ctot_ener_proj144;
const double ctot_ener_proj265 =             0.1381966*ctot_ener_proj156;
const double ctot_ener_proj266 =             -0.1381966*ctot_ener_proj172;
const double ctot_ener_proj267 =             -0.1381966*ctot_ener_proj173;
const double ctot_ener_proj268 =             -0.1381966*ctot_ener_proj174;
const double ctot_ener_proj269 =             0.1381966*ctot_ener_proj202;
const double ctot_ener_proj270 =             0.1381966*ctot_ener_proj203;
const double ctot_ener_proj271 =             0.1381966*ctot_ener_proj204;
const double ctot_ener_proj272 =             0.1381966*ctot_ener_proj205;
const double ctot_ener_proj273 =             0.1381966*ctot_ener_proj206;
const double ctot_ener_proj274 =             0.1381966*ctot_ener_proj207;
const double ctot_ener_proj275 =             0.1381966*ctot_ener_proj241;
const double ctot_ener_proj276 =             0.1381966*ctot_ener_proj242;
const double ctot_ener_proj277 =             0.1381966*ctot_ener_proj243;
const double ctot_ener_proj278 =             -0.1381966*ctot_ener_proj259;
const double ctot_ener_proj279 =             -0.1381966*ctot_ener_proj260;
const double ctot_ener_proj280 =             -0.1381966*ctot_ener_proj261;
const double ctot_ener_proj281 =             ctot_ener_proj10 + ctot_ener_proj100 + ctot_ener_proj103 + ctot_ener_proj106 + ctot_ener_proj21 + ctot_ener_proj262 + ctot_ener_proj263 + ctot_ener_proj264 + ctot_ener_proj265 + ctot_ener_proj266 + ctot_ener_proj267 + ctot_ener_proj268 + ctot_ener_proj269 + ctot_ener_proj270 + ctot_ener_proj271 + ctot_ener_proj272 + ctot_ener_proj273 + ctot_ener_proj274 + ctot_ener_proj275 + ctot_ener_proj276 + ctot_ener_proj277 + ctot_ener_proj278 + ctot_ener_proj279 + ctot_ener_proj280 + ctot_ener_proj32 + ctot_ener_proj43 + ctot_ener_proj47 + ctot_ener_proj50 + ctot_ener_proj53 + ctot_ener_proj58 + ctot_ener_proj61 + ctot_ener_proj64 + ctot_ener_proj67 + ctot_ener_proj70 + ctot_ener_proj73 + ctot_ener_proj89 + ctot_ener_proj92 + ctot_ener_proj95 - 0.19999999899376*dUdt_0_4 - 0.19999999899376*dUdt_1_4;
            tot_ener_proj[0]=ctot_ener_proj10 + ctot_ener_proj100 + ctot_ener_proj103 + ctot_ener_proj106 + ctot_ener_proj107 + ctot_ener_proj108 + ctot_ener_proj115 + ctot_ener_proj118 + 0.5854102*ctot_ener_proj120 + ctot_ener_proj127 + ctot_ener_proj130 + 0.5854102*ctot_ener_proj132 + ctot_ener_proj139 + ctot_ener_proj142 + 0.5854102*ctot_ener_proj144 + ctot_ener_proj151 + ctot_ener_proj154 + 0.5854102*ctot_ener_proj156 + ctot_ener_proj159 + ctot_ener_proj161 + ctot_ener_proj163 + ctot_ener_proj166 + ctot_ener_proj168 + ctot_ener_proj170 - 0.5854102*ctot_ener_proj172 - 0.5854102*ctot_ener_proj173 - 0.5854102*ctot_ener_proj174 + ctot_ener_proj177 + ctot_ener_proj179 + ctot_ener_proj181 + ctot_ener_proj183 + ctot_ener_proj185 + ctot_ener_proj187 + ctot_ener_proj190 + ctot_ener_proj192 + ctot_ener_proj194 + ctot_ener_proj196 + ctot_ener_proj198 + ctot_ener_proj200 + 0.5854102*ctot_ener_proj202 + 0.5854102*ctot_ener_proj203 + 0.5854102*ctot_ener_proj204 + 0.5854102*ctot_ener_proj205 + 0.5854102*ctot_ener_proj206 + 0.5854102*ctot_ener_proj207 + ctot_ener_proj21 + ctot_ener_proj218 + ctot_ener_proj220 + ctot_ener_proj222 + ctot_ener_proj230 + ctot_ener_proj232 + ctot_ener_proj234 + 0.5854102*ctot_ener_proj241 + 0.5854102*ctot_ener_proj242 + 0.5854102*ctot_ener_proj243 + ctot_ener_proj246 + ctot_ener_proj248 + ctot_ener_proj250 + ctot_ener_proj253 + ctot_ener_proj255 + ctot_ener_proj257 - 0.5854102*ctot_ener_proj259 - 0.5854102*ctot_ener_proj260 - 0.5854102*ctot_ener_proj261 + ctot_ener_proj32 + ctot_ener_proj43 + ctot_ener_proj47 + ctot_ener_proj50 + ctot_ener_proj53 + ctot_ener_proj58 + ctot_ener_proj61 + ctot_ener_proj64 + ctot_ener_proj67 + ctot_ener_proj70 + ctot_ener_proj73 + ctot_ener_proj89 + ctot_ener_proj92 + ctot_ener_proj95 - 0.40000000301872*dUdt_0_4 - 0.19999999899376*dUdt_1_4;
            tot_ener_proj[1]=-0.5854102*ctot_ener_proj102 - 0.5854102*ctot_ener_proj105 + ctot_ener_proj107 + ctot_ener_proj108 + ctot_ener_proj115 + ctot_ener_proj118 + ctot_ener_proj127 + ctot_ener_proj130 + ctot_ener_proj139 + ctot_ener_proj142 + ctot_ener_proj151 + ctot_ener_proj154 + ctot_ener_proj159 + ctot_ener_proj161 + ctot_ener_proj163 + ctot_ener_proj166 + ctot_ener_proj168 + ctot_ener_proj170 + ctot_ener_proj177 + ctot_ener_proj179 + ctot_ener_proj181 + ctot_ener_proj183 + ctot_ener_proj185 + ctot_ener_proj187 + ctot_ener_proj190 + ctot_ener_proj192 + ctot_ener_proj194 + ctot_ener_proj196 + ctot_ener_proj198 + 0.5854102*ctot_ener_proj20 + ctot_ener_proj200 + ctot_ener_proj218 + ctot_ener_proj220 + ctot_ener_proj222 + ctot_ener_proj230 + ctot_ener_proj232 + ctot_ener_proj234 + ctot_ener_proj246 + ctot_ener_proj248 + ctot_ener_proj250 + ctot_ener_proj253 + ctot_ener_proj255 + ctot_ener_proj257 + ctot_ener_proj262 + ctot_ener_proj263 + ctot_ener_proj264 + ctot_ener_proj265 + ctot_ener_proj266 + ctot_ener_proj267 + ctot_ener_proj268 + ctot_ener_proj269 + ctot_ener_proj270 + ctot_ener_proj271 + ctot_ener_proj272 + ctot_ener_proj273 + ctot_ener_proj274 + ctot_ener_proj275 + ctot_ener_proj276 + ctot_ener_proj277 + ctot_ener_proj278 + ctot_ener_proj279 + ctot_ener_proj280 + 0.5854102*ctot_ener_proj31 + 0.5854102*ctot_ener_proj42 - 0.5854102*ctot_ener_proj46 - 0.5854102*ctot_ener_proj49 - 0.5854102*ctot_ener_proj52 + 0.5854102*ctot_ener_proj57 + 0.5854102*ctot_ener_proj60 + 0.5854102*ctot_ener_proj63 + 0.5854102*ctot_ener_proj66 + 0.5854102*ctot_ener_proj69 + 0.5854102*ctot_ener_proj72 + 0.5854102*ctot_ener_proj88 + 0.5854102*ctot_ener_proj9 + 0.5854102*ctot_ener_proj91 + 0.5854102*ctot_ener_proj94 - 0.5854102*ctot_ener_proj99 - 0.19999999899376*dUdt_0_4 - 0.40000000301872*dUdt_1_4;
            tot_ener_proj[2]=ctot_ener_proj108 + ctot_ener_proj115 + 0.5854102*ctot_ener_proj117 + ctot_ener_proj127 + 0.5854102*ctot_ener_proj129 + ctot_ener_proj139 + 0.5854102*ctot_ener_proj141 + ctot_ener_proj151 + 0.5854102*ctot_ener_proj153 + ctot_ener_proj159 + ctot_ener_proj161 + ctot_ener_proj163 - 0.5854102*ctot_ener_proj165 - 0.5854102*ctot_ener_proj167 - 0.5854102*ctot_ener_proj169 + ctot_ener_proj177 + ctot_ener_proj179 + ctot_ener_proj181 + ctot_ener_proj183 + ctot_ener_proj185 + ctot_ener_proj187 + 0.5854102*ctot_ener_proj189 + 0.5854102*ctot_ener_proj191 + 0.5854102*ctot_ener_proj193 + 0.5854102*ctot_ener_proj195 + 0.5854102*ctot_ener_proj197 + 0.5854102*ctot_ener_proj199 + ctot_ener_proj218 + ctot_ener_proj220 + ctot_ener_proj222 + 0.5854102*ctot_ener_proj229 + 0.5854102*ctot_ener_proj231 + 0.5854102*ctot_ener_proj233 + ctot_ener_proj246 + ctot_ener_proj248 + ctot_ener_proj250 - 0.5854102*ctot_ener_proj252 - 0.5854102*ctot_ener_proj254 - 0.5854102*ctot_ener_proj256 + ctot_ener_proj281 - 0.40000000301872*dUdt_2_4;
            tot_ener_proj[3]=ctot_ener_proj107 + 0.5854102*ctot_ener_proj114 + ctot_ener_proj118 + 0.5854102*ctot_ener_proj126 + ctot_ener_proj130 + 0.5854102*ctot_ener_proj138 + ctot_ener_proj142 + 0.5854102*ctot_ener_proj150 + ctot_ener_proj154 - 0.5854102*ctot_ener_proj158 - 0.5854102*ctot_ener_proj160 - 0.5854102*ctot_ener_proj162 + ctot_ener_proj166 + ctot_ener_proj168 + ctot_ener_proj170 + 0.5854102*ctot_ener_proj176 + 0.5854102*ctot_ener_proj178 + 0.5854102*ctot_ener_proj180 + 0.5854102*ctot_ener_proj182 + 0.5854102*ctot_ener_proj184 + 0.5854102*ctot_ener_proj186 + ctot_ener_proj190 + ctot_ener_proj192 + ctot_ener_proj194 + ctot_ener_proj196 + ctot_ener_proj198 + ctot_ener_proj200 + 0.5854102*ctot_ener_proj217 + 0.5854102*ctot_ener_proj219 + 0.5854102*ctot_ener_proj221 + ctot_ener_proj230 + ctot_ener_proj232 + ctot_ener_proj234 - 0.5854102*ctot_ener_proj245 - 0.5854102*ctot_ener_proj247 - 0.5854102*ctot_ener_proj249 + ctot_ener_proj253 + ctot_ener_proj255 + ctot_ener_proj257 + ctot_ener_proj281 - 0.40000000301872*dUdt_3_4;

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    tot_ener_proj *= data.volume / static_cast<double>(n_nodes);

    // Assembly the projection contributions
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
#pragma omp atomic
        r_geometry[i_node].GetValue(TOTAL_ENERGY_PROJECTION) += tot_ener_proj[i_node];
    }

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::CalculateRightHandSideInternal(
    BoundedVector<double, 15> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned  int nodesElement = 3;
    constexpr unsigned  int SpaceDimension = 2;
    constexpr unsigned  int nScalarVariables = SpaceDimension + 3;
    constexpr unsigned  int nNodalVariables = nScalarVariables*nodesElement;

    unsigned int i, j, k, l, m, s, t, tt, p, pp;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    const unsigned int size3 = nScalarVariables * SpaceDimension * nScalarVariables;
    const unsigned int size4 = nScalarVariables * SpaceDimension * nScalarVariables * nScalarVariables;
    const unsigned int sizeK = SpaceDimension * SpaceDimension * nScalarVariables * SpaceDimension;
    const unsigned int sizeKT = nScalarVariables;

    const BoundedMatrix<double,nodesElement,nScalarVariables>& UU = data.U;			
    //const BoundedMatrix<double,nodesElement,nScalarVariables>& UUn = data.Un;
    const BoundedMatrix<double,nodesElement,nScalarVariables>& Up = data.dUdt;    // Useful for the stabilizing part.
    BoundedMatrix<double,nodesElement,nScalarVariables> UUp = Up;

    array_1d<double, size3>     A;
    array_1d<double, size4>     dAdU;
    array_1d<double, sizeK>     K;
    array_1d<double, sizeKT>    KT;

    array_1d<double,nNodalVariables>                     LumpedMassMatrix;
    array_1d<double,nScalarVariables>                    U_gauss;
    array_1d<double,nNodalVariables>                     U;
    array_1d<double,nNodalVariables>                     Un;
    array_1d<double,nNodalVariables>                     up;
    array_1d<double,nScalarVariables>                    L;
    array_1d<double,nScalarVariables>                    Residual;
    array_1d<double,nNodalVariables*nScalarVariables>    Lstar;
    array_1d<double,nScalarVariables*SpaceDimension>     G;
    array_1d<double,nScalarVariables*nScalarVariables>   S;
    array_1d<double,nScalarVariables*nScalarVariables>   B;
    array_1d<double,SpaceDimension>                      dt_diff;
	array_1d<double,SpaceDimension>                      ds_diff;
	array_1d<double,SpaceDimension*SpaceDimension>       tau;
    array_1d<double,SpaceDimension>                      q;
	array_1d<double,SpaceDimension>                      a_el;
    array_1d<double,nScalarVariables*nNodalVariables>    NN;
    array_1d<double,nScalarVariables*SpaceDimension>     gradU;
    array_1d<double,(nScalarVariables*SpaceDimension)*nNodalVariables>   gradV;
    array_1d<double,nScalarVariables>                    invtauStab;
	array_1d<double,nScalarVariables>                    switchStab;

    array_1d<double,nNodalVariables>     FConv;
    array_1d<double,nNodalVariables>     FStab;
    array_1d<double,nNodalVariables>     FDiff;
    array_1d<double,nNodalVariables>     F;
    
	const double& ctau = 0.5;   // This coefficient multiplies the divergence of the velocity in the calculation of tau. 
                                // In 3d would be 0.66667
    
    const double& dt = data.time_step;


    const double h = data.h;
    
    const BoundedMatrix<double,nodesElement,SpaceDimension>& f_ext = data.f_ext;			
    const array_1d<double, nodesElement>& r_ext = data.r_ext;
    const array_1d<double, nodesElement>& m_ext = data.m_ext;
    
    const array_1d<double,nodesElement>& r = data.r_ext;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double Cv = data.c_v;
    const double gamma = data.gamma;
    const double Cp = Cv * gamma;
    const double R = Cp - Cv;
    const double ros = data.ros;
    const double Cs = data.c_s;

	const double sw_conv = 1.0;
    const double sw_diff = 1.0;
    const double sw_stab = 1.0;


    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    switchStab[0] = 1.0;
	switchStab[1] = 1.0;
	switchStab[2] = 1.0;
	switchStab[3] = 1.0;
	switchStab[4] = 1.0;    

    const unsigned int NGaussPoints = 1;

    const array_1d<double,nodesElement>& N = data.N;	
    const BoundedMatrix<double,nodesElement,SpaceDimension>& DN = data.DN_DX;	
    
    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,SpaceDimension> f_gauss = prod(trans(f_ext), N);      
    const double r_gauss = N(0)*r_ext(0) + N(1)*r_ext(1) + N(2)*r_ext(2);
    

    // Stabilization parameters
   for (i = 0; i < nScalarVariables; i++){
        U_gauss(i) = 0;
        for (j = 0; j < nodesElement; j++){
            U_gauss(i) += N(j)*UU(j,i);

        }
    }

    for (i = 0; i < nScalarVariables*SpaceDimension; i++)       gradU(i) = 0.0;
// Compute the gradient of the variables /dro/dx,dro/dy, dm1/dx dm1/dy dm2/dx dm2/dy de/dx de/dy)
    for (i = 0; i < nScalarVariables; i++){         // Controllare
        for (k = 0; k < nodesElement; k++){

            gradU(i*SpaceDimension + 0) += DN(k,0)*UU(k,i);
            gradU(i*SpaceDimension + 1) += DN(k,1)*UU(k,i);

        }
    }

    // This is convenient during the implementation but has to be removed 
    // after some modification of the remainding part of the file
   
    // This is convenient during the implementation but has to be removed 
    // after some modification of the remainding part of the file    
    for ( i = 0; i < nScalarVariables * SpaceDimension * nNodalVariables; i++)  gradV[i] = 0.0;

	for (i = 0; i < nodesElement; i++ ){
        
        gradV[0*nNodalVariables + nScalarVariables*i    ] = DN(i,0);
        gradV[1*nNodalVariables + nScalarVariables*i    ] = DN(i,1);

        gradV[2*nNodalVariables + nScalarVariables*i + 1] = DN(i,0);
        gradV[3*nNodalVariables + nScalarVariables*i + 1] = DN(i,1);

        gradV[4*nNodalVariables + nScalarVariables*i + 2] = DN(i,0);
        gradV[5*nNodalVariables + nScalarVariables*i + 2] = DN(i,1);

        gradV[6*nNodalVariables + nScalarVariables*i + 3] = DN(i,0);
        gradV[7*nNodalVariables + nScalarVariables*i + 3] = DN(i,1);

        gradV[8*nNodalVariables + nScalarVariables*i + 4] = DN(i,0);  // Verify these lines
        gradV[9*nNodalVariables + nScalarVariables*i + 4] = DN(i,1);

	}


    const double DTOT_el = U_gauss(0);
    const double DS_el = U_gauss(1);
    const double m1_el = U_gauss(2);
    const double m2_el = U_gauss(3);
    const double etot_el = U_gauss(4);

    const double DG_el = DTOT_el - DS_el; 

    const double  kk = DS_el/DG_el;

    // HERE DEFINE MIXTURE COEFFICIENTS CV_M, ECC.

    const double mu_mixture = mu/(1 + kk);        // Controllare questo spaccimmo di termine (trovare su Neri et al.)

	const double Cp_mixture = (Cp + kk*Cs)/(1.0 + kk);
	const double Cv_mixture = (Cv + kk*Cs)/(1.0 + kk);

	const double Prandtl = 0.5;

//	const double lambda_mixture = Cp_mixture*mu_mixture/Prandtl;
	const double lambda_mixture = lambda;

	// Fine variazione

    const double u1_el = m1_el/DTOT_el;
	const double u2_el = m2_el/DTOT_el;

	const double norm2u = u1_el*u1_el + u2_el*u2_el;
	const double norm_u = sqrt(norm2u);
	const double norm2m = m1_el*m1_el + m2_el*m2_el;

	a_el(0) = u1_el/norm_u;
	a_el(1) = u2_el/norm_u;

	if (norm_u == 0){
		a_el(0) = 0.0;
		a_el(1) = 0.0;	 
	}	



    const double	DG_el2  = DG_el*DG_el;
	const double	DTOT2	= DTOT_el*DTOT_el;
	const double	DTOT3	= DTOT_el*DTOT2;
	
	const double	Cmixed  = Cs*DS_el + Cv*DG_el;
	const double	Cmixed2 = Cmixed*Cmixed;
	const double	Cmixed3 = Cmixed*Cmixed2;
	const double	Cpmixed = Cs*DS_el + (Cv + R)*DG_el;
	
	const double 	Cs2Ds3	= Cs*Cs*DS_el*DS_el*DS_el;
	const double 	Cv2Dg3	= Cv*Cv*DG_el*DG_el*DG_el;
	
	const double 	a_d 	= 0;
	
	const double	p_el   	= (DG_el*(2*DTOT_el*etot_el - norm2m)*R)/(2*DTOT_el*Cmixed);
		
	const double	pdg 	=  ((Cv*DG_el2*norm2m - Cs*DS_el*(-2*DTOT2*etot_el + DS_el*norm2m))*R)/(2*DTOT2*Cmixed2);
	const double	pds 	=  (Cs*(-2*DTOT_el*etot_el + norm2m)*R)/(2*Cmixed2);
	const double	pm1 	= -(DG_el*m1_el*R)/(DTOT_el*Cmixed);
	const double	pm2 	= -(DG_el*m2_el*R)/(DTOT_el*Cmixed);
	const double	pet 	= DG_el*R/Cmixed;
	
	const double	pdgdg 	=  -(2*Cs*Cv*DS_el*etot_el*R)/Cmixed3 + ((Cs2Ds3 + Cv2Dg3 + Cs*Cv*DS_el*DS_el*(-2*DS_el + 3*DTOT_el))*norm2m*R)/(DTOT3*Cmixed3);
	const double	pdgds 	=  (Cs*(-Cs*DS_el + Cv*(DS_el + DTOT_el))*etot_el*R)/Cmixed3 - (Cs*Cv*norm2m*R)/Cmixed3;
	const double	pdgm1 	=  ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*m1_el*R)/(DTOT2*Cmixed2);
	const double	pdgm2 	=  ((Cv*DG_el*DG_el - Cs*DS_el*DS_el)*m2_el*R)/(DTOT2*Cmixed2);
	const double	pdget 	=  Cs*DS_el*R/Cmixed2;
	
	const double	pdsds 	=  (2*Cs*(Cs - Cv)*DTOT_el*etot_el*R)/Cmixed3 + (Cs*(-Cs + Cv)*norm2m*R)/Cmixed3;
	const double	pdsm1 	=  Cs*m1_el*R/Cmixed2;
	const double	pdsm2 	=  Cs*m2_el*R/Cmixed2;
	const double	pdset 	= -Cs*DTOT_el*R/Cmixed2;
	
	const double	pm1m1 	= -DG_el*R/(DTOT_el*Cmixed);
	
	const double	pm2m2 	= -DG_el*R/(DTOT_el*Cmixed);
	
	const double    Temperature = p_el/(DG_el*R);

	const double gas_concentration = 1 - DS_el/ros;
	const double gas_density = DG_el/gas_concentration;

	// printf("eps_g = %.3e\n", gas_concentration);

    double SpeedSound2  = DG_el*R*Cpmixed*(2*etot_el - DTOT_el*norm2u)/(2*DTOT_el*Cmixed2);  // verificare sound_speed
	double SpeedSound;

	SpeedSound = sqrt(SpeedSound2);
		 

	// printf("gamma = %.3e - SpeedSound = %.3e\n", gamma, SpeedSound); 

	if (dt > h/(SpeedSound + norm_u))   printf("dt = %.3e dt_crit = %.3e\n", dt, h/(SpeedSound + norm_u));

	if (DTOT_el < 0)    				printf("dtot_el = %.3e \n", DTOT_el);

    for (i = 0; i < size3; i++)     A(i) = 0.0;
	for (i = 0; i < size4; i++)     dAdU(i) = 0.0;
			
	// Build A

	A[conta_new(0,0,2,0,5,2,5,1)] = 1.0;
	A[conta_new(0,1,3,0,5,2,5,1)] = 1.0;
	
	A[conta_new(1,0,0,0,5,2,5,1)] = -DS_el*m1_el/DTOT2; 
	A[conta_new(1,0,1,0,5,2,5,1)] =  m1_el/DTOT_el;
	A[conta_new(1,0,2,0,5,2,5,1)] =  DS_el/DTOT_el;
	
	A[conta_new(1,1,0,0,5,2,5,1)] = -DS_el*m2_el/DTOT2; 
	A[conta_new(1,1,1,0,5,2,5,1)] =  m2_el/DTOT_el;
	A[conta_new(1,1,3,0,5,2,5,1)] =  DS_el/DTOT_el;
	
	A[conta_new(2,0,0,0,5,2,5,1)] =  -m1_el*m1_el/DTOT2 + pdg;
	A[conta_new(2,0,1,0,5,2,5,1)] =  pds;
	A[conta_new(2,0,2,0,5,2,5,1)] =  2*m1_el/DTOT_el + pm1;
	A[conta_new(2,0,3,0,5,2,5,1)] =  pm2;
	A[conta_new(2,0,4,0,5,2,5,1)] =  pet;
	
	A[conta_new(2,1,0,0,5,2,5,1)] =  -m1_el*m2_el/DTOT2;
	A[conta_new(2,1,2,0,5,2,5,1)] =  m2_el/DTOT_el;
	A[conta_new(2,1,3,0,5,2,5,1)] =  m1_el/DTOT_el;
	
	A[conta_new(3,0,0,0,5,2,5,1)] =  -m1_el*m2_el/DTOT2;
	A[conta_new(3,0,2,0,5,2,5,1)] =  m2_el/DTOT_el;
	A[conta_new(3,0,3,0,5,2,5,1)] =  m1_el/DTOT_el;
		
	A[conta_new(3,1,0,0,5,2,5,1)] =  -m2_el*m2_el/DTOT2 + pdg;
	A[conta_new(3,1,1,0,5,2,5,1)] =  pds;
	A[conta_new(3,1,2,0,5,2,5,1)] =  pm1;
	A[conta_new(3,1,3,0,5,2,5,1)] =  2*m2_el/DTOT_el + pm2;
	A[conta_new(3,1,4,0,5,2,5,1)] =  pet;
	
	A[conta_new(4,0,0,0,5,2,5,1)] =  -(m1_el*(etot_el - DTOT_el*pdg + p_el))/DTOT2;
	A[conta_new(4,0,1,0,5,2,5,1)] =  m1_el*pds/DTOT_el;
	A[conta_new(4,0,2,0,5,2,5,1)] =  (etot_el + m1_el*pm1 + p_el)/DTOT_el;
	A[conta_new(4,0,3,0,5,2,5,1)] =  m1_el*pm2/DTOT_el;
	A[conta_new(4,0,4,0,5,2,5,1)] =	m1_el*(1.0 + pet)/DTOT_el;
	
	A[conta_new(4,1,0,0,5,2,5,1)] =  -(m2_el*(etot_el - DTOT_el*pdg + p_el))/DTOT2;
	A[conta_new(4,1,1,0,5,2,5,1)] =  m2_el*pds/DTOT_el;
	A[conta_new(4,1,2,0,5,2,5,1)] =	m2_el*pm1/DTOT_el;  
	A[conta_new(4,1,3,0,5,2,5,1)] =  (etot_el + m2_el*pm2 + p_el)/DTOT_el;
	A[conta_new(4,1,4,0,5,2,5,1)] =  m2_el*(1.0 + pet)/DTOT_el;


    //	Build dAdU

// 1  ////////////////////////////////////////////////////////
			
	dAdU[conta_new(1,0,0,0,5,2,5,5)] = 2*DS_el*m1_el/DTOT3; 
	dAdU[conta_new(1,0,0,1,5,2,5,5)] = -m1_el/DTOT2;
	dAdU[conta_new(1,0,0,2,5,2,5,5)] = -DS_el/DTOT2;
	
	dAdU[conta_new(1,0,1,0,5,2,5,5)] = -m1_el/DTOT2;  
	dAdU[conta_new(1,0,1,2,5,2,5,5)] = 1.0/DTOT_el;
	 
	dAdU[conta_new(1,0,2,0,5,2,5,5)] = -DS_el/DTOT2; 
	dAdU[conta_new(1,0,2,1,5,2,5,5)] =  1.0/DTOT_el;
	
/////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(1,1,0,0,5,2,5,5)] =  2*DS_el*m2_el/DTOT3;
	dAdU[conta_new(1,1,0,1,5,2,5,5)] =  -m2_el/DTOT2;
	dAdU[conta_new(1,1,0,3,5,2,5,5)] = -DS_el/DTOT2; 
	
	dAdU[conta_new(1,1,1,0,5,2,5,5)] =  -m2_el/DTOT2;
	dAdU[conta_new(1,1,1,3,5,2,5,5)] = 1.0/DTOT_el; 
	 
	dAdU[conta_new(1,1,3,0,5,2,5,5)] = -DS_el/DTOT2;
	dAdU[conta_new(1,1,3,1,5,2,5,5)] =  1.0/DTOT_el;
				
// 2  /////////////////////////////////////////////////////////////
			
	dAdU[conta_new(2,0,0,0,5,2,5,5)] =  2*m1_el*m1_el/DTOT3 + pdgdg;	 
	dAdU[conta_new(2,0,0,1,5,2,5,5)] =  pdgds;
	dAdU[conta_new(2,0,0,2,5,2,5,5)] = -2*m1_el/DTOT2 + pdgm1;
	dAdU[conta_new(2,0,0,3,5,2,5,5)] =  pdgm2;
	dAdU[conta_new(2,0,0,4,5,2,5,5)] =  pdget;
	
	dAdU[conta_new(2,0,1,0,5,2,5,5)] = pdgds;
	dAdU[conta_new(2,0,1,1,5,2,5,5)] = pdsds;
	dAdU[conta_new(2,0,1,2,5,2,5,5)] = pdsm1;
	dAdU[conta_new(2,0,1,3,5,2,5,5)] = pdsm2;
	dAdU[conta_new(2,0,1,4,5,2,5,5)] = pdset;
	 
	dAdU[conta_new(2,0,2,0,5,2,5,5)] = -2*m1_el/DTOT2 + pdgm1;
	dAdU[conta_new(2,0,2,1,5,2,5,5)] = pdsm1;
	dAdU[conta_new(2,0,2,2,5,2,5,5)] = 2.0/DTOT_el + pm1m1;
	 
	dAdU[conta_new(2,0,3,0,5,2,5,5)] = pdgm2; 
	dAdU[conta_new(2,0,3,1,5,2,5,5)] = pdsm2;
	dAdU[conta_new(2,0,3,3,5,2,5,5)] = pm2m2; 
				
	dAdU[conta_new(2,0,4,0,5,2,5,5)] = pdget;
	dAdU[conta_new(2,0,4,1,5,2,5,5)] = pdset;
	 
	dAdU[conta_new(2,1,0,0,5,2,5,5)] = 2*m1_el*m2_el/DTOT3;	 
	dAdU[conta_new(2,1,0,2,5,2,5,5)] = -m2_el/DTOT2;
	dAdU[conta_new(2,1,0,3,5,2,5,5)] = -m1_el/DTOT2;
	
	dAdU[conta_new(2,1,2,0,5,2,5,5)] = -m2_el/DTOT2; 
	dAdU[conta_new(2,1,2,3,5,2,5,5)] = 1.0/DTOT_el; 
			
	dAdU[conta_new(2,1,3,0,5,2,5,5)] = -m1_el/DTOT2; 
	dAdU[conta_new(2,1,3,2,5,2,5,5)] = 1.0/DTOT_el;
	
// ARRIVATO QUA
// 3 /////////////////////////////////////////////////////////////////////////////
			
	dAdU[conta_new(3,0,0,0,5,2,5,5)] = 2*m1_el*m2_el/DTOT3;	 
	dAdU[conta_new(3,0,0,2,5,2,5,5)] = -m2_el/DTOT2;
	dAdU[conta_new(3,0,0,3,5,2,5,5)] = -m1_el/DTOT2;
		
	dAdU[conta_new(3,0,2,0,5,2,5,5)] = -m2_el/DTOT2; 
	dAdU[conta_new(3,0,2,3,5,2,5,5)] = 1.0/DTOT_el; 
			
	dAdU[conta_new(3,0,3,0,5,2,5,5)] = -m1_el/DTOT2; 
	dAdU[conta_new(3,0,3,2,5,2,5,5)] = 1.0/DTOT_el;
	
//////////////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(3,1,0,0,5,2,5,5)] = 2*m2_el*m2_el/DTOT3 + pdgdg;	 
	dAdU[conta_new(3,1,0,1,5,2,5,5)] = pdgds;
	dAdU[conta_new(3,1,0,2,5,2,5,5)] = pdgm1;
	dAdU[conta_new(3,1,0,3,5,2,5,5)] = -2*m2_el/DTOT2 + pdgm2;
	dAdU[conta_new(3,1,0,4,5,2,5,5)] = pdget;
	
	dAdU[conta_new(3,1,1,0,5,2,5,5)] = pdgds;
	dAdU[conta_new(3,1,1,1,5,2,5,5)] = pdsds;
	dAdU[conta_new(3,1,1,2,5,2,5,5)] = pdsm1;
	dAdU[conta_new(3,1,1,3,5,2,5,5)] = pdsm2;
	dAdU[conta_new(3,1,1,4,5,2,5,5)] = pdset;
	 
	dAdU[conta_new(3,1,2,0,5,2,5,5)] = pdgm1;
	dAdU[conta_new(3,1,2,1,5,2,5,5)] = pdsm1;
	dAdU[conta_new(3,1,2,2,5,2,5,5)] = pm1m1;
	
	dAdU[conta_new(3,1,3,0,5,2,5,5)] = -2*m2_el/DTOT2 + pdgm2;
	dAdU[conta_new(3,1,3,1,5,2,5,5)] = pdsm2;
	dAdU[conta_new(3,1,3,3,5,2,5,5)] = 2.0/DTOT_el + pm2m2; 
	
	dAdU[conta_new(3,1,4,0,5,2,5,5)] = pdget;
	dAdU[conta_new(3,1,4,1,5,2,5,5)] = pdset;
	dAdU[conta_new(3,1,2,1,5,2,5,5)] = pdsm1;
	dAdU[conta_new(3,1,2,2,5,2,5,5)] = pm1m1;
	
	dAdU[conta_new(3,1,3,0,5,2,5,5)] = -2*m2_el/DTOT2 + pdgm2;
	dAdU[conta_new(3,1,3,1,5,2,5,5)] = pdsm2;
	dAdU[conta_new(3,1,3,3,5,2,5,5)] = 2.0/DTOT_el + pm2m2; 
	
	dAdU[conta_new(3,1,4,0,5,2,5,5)] = pdget;
	dAdU[conta_new(3,1,4,1,5,2,5,5)] = pdset;
			
// 4 //////////////////////////////////////////////////////////////////////
			
	dAdU[conta_new(4,0,0,0,5,2,5,5)] = (m1_el*(2*etot_el - 2*DTOT_el*pdg + DTOT_el*DTOT_el*pdgdg + 2*p_el))/DTOT3;
	dAdU[conta_new(4,0,0,1,5,2,5,5)] = (m1_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(4,0,0,2,5,2,5,5)] = -((etot_el - DTOT_el*(pdg + m1_el*pdgm1) + m1_el*pm1 + p_el)/DTOT2);
	dAdU[conta_new(4,0,0,3,5,2,5,5)] = (m1_el*(DTOT_el*pdgm2 - pm2))/DTOT2;
	dAdU[conta_new(4,0,0,4,5,2,5,5)] = (m1_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	
	dAdU[conta_new(4,0,1,0,5,2,5,5)] = (m1_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(4,0,1,1,5,2,5,5)] = (m1_el*pdsds)/DTOT_el;
	dAdU[conta_new(4,0,1,2,5,2,5,5)] = (pds + m1_el*pdsm1)/DTOT_el;
	dAdU[conta_new(4,0,1,3,5,2,5,5)] = (m1_el*pdsm2)/DTOT_el;
	dAdU[conta_new(4,0,1,4,5,2,5,5)] = (m1_el*pdset)/DTOT_el;
	 
	dAdU[conta_new(4,0,2,0,5,2,5,5)] = -((etot_el - DTOT_el*(pdg + m1_el*pdgm1) + m1_el*pm1 + p_el)/DTOT2);
	dAdU[conta_new(4,0,2,1,5,2,5,5)] = (pds + m1_el*pdsm1)/DTOT_el;
	dAdU[conta_new(4,0,2,2,5,2,5,5)] = (2*pm1 + m1_el*pm1m1)/DTOT_el;
	dAdU[conta_new(4,0,2,3,5,2,5,5)] = pm2/DTOT_el;
	dAdU[conta_new(4,0,2,4,5,2,5,5)] = (1 + pet)/DTOT_el;
			
	dAdU[conta_new(4,0,3,0,5,2,5,5)] = (m1_el*(DTOT_el*pdgm2 - pm2))/DTOT2;
	dAdU[conta_new(4,0,3,1,5,2,5,5)] = (m1_el*pdsm2)/DTOT_el;
	dAdU[conta_new(4,0,3,2,5,2,5,5)] = pm2/DTOT_el;
	dAdU[conta_new(4,0,3,3,5,2,5,5)] = (m1_el*pm2m2)/DTOT_el;
			
	dAdU[conta_new(4,0,4,0,5,2,5,5)] = (m1_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	dAdU[conta_new(4,0,4,1,5,2,5,5)] = (m1_el*pdset)/DTOT_el;
	dAdU[conta_new(4,0,4,2,5,2,5,5)] = (1 + pet)/DTOT_el;
	 
	dAdU[conta_new(4,1,0,0,5,2,5,5)] = (m2_el*(2*etot_el - 2*DTOT_el*pdg + DTOT2*pdgdg + 2*p_el))/DTOT3;	 
	dAdU[conta_new(4,1,0,1,5,2,5,5)] = (m2_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(4,1,0,2,5,2,5,5)] = (m2_el*(DTOT_el*pdgm1 - pm1))/DTOT2; 
	dAdU[conta_new(4,1,0,3,5,2,5,5)] = -((etot_el - DTOT_el*(pdg + m2_el*pdgm2) + m2_el*pm2 + p_el)/DTOT2);
	dAdU[conta_new(4,1,0,4,5,2,5,5)] = (m2_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	
	dAdU[conta_new(4,1,1,0,5,2,5,5)] = (m2_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(4,1,1,1,5,2,5,5)] = (m2_el*pdsds)/DTOT_el;
	dAdU[conta_new(4,1,1,2,5,2,5,5)] = (m2_el*pdsm1)/DTOT_el;
	dAdU[conta_new(4,1,1,3,5,2,5,5)] = (pds + m2_el*pdsm2)/DTOT_el;
	dAdU[conta_new(4,1,1,4,5,2,5,5)] = (m2_el*pdset)/DTOT_el;
	 
	dAdU[conta_new(4,1,2,0,5,2,5,5)] = (m2_el*(DTOT_el*pdgm1 - pm1))/DTOT2;
	dAdU[conta_new(4,1,2,1,5,2,5,5)] = (m2_el*pdsm1)/DTOT_el;
	dAdU[conta_new(4,1,2,2,5,2,5,5)] = (m2_el*pm1m1)/DTOT_el;
	dAdU[conta_new(4,1,2,3,5,2,5,5)] = pm1/DTOT_el;
			
	dAdU[conta_new(4,1,3,0,5,2,5,5)] = -((etot_el - DTOT_el*(pdg + m2_el*pdgm2) + m2_el*pm2 + p_el)/DTOT2);
	dAdU[conta_new(4,1,3,1,5,2,5,5)] = (pds + m2_el*pdsm2)/DTOT_el;  
	dAdU[conta_new(4,1,3,2,5,2,5,5)] = pm1/DTOT_el;
	dAdU[conta_new(4,1,3,3,5,2,5,5)] = (2*pm2 + m2_el*pm2m2)/DTOT_el;
	dAdU[conta_new(4,1,3,4,5,2,5,5)] = (1 + pet)/DTOT_el;
			
	dAdU[conta_new(4,1,4,0,5,2,5,5)] = (m2_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	dAdU[conta_new(4,1,4,1,5,2,5,5)] = (m2_el*pdset)/DTOT_el;
	dAdU[conta_new(4,1,4,3,5,2,5,5)] = (1 + pet)/DTOT_el;

    for (i = 0; i < nScalarVariables*nScalarVariables; i++)   S(i) = 0.0;
																				// CONTROLLARE S					
    S(2*nScalarVariables + 0) = f_gauss(0);
    S(3*nScalarVariables + 0) = f_gauss(1);
    S(4*nScalarVariables + 0) = r_gauss;
    S(4*nScalarVariables + 2) = f_gauss(0);
    S(4*nScalarVariables + 3) = f_gauss(1);

    for (i = 0; i < nodesElement*nScalarVariables*nScalarVariables; i++)     Lstar[i] = 0.0;

    for (i = 0; i < nodesElement; i++){
			
		pp = i*nScalarVariables*nScalarVariables;
		
		Lstar[pp + 2*nScalarVariables + 0] = N[i]*S[2*nScalarVariables + 0];
		Lstar[pp + 3*nScalarVariables + 0] = N[i]*S[3*nScalarVariables + 0];
				
		for (k = 0; k < nScalarVariables - 1; k++){
			Lstar[pp + 4*nScalarVariables + k] = N[i]*S[4*nScalarVariables + k];
		}
    }

    for (k = 0; k < nScalarVariables; k++){
		for ( s = 0; s < nNodalVariables; s++){
			
			p = s*nScalarVariables + k;
			
			for (i = 0; i < nScalarVariables; i++){
				
				pp = i*nNodalVariables + s;

				for (j = 0; j < SpaceDimension; j++){

					t = conta_new(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

					Lstar[p] += (-A[t]*gradV[(SpaceDimension*i+j)*nNodalVariables + s]);
                    
	
				}
			}
		}
	}

    for (i = 0; i < nScalarVariables*nScalarVariables; i++) B[i] = 0.0;

    for (i = 0; i < nScalarVariables; i++){
        for (j = 0; j < SpaceDimension; j++){
            for (k = 0; k < nScalarVariables; k++){
                pp = i * nScalarVariables + k;
                for (m = 0; m < nScalarVariables; m++){
                    tt = conta_new(i,j,k,m,nScalarVariables,SpaceDimension,nScalarVariables,nScalarVariables);

                    B[pp]  += dAdU[tt]*gradU[m*SpaceDimension + j];

                }
            }
        }
	}

    for (j = 0; j < nodesElement; j++){
        pp = j*nScalarVariables*nScalarVariables;
        for (i = 0; i < nScalarVariables; i++){
            p = i*nScalarVariables;
            for (k = 0; k < nScalarVariables; k++){
                Lstar[pp + p + k] -= N[j]*B[p + k];
            }
        }
    }


	invtauStab[0] =	stab_c2*(norm_u + SpeedSound)/h; 
    invtauStab[1] =	stab_c2*(norm_u + SpeedSound)/h;
	invtauStab[2] =	stab_c1*mu_mixture/(DTOT_el*h*h) + invtauStab[0];
	invtauStab[3] =	invtauStab[2];
	invtauStab[4] = stab_c1*lambda_mixture/(DTOT_el*Cp_mixture*h*h) + invtauStab[0];  

	
// controllare L
    L[0] = 0.0;
    L[1] = 0.0;
    L[2] = -S[2*nScalarVariables + 0]*U_gauss[0] - S[2*nScalarVariables + 1]*U_gauss[1];
    L[3] = -S[3*nScalarVariables + 0]*U_gauss[0] - S[3*nScalarVariables + 1]*U_gauss[1];
    L[4] = 0.0;
    

    for (k = 0; k < nScalarVariables - 1 ; k++){
        L[4] -= S[4*nScalarVariables + k]*U_gauss[k];
    }

    for (i = 0; i < nScalarVariables; i++ ){
		for (k = 0; k < nScalarVariables; k++){
			for (j = 0; j < SpaceDimension; j++){

				s = conta_new(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

				L[i] += A[s]*gradU[k*SpaceDimension + j];
			}
		}

		Residual[i] = -L[i];

		for (k = 0; k < nodesElement; k++){
            Residual[i] -= N[k]*UUp(k,i);
        }
	}

	for (i = 0; i < nScalarVariables; i++){
        for (k = 0; k < nodesElement; k++){
			FConv[i + k*nScalarVariables] = N[k]*L[i];
		}
    }

    // Build diffusive term: stress tensor and thermal diffusion

    
    // HERE FiND TAU AND q

    for (i = 0; i < sizeK;  i++)    K[i]  = 0.0;
    for (i = 0; i < sizeKT; i++)    KT[i] = 0.0;
			
	K[conta_new(0,0,0,0,2,2,5,2)] =  (-2.0 + ctau)*m1_el/DTOT2;
	K[conta_new(0,0,0,1,2,2,5,2)] =  ctau*m2_el/DTOT2;
	K[conta_new(0,0,2,0,2,2,5,2)] =  (2.0 - ctau)/DTOT_el;
	K[conta_new(0,0,3,1,2,2,5,2)] = -ctau/DTOT_el;  
	
	K[conta_new(0,1,0,0,2,2,5,2)] = -m2_el/DTOT2;
	K[conta_new(0,1,0,1,2,2,5,2)] = -m1_el/DTOT2;
	K[conta_new(0,1,2,1,2,2,5,2)] =  1.0/DTOT_el; 
	K[conta_new(0,1,3,0,2,2,5,2)] =  1.0/DTOT_el;
	
////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	K[conta_new(1,0,0,0,2,2,5,2)] = -m2_el/DTOT2;
	K[conta_new(1,0,0,1,2,2,5,2)] = -m1_el/DTOT2;
	K[conta_new(1,0,2,1,2,2,5,2)] =  1.0/DTOT_el;
	K[conta_new(1,0,3,0,2,2,5,2)] =  1.0/DTOT_el; 
	
	K[conta_new(1,1,0,0,2,2,5,2)] =  ctau*m1_el/DTOT2;
	K[conta_new(1,1,0,1,2,2,5,2)] =  (ctau - 2.0)*m2_el/DTOT2;
	K[conta_new(1,1,2,0,2,2,5,2)] = -ctau/DTOT_el;
	K[conta_new(1,1,3,1,2,2,5,2)] =  (2.0 - ctau)/DTOT_el;


	
	KT[0] =  pdg/(DG_el*R) - p_el/(DG_el*DG_el*R);
	KT[1] =  pds/(DG_el*R) + p_el/(DG_el*DG_el*R);
	KT[2] =  pm1/(DG_el*R);
	KT[3] =  pm2/(DG_el*R);
	KT[4] =  pet/(DG_el*R);


    for (i = 0; i < SpaceDimension; i++){
        
        q(i) = 0.0;
		dt_diff(i) = 0.0;
		ds_diff(i) = 0.0;

        for (j = 0; j < SpaceDimension; j++){

            tau(i*SpaceDimension + j) = 0.0;
            
            for (l = 0; l < nScalarVariables; l++){
                for (m = 0; m < SpaceDimension; m++){

                    tau(i*SpaceDimension + j) += mu_mixture*K[conta_new(i,j,l,m,2,2,5,2)]*gradU[l*SpaceDimension + m];

                }
            }
        }
        for (l = 0; l < nScalarVariables; l++){
            q(i) -= lambda_mixture*KT(l)*gradU(l*SpaceDimension + i);
        }
    }

    ShockCapturing2d_new(mu_mixture, lambda_mixture, Cv_mixture, ros, h, dt_diff, ds_diff,tau, q, DTOT_el, gradU, Residual,U_gauss,a_el,norm_u,SpeedSound);

    // Build diffusive term: Diffusion tensor

	for ( i = 0; i < nScalarVariables*SpaceDimension; i++ )    G[i] = 0.0;
		
	for (i = 0; i < SpaceDimension; i++){

		for (j = 0; j < SpaceDimension; j++)
			G[(i + 2)*SpaceDimension + j] = -tau[i*SpaceDimension + j];
    }

	for (j = 0; j < SpaceDimension; j++){

		G[0*SpaceDimension + j] = -dt_diff[j];
		G[1*SpaceDimension + j] = -ds_diff[j];   // Numerical diffusivity for dust
		G[4*SpaceDimension + j] = q[j];

		for (i = 0; i < SpaceDimension; i++)
            G[4*SpaceDimension + j] += (-U_gauss[i + 2]/DTOT_el*tau[i*SpaceDimension + j]);

	}


    // Build diffusive term: Diffusion force

	for (s = 0; s < nNodalVariables; s++){

		FDiff[s] = 0.0;

		for (i = 0; i < nScalarVariables; i++){
			for ( j = 0; j < SpaceDimension; j++){
				FDiff[s] -= gradV[(SpaceDimension*i + j)*nNodalVariables + s]*G[i*SpaceDimension + j];
            }
        }
 	}

    // Stabilizing residual part

	for (s = 0; s < nNodalVariables; s++){

		FStab[s] = 0.0;

		for (k = 0; k < nScalarVariables; k++){
			FStab[s] += Lstar[s*nScalarVariables + k]*Residual[k]/invtauStab[k]*switchStab[k];
		}
	}

    // Force contribuution at the Gauss Point   

	int check = 1;

    for (i = 0; i < nNodalVariables; i++){

		F[i] = - sw_conv*FConv[i] - sw_diff*FDiff[i] - sw_stab*FStab[i];

		if (std::isnan(FConv[i]) == 1 || std::isnan(FDiff[i]) == 1 || std::isnan(FStab[i]) == 1){
			printf("%d %.3e %.3e %.3e\n", i, FConv[i], FDiff[i], FStab[i]);
			check = 0;
		}

	}
    
   
	if (check == 0)	{

        printf("%.3e %.3e %.3e %.3e %.3e \n\n", U_gauss(0), U_gauss(1), U_gauss(2), U_gauss(3), U_gauss(4));

		for (i = 0; i < nodesElement; i++){
            for (j = 0; j < nScalarVariables; j++){

			    printf("%.3e ",UU(i,j));
            }
            printf("\n");
		}
        printf("\n");
		printf("stab_c2 = %.3e - norm_u = %.3e - SpeedSound = %.3e - stab_c1 = %.3e - mu_mixture = %.3e - lambda_mixture = %.3e"
			"Cp_mixture = %.3e \n", stab_c2, norm_u, SpeedSound, stab_c1, mu_mixture, lambda_mixture, Cp_mixture);

		double sps2 = DG_el*R*Cpmixed*(2*etot_el - DTOT_el*norm2u)/(2*DTOT_el*Cmixed2);
		
		printf("spsound2 = %.3e\n",sps2);

		printf("DG_el = %.3e - Cpmixed = %.3e - en-u2 = %.3e - DTOT_el = %.3e - Cmixed2 = %.3e\n",
		DG_el, Cpmixed, 2*etot_el - DTOT_el*norm2u, DTOT_el, Cmixed2);
		
		printf("%.3e %.3e \n",2*etot_el, DTOT_el*norm2u);

		printf("cos cos cos \n\n");
		


		abort();
	}


    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector =  F*data.volume/NGaussPoints;        // Da controllare, dovrebbe essere Gauss Points  

    KRATOS_CATCH("")
}

template<>
void CompressibleNSBiphaseExplicit<3,4>::CalculateRightHandSideInternal(
    BoundedVector<double, 24> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int NGaussPoints = 1;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double h = data.h;
    const unsigned  int nodesElement = 4;   //OK
    const unsigned  int SpaceDimension = 3; //OK
    const unsigned  int nScalarVariables = SpaceDimension + 3;  //OK
    const unsigned  int nNodalVariables = nScalarVariables*nodesElement; //OK
    

    unsigned int i, j, k, l, m, s, t, tt, p, pp;
    
    const unsigned int size3 = nScalarVariables * SpaceDimension * nScalarVariables;
    const unsigned int size4 = nScalarVariables * SpaceDimension * nScalarVariables * nScalarVariables;
    const unsigned int sizeK = SpaceDimension * SpaceDimension * nScalarVariables * SpaceDimension;
    const unsigned int sizeKT = nScalarVariables;


    array_1d<double, size3>     A;
    array_1d<double, size4>     dAdU;
    array_1d<double, sizeK>     K;
    array_1d<double, sizeKT>    KT;

    array_1d<double,nNodalVariables>                     LumpedMassMatrix;
    array_1d<double,nScalarVariables>                    U_gauss;
    array_1d<double,nNodalVariables>                     U;
    array_1d<double,nNodalVariables>                     Un;
    array_1d<double,nNodalVariables>                     up;
    array_1d<double,nScalarVariables>                    L;
    array_1d<double,nScalarVariables>                    Residual;
    array_1d<double,nNodalVariables*nScalarVariables>    Lstar;
    array_1d<double,nScalarVariables*SpaceDimension>     G;
    array_1d<double,nScalarVariables*nScalarVariables>   S;
    array_1d<double,nScalarVariables*nScalarVariables>   B;
    array_1d<double,SpaceDimension>                      dt_diff;
	array_1d<double,SpaceDimension>                      ds_diff;
	array_1d<double,SpaceDimension*SpaceDimension>       tau;
    array_1d<double,SpaceDimension>                      q;
	array_1d<double,SpaceDimension>                      a_el;
    array_1d<double,nScalarVariables*nNodalVariables>    NN;
    array_1d<double,nScalarVariables*SpaceDimension>     gradU;
    array_1d<double,(nScalarVariables*SpaceDimension)*nNodalVariables>   gradV;
    array_1d<double,nScalarVariables>                    invtauStab;
	array_1d<double,nScalarVariables>                    switchStab;

    array_1d<double,nNodalVariables>     FConv;
    array_1d<double,nNodalVariables>     FStab;
    array_1d<double,nNodalVariables>     FDiff;
    array_1d<double,nNodalVariables>     F;
    
	
    
    const double& ctau = 2.0/3.0;   // This coefficient multiplies the divergence of the velocity in the calculation of tau. 
                                // In 3d would be 0.66667   //OK
    
    const double& dt = data.time_step;

    // In this implementation this function returns only nodal forces, rgardless of the time integration scheme used.
    // const double& bdf0 = data.bdf0;
    // const double& bdf1 = data.bdf1;
    // const double& bdf2 = data.bdf2;

    const BoundedMatrix<double,nodesElement,nScalarVariables>& UU = data.U;			
    //const BoundedMatrix<double,nodesElement,nScalarVariables>& UUn = data.Un;
    const BoundedMatrix<double,nodesElement,nScalarVariables>& Up = data.dUdt;    // Useful for the stabilizing part.
    BoundedMatrix<double,nodesElement,nScalarVariables> UUp = Up;

    const BoundedMatrix<double,nodesElement,SpaceDimension>& f_ext = data.f_ext;			
    const array_1d<double,nodesElement>& r = data.r_ext;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double Cv = data.c_v;
    const double gamma = data.gamma;
    const double Cp = Cv * gamma;
    const double R = Cp - Cv;
    const double ros = data.ros;
    const double Cs = data.c_s;

	const double sw_conv = 1.0;
    const double sw_diff = 1.0;
    const double sw_stab = 1.0;


    const double stab_c1 = 4;
    const double stab_c2 = 2;

	switchStab[0] = 1.0;
	switchStab[1] = 1.0;
	switchStab[2] = 1.0;
	switchStab[3] = 1.0;
	switchStab[4] = 1.0;
    switchStab[5] = 1.0;    //OK

    // Get shape function values
    array_1d<double,nodesElement> N;

    N[0] = 0.25;
    N[1] = 0.25;
    N[2] = 0.25;
    N[3] = 0.25;

    const BoundedMatrix<double,nodesElement,SpaceDimension>& DN = data.DN_DX;	

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,SpaceDimension> f_gauss = prod(trans(f_ext), N);      
    const double r_gauss = N(0)*r(0) + N(1)*r(1) + N(2)*r(2) + N(3)*r(3);   //OK
    

    // Define U and Udot
    for (i = 0; i < nodesElement; i++){
        for (j = 0; j < nScalarVariables; j++){
            
    //        UUp(i,j) = (UU(i,j) - UUn(i,j))/dt;   // This could be done better, since we have an estimation of the udot from the previous substep
            
        }
    }

    for (i = 0; i < nScalarVariables; i++){
        U_gauss(i) = 0;
        for (j = 0; j < nodesElement; j++){
            U_gauss(i) += N(j)*UU(j,i);

        }
    }

    for (i = 0; i < nScalarVariables*SpaceDimension; i++)       gradU(i) = 0.0;
// Compute the gradient of the variables /dro/dx,dro/dy, dm1/dx dm1/dy dm2/dx dm2/dy de/dx de/dy)
    for (i = 0; i < nScalarVariables; i++){         // Controllare
        for (k = 0; k < nodesElement; k++){

            gradU(i*SpaceDimension + 0) += DN(k,0)*UU(k,i);
            gradU(i*SpaceDimension + 1) += DN(k,1)*UU(k,i);
            gradU(i*SpaceDimension + 2) += DN(k,2)*UU(k,i);         //OK

        }
    }

    // for (i = 0; i < nNodalVariables*nScalarVariables; i++ )  NN(i) = 0.0;

    // This is convenient during the implementation but has to be removed 
    // after some modification of the remainding part of the file
   
    // This is convenient during the implementation but has to be removed 
    // after some modification of the remainding part of the file    
    for ( i = 0; i < nScalarVariables * SpaceDimension * nNodalVariables; i++)  gradV[i] = 0.0;

	for (i = 0; i < nodesElement; i++ ){
        
        gradV[ 0*nNodalVariables + nScalarVariables*i    ] = DN(i,0);
        gradV[ 1*nNodalVariables + nScalarVariables*i    ] = DN(i,1);
		gradV[ 2*nNodalVariables + nScalarVariables*i    ] = DN(i,2);

        gradV[ 3*nNodalVariables + nScalarVariables*i + 1] = DN(i,0);
        gradV[ 4*nNodalVariables + nScalarVariables*i + 1] = DN(i,1);
		gradV[ 5*nNodalVariables + nScalarVariables*i + 1] = DN(i,2);

        gradV[ 6*nNodalVariables + nScalarVariables*i + 2] = DN(i,0);
        gradV[ 7*nNodalVariables + nScalarVariables*i + 2] = DN(i,1);
		gradV[ 8*nNodalVariables + nScalarVariables*i + 2] = DN(i,2);

        gradV[ 9*nNodalVariables + nScalarVariables*i + 3] = DN(i,0);
        gradV[10*nNodalVariables + nScalarVariables*i + 3] = DN(i,1);
		gradV[11*nNodalVariables + nScalarVariables*i + 3] = DN(i,2);

        gradV[12*nNodalVariables + nScalarVariables*i + 4] = DN(i,0);  // Verify these lines
        gradV[13*nNodalVariables + nScalarVariables*i + 4] = DN(i,1);
		gradV[14*nNodalVariables + nScalarVariables*i + 4] = DN(i,2);

        gradV[15*nNodalVariables + nScalarVariables*i + 5] = DN(i,0);  // Verify these lines  //OK //but still verify
        gradV[16*nNodalVariables + nScalarVariables*i + 5] = DN(i,1);
		gradV[17*nNodalVariables + nScalarVariables*i + 5] = DN(i,2);

	}


    const double DTOT_el 	= U_gauss(0);
    const double DS_el 		= U_gauss(1);
    const double m1_el 		= U_gauss(2);
    const double m2_el 		= U_gauss(3);
    const double m3_el 		= U_gauss(4);
    const double etot_el 	= U_gauss(5);  //OK

    const double DG_el = DTOT_el - DS_el; 

    const double  kk = DS_el/DG_el;

    // HERE DEFINE MIXTURE COEFFICIENTS CV_M, ECC.

    const double mu_mixture = mu/(1 + kk);        // Controllare questo spaccimmo di termine (trovare su Neri et al.)

	const double Cp_mixture = (Cp + kk*Cs)/(1.0 + kk);
	const double Cv_mixture = (Cv + kk*Cs)/(1.0 + kk);

	const double Prandtl = 0.5;

//	const double lambda_mixture = Cp_mixture*mu_mixture/Prandtl;
	const double lambda_mixture = lambda;

	// Fine variazione

    const double u1_el = m1_el/DTOT_el;
	const double u2_el = m2_el/DTOT_el;
    const double u3_el = m3_el/DTOT_el; //OK

	const double norm2u = u1_el*u1_el + u2_el*u2_el + u3_el*u3_el;  //OK
	const double norm_u = sqrt(norm2u);
	const double norm2m = m1_el*m1_el + m2_el*m2_el + m3_el*m3_el;  //OK

	a_el(0) = u1_el/norm_u;
	a_el(1) = u2_el/norm_u;
    a_el(2) = u3_el/norm_u; // OK

	if (norm_u == 0){
		a_el(0) = 0.0;
		a_el(1) = 0.0;	 
        a_el(2) = 0.0;	 //OK
	}	


    const   double	DG_el2    	= DG_el*DG_el;
	const   double	DS_el2    	= DS_el*DS_el;
	const   double	DTOT2		= DTOT_el*DTOT_el;
	const   double  DTOT3		= DTOT_el*DTOT2;
	
	const   double  Cmixed      = (Cs*DS_el + Cv*(-DS_el + DTOT_el));
	const   double  Cmixed2     = Cmixed*Cmixed;
	const   double  Cmixed3     = Cmixed*Cmixed2;
    const   double	Cpmixed     = Cs*DS_el + (Cv + R)*DG_el;
	
	const   double  Cs2Ds3		= Cs*Cs*DS_el*DS_el*DS_el;
	const   double  Cv2Dg3		= Cv*Cv*DG_el*DG_el*DG_el;

	
	const   double 	a_d 	= 0;
		
	const   double  p_el   	= (DG_el*(2*DTOT_el*etot_el - norm2m)*R)/(2*DTOT_el*Cmixed);
		
	const   double  pdg 	=  ((Cv*DG_el2*norm2m + Cs*DS_el*(2*DTOT2*etot_el - DS_el*norm2m))*R)/(2*DTOT2*Cmixed2);
	const   double  pds 	=  (Cs*(-2*DTOT_el*etot_el + norm2m)*R)/(2*Cmixed2);
	const   double  pm1 	= -(DG_el*m1_el*R)/(DTOT_el*Cmixed);
	const   double  pm2 	= -(DG_el*m2_el*R)/(DTOT_el*Cmixed);
	const   double  pm3 	= -(DG_el*m3_el*R)/(DTOT_el*Cmixed);
	const   double  pet 	= DG_el*R/Cmixed;
	
	const   double  pdgdg 	=  ((Cs2Ds3 + Cv2Dg3)*norm2m + Cs*Cv*DS_el*(-2*DTOT3*etot_el + DS_el*(-2*DS_el + 3*DTOT_el)*norm2m))/(DTOT3*Cmixed3/R);
	const   double  pdgds 	=  (Cs*(-Cs*DS_el*etot_el + Cv*(DS_el + DTOT_el)*etot_el -Cv*norm2m)*R)/Cmixed3;
	const   double  pdgm1 	=  ((-Cs*DS_el2 - Cv*DG_el2)*m1_el*R)/(DTOT2*Cmixed2);
	const   double  pdgm2 	=  ((-Cs*DS_el2 - Cv*DG_el2)*m2_el*R)/(DTOT2*Cmixed2);
	const   double  pdgm3 	=  ((-Cs*DS_el2 - Cv*DG_el2)*m3_el*R)/(DTOT2*Cmixed2);
	const   double  pdget 	=  Cs*DS_el*R/Cmixed2;
	
	const   double  pdsds 	=  (Cs*(Cs - Cv)*(2*DTOT_el*etot_el - norm2m)*R)/Cmixed3;
	const   double  pdsm1 	=  Cs*m1_el*R/Cmixed2;
	const   double  pdsm2 	=  Cs*m2_el*R/Cmixed2;
	const   double  pdsm3 	=  Cs*m3_el*R/Cmixed2;
	const   double  pdset 	= -Cs*DTOT_el*R/Cmixed2;
	
	const   double  pm1m1 	= -DG_el*R/(DTOT_el*Cmixed);
	
	const   double  pm2m2 	= -DG_el*R/(DTOT_el*Cmixed);

	const   double  pm3m3 	= -DG_el*R/(DTOT_el*Cmixed);
	
	const double    Temperature = p_el/(DG_el*R);

	const double gas_concentration = 1 - DS_el/ros;
	const double gas_density = DG_el/gas_concentration;

	// printf("eps_g = %.3e\n", gas_concentration);

    double SpeedSound2  = DG_el*R*Cpmixed*(2*etot_el - DTOT_el*norm2u)/(2*DTOT_el*Cmixed2);  // verificare sound_speed
	double SpeedSound;

	SpeedSound = sqrt(SpeedSound2);
		 

	// printf("gamma = %.3e - SpeedSound = %.3e\n", gamma, SpeedSound); 

	if (dt > h/(SpeedSound + norm_u))   printf("dt = %.3e dt_crit = %.3e\n", dt, h/(SpeedSound + norm_u));

	if (DTOT_el < 0)    				printf("dtot_el = %.3e \n", DTOT_el);

    for (i = 0; i < size3; i++)     A(i) = 0.0;
	for (i = 0; i < size4; i++)     dAdU(i) = 0.0;
			
	// Build A

	A[conta_new(0,0,2,0,6,3,6,1)] = 1.0;
	A[conta_new(0,1,3,0,6,3,6,1)] = 1.0;
	A[conta_new(0,2,4,0,6,3,6,1)] = 1.0;
	
	A[conta_new(1,0,0,0,6,3,6,1)] = -DS_el*m1_el/DTOT2; 
	A[conta_new(1,0,1,0,6,3,6,1)] =  m1_el/DTOT_el;
	A[conta_new(1,0,2,0,6,3,6,1)] =  DS_el/DTOT_el;
	
	A[conta_new(1,1,0,0,6,3,6,1)] = -DS_el*m2_el/DTOT2; 
	A[conta_new(1,1,1,0,6,3,6,1)] =  m2_el/DTOT_el;
	A[conta_new(1,1,3,0,6,3,6,1)] =  DS_el/DTOT_el;
	
	A[conta_new(1,2,0,0,6,3,6,1)] = -DS_el*m3_el/DTOT2; 
	A[conta_new(1,2,1,0,6,3,6,1)] =  m3_el/DTOT_el;
	A[conta_new(1,2,4,0,6,3,6,1)] =  DS_el/DTOT_el;
	
	A[conta_new(2,0,0,0,6,3,6,1)] =  -m1_el*m1_el/DTOT2 + pdg;
	A[conta_new(2,0,1,0,6,3,6,1)] =  pds;
	A[conta_new(2,0,2,0,6,3,6,1)] =  2*m1_el/DTOT_el + pm1;
	A[conta_new(2,0,3,0,6,3,6,1)] =  pm2;
	A[conta_new(2,0,4,0,6,3,6,1)] =  pm3;
	A[conta_new(2,0,5,0,6,3,6,1)] =  pet;
	
	A[conta_new(2,1,0,0,6,3,6,1)] =  -m1_el*m2_el/DTOT2;
	A[conta_new(2,1,2,0,6,3,6,1)] =  m2_el/DTOT_el;
	A[conta_new(2,1,3,0,6,3,6,1)] =  m1_el/DTOT_el;
	
	A[conta_new(2,2,0,0,6,3,6,1)] =  -m1_el*m3_el/DTOT2;
	A[conta_new(2,2,2,0,6,3,6,1)] =  m3_el/DTOT_el;
	A[conta_new(2,2,4,0,6,3,6,1)] =  m1_el/DTOT_el;
	
	A[conta_new(3,0,0,0,6,3,6,1)] =  -m1_el*m2_el/DTOT2;
	A[conta_new(3,0,2,0,6,3,6,1)] =  m2_el/DTOT_el;
	A[conta_new(3,0,3,0,6,3,6,1)] =  m1_el/DTOT_el;
		
	A[conta_new(3,1,0,0,6,3,6,1)] =  -m2_el*m2_el/DTOT2 + pdg;
	A[conta_new(3,1,1,0,6,3,6,1)] =  pds;
	A[conta_new(3,1,2,0,6,3,6,1)] =  pm1;
	A[conta_new(3,1,3,0,6,3,6,1)] =  2*m2_el/DTOT_el + pm2;
	A[conta_new(3,1,4,0,6,3,6,1)] =  pm3;
	A[conta_new(3,1,5,0,6,3,6,1)] =  pet;
	
	A[conta_new(3,2,0,0,6,3,6,1)] =  -m2_el*m3_el/DTOT2;
	A[conta_new(3,2,3,0,6,3,6,1)] =  m3_el/DTOT_el;
	A[conta_new(3,2,4,0,6,3,6,1)] =  m2_el/DTOT_el;
	
	A[conta_new(4,0,0,0,6,3,6,1)] =  -m1_el*m3_el/DTOT2;
	A[conta_new(4,0,2,0,6,3,6,1)] =  m3_el/DTOT_el;
	A[conta_new(4,0,4,0,6,3,6,1)] =  m1_el/DTOT_el;
		
	A[conta_new(4,1,0,0,6,3,6,1)] =  -m2_el*m3_el/DTOT2;
	A[conta_new(4,1,3,0,6,3,6,1)] =  m3_el/DTOT_el;
	A[conta_new(4,1,4,0,6,3,6,1)] =  m2_el/DTOT_el;
	
	A[conta_new(4,2,0,0,6,3,6,1)] =  -m3_el*m3_el/DTOT2 + pdg;
	A[conta_new(4,2,1,0,6,3,6,1)] =  pds;
	A[conta_new(4,2,2,0,6,3,6,1)] =  pm1;
	A[conta_new(4,2,3,0,6,3,6,1)] =  pm2;
	A[conta_new(4,2,4,0,6,3,6,1)] =  2*m3_el/DTOT_el + pm3;
	A[conta_new(4,2,5,0,6,3,6,1)] =  pet;
	
	A[conta_new(5,0,0,0,6,3,6,1)] =  -(m1_el*(etot_el - DTOT_el*pdg + p_el))/DTOT2;
	A[conta_new(5,0,1,0,6,3,6,1)] =  m1_el*pds/DTOT_el;
	A[conta_new(5,0,2,0,6,3,6,1)] =  (etot_el + m1_el*pm1 + p_el)/DTOT_el;
	A[conta_new(5,0,3,0,6,3,6,1)] =  m1_el*pm2/DTOT_el;
	A[conta_new(5,0,4,0,6,3,6,1)] =  m1_el*pm3/DTOT_el;
	A[conta_new(5,0,5,0,6,3,6,1)] =	 m1_el*(1.0 + pet)/DTOT_el;
	
	A[conta_new(5,1,0,0,6,3,6,1)] =  -(m2_el*(etot_el - DTOT_el*pdg + p_el))/DTOT2;
	A[conta_new(5,1,1,0,6,3,6,1)] =  m2_el*pds/DTOT_el;
	A[conta_new(5,1,2,0,6,3,6,1)] =	 m2_el*pm1/DTOT_el;  
	A[conta_new(5,1,3,0,6,3,6,1)] =  (etot_el + m2_el*pm2 + p_el)/DTOT_el;
	A[conta_new(5,1,4,0,6,3,6,1)] =	 m2_el*pm3/DTOT_el;
	A[conta_new(5,1,5,0,6,3,6,1)] =  m2_el*(1.0 + pet)/DTOT_el;

	A[conta_new(5,2,0,0,6,3,6,1)] =  -(m3_el*(etot_el - DTOT_el*pdg + p_el))/DTOT2;
	A[conta_new(5,2,1,0,6,3,6,1)] =  m3_el*pds/DTOT_el;
	A[conta_new(5,2,2,0,6,3,6,1)] =	m3_el*pm1/DTOT_el;  
	A[conta_new(5,2,3,0,6,3,6,1)] =  m3_el*pm2/DTOT_el;
	A[conta_new(5,2,4,0,6,3,6,1)] =	(etot_el + m3_el*pm3 + p_el)/DTOT_el;
	A[conta_new(5,2,5,0,6,3,6,1)] =  m3_el*(1.0 + pet)/DTOT_el;

    //	Build dAdU

// 1  ////////////////////////////////////////////////////////
			
	dAdU[conta_new(1,0,0,0,6,3,6,6)] = 2*DS_el*m1_el/DTOT3; 
	dAdU[conta_new(1,0,0,1,6,3,6,6)] = -m1_el/DTOT2;
	dAdU[conta_new(1,0,0,2,6,3,6,6)] = -DS_el/DTOT2;
	
	dAdU[conta_new(1,0,1,0,6,3,6,6)] = -m1_el/DTOT2;  
	dAdU[conta_new(1,0,1,2,6,3,6,6)] = 1.0/DTOT_el;
	 
	dAdU[conta_new(1,0,2,0,6,3,6,6)] = -DS_el/DTOT2; 
	dAdU[conta_new(1,0,2,1,6,3,6,6)] =  1.0/DTOT_el;
	
/////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(1,1,0,0,6,3,6,6)] =  2*DS_el*m2_el/DTOT3;
	dAdU[conta_new(1,1,0,1,6,3,6,6)] =  -m2_el/DTOT2;
	dAdU[conta_new(1,1,0,3,6,3,6,6)] = -DS_el/DTOT2; 
	
	dAdU[conta_new(1,1,1,0,6,3,6,6)] =  -m2_el/DTOT2;
	dAdU[conta_new(1,1,1,3,6,3,6,6)] = 1.0/DTOT_el; 
	 
	dAdU[conta_new(1,1,3,0,6,3,6,6)] = -DS_el/DTOT2;
	dAdU[conta_new(1,1,3,1,6,3,6,6)] =  1.0/DTOT_el;

/////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(1,2,0,0,6,3,6,6)] =  2*DS_el*m3_el/DTOT3;
	dAdU[conta_new(1,2,0,1,6,3,6,6)] =  -m3_el/DTOT2;
	dAdU[conta_new(1,2,0,4,6,3,6,6)] = -DS_el/DTOT2; 
	
	dAdU[conta_new(1,2,1,0,6,3,6,6)] =  -m3_el/DTOT2;
	dAdU[conta_new(1,2,1,4,6,3,6,6)] = 1.0/DTOT_el; 
	 
	dAdU[conta_new(1,2,4,0,6,3,6,6)] = -DS_el/DTOT2;
	dAdU[conta_new(1,2,4,1,6,3,6,6)] =  1.0/DTOT_el;
	
// 2  /////////////////////////////////////////////////////////////
			
	dAdU[conta_new(2,0,0,0,6,3,6,6)] =  2*m1_el*m1_el/DTOT3 + pdgdg;	 
	dAdU[conta_new(2,0,0,1,6,3,6,6)] =  pdgds;
	dAdU[conta_new(2,0,0,2,6,3,6,6)] = -2*m1_el/DTOT2 + pdgm1;
	dAdU[conta_new(2,0,0,3,6,3,6,6)] =  pdgm2;
	dAdU[conta_new(2,0,0,4,6,3,6,6)] =  pdgm3;
	dAdU[conta_new(2,0,0,5,6,3,6,6)] =  pdget;
	
	dAdU[conta_new(2,0,1,0,6,3,6,6)] = pdgds;
	dAdU[conta_new(2,0,1,1,6,3,6,6)] = pdsds;
	dAdU[conta_new(2,0,1,2,6,3,6,6)] = pdsm1;
	dAdU[conta_new(2,0,1,3,6,3,6,6)] = pdsm2;
	dAdU[conta_new(2,0,1,4,6,3,6,6)] = pdsm3;
	dAdU[conta_new(2,0,1,5,6,3,6,6)] = pdset;
	 
	dAdU[conta_new(2,0,2,0,6,3,6,6)] = -2*m1_el/DTOT2 + pdgm1;
	dAdU[conta_new(2,0,2,1,6,3,6,6)] = pdsm1;
	dAdU[conta_new(2,0,2,2,6,3,6,6)] = 2.0/DTOT_el + pm1m1;
	 
	dAdU[conta_new(2,0,3,0,6,3,6,6)] = pdgm2; 
	dAdU[conta_new(2,0,3,1,6,3,6,6)] = pdsm2;
	dAdU[conta_new(2,0,3,3,6,3,6,6)] = pm2m2; 
	
	dAdU[conta_new(2,0,4,0,6,3,6,6)] = pdgm3;
	dAdU[conta_new(2,0,4,1,6,3,6,6)] = pdsm3;
	dAdU[conta_new(2,0,4,4,6,3,6,6)] = pm3m3;
		
	dAdU[conta_new(2,0,5,0,6,3,6,6)] = pdget;
	dAdU[conta_new(2,0,5,1,6,3,6,6)] = pdset;
	 
	dAdU[conta_new(2,1,0,0,6,3,6,6)] = 2*m1_el*m2_el/DTOT3;	 
	dAdU[conta_new(2,1,0,2,6,3,6,6)] = -m2_el/DTOT2;
	dAdU[conta_new(2,1,0,3,6,3,6,6)] = -m1_el/DTOT2;
	
	dAdU[conta_new(2,1,2,0,6,3,6,6)] = -m2_el/DTOT2; 
	dAdU[conta_new(2,1,2,3,6,3,6,6)] = 1.0/DTOT_el; 
			
	dAdU[conta_new(2,1,3,0,6,3,6,6)] = -m1_el/DTOT2; 
	dAdU[conta_new(2,1,3,2,6,3,6,6)] = 1.0/DTOT_el;
	
	dAdU[conta_new(2,2,0,0,6,3,6,6)] = 2*m1_el*m3_el/DTOT3;	 
	dAdU[conta_new(2,2,0,2,6,3,6,6)] = -m3_el/DTOT2;
	dAdU[conta_new(2,2,0,4,6,3,6,6)] = -m1_el/DTOT2;
	
	dAdU[conta_new(2,2,2,0,6,3,6,6)] = -m3_el/DTOT2; 
	dAdU[conta_new(2,2,2,4,6,3,6,6)] = 1.0/DTOT_el; 
			
	dAdU[conta_new(2,2,4,0,6,3,6,6)] = -m1_el/DTOT2; 
	dAdU[conta_new(2,2,4,2,6,3,6,6)] = 1.0/DTOT_el;
	

// 3 /////////////////////////////////////////////////////////////////////////////
			
	dAdU[conta_new(3,0,0,0,6,3,6,6)] = 2*m1_el*m2_el/DTOT3;	 
	dAdU[conta_new(3,0,0,2,6,3,6,6)] = -m2_el/DTOT2;
	dAdU[conta_new(3,0,0,3,6,3,6,6)] = -m1_el/DTOT2;
		
	dAdU[conta_new(3,0,2,0,6,3,6,6)] = -m2_el/DTOT2; 
	dAdU[conta_new(3,0,2,3,6,3,6,6)] = 1.0/DTOT_el; 
			
	dAdU[conta_new(3,0,3,0,6,3,6,6)] = -m1_el/DTOT2; 
	dAdU[conta_new(3,0,3,2,6,3,6,6)] = 1.0/DTOT_el;
	
//////////////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(3,1,0,0,6,3,6,6)] = 2*m2_el*m2_el/DTOT3 + pdgdg;	 
	dAdU[conta_new(3,1,0,1,6,3,6,6)] = pdgds;
	dAdU[conta_new(3,1,0,2,6,3,6,6)] = pdgm1;
	dAdU[conta_new(3,1,0,3,6,3,6,6)] = -2*m2_el/DTOT2 + pdgm2;
	dAdU[conta_new(3,1,0,4,6,3,6,6)] = pdgm3;
	dAdU[conta_new(3,1,0,5,6,3,6,6)] = pdget;
	
	dAdU[conta_new(3,1,1,0,6,3,6,6)] = pdgds;
	dAdU[conta_new(3,1,1,1,6,3,6,6)] = pdsds;
	dAdU[conta_new(3,1,1,2,6,3,6,6)] = pdsm1;
	dAdU[conta_new(3,1,1,3,6,3,6,6)] = pdsm2;
	dAdU[conta_new(3,1,1,4,6,3,6,6)] = pdsm3;
	dAdU[conta_new(3,1,1,5,6,3,6,6)] = pdset;
	 
	dAdU[conta_new(3,1,2,0,6,3,6,6)] = pdgm1;
	dAdU[conta_new(3,1,2,1,6,3,6,6)] = pdsm1;
	dAdU[conta_new(3,1,2,2,6,3,6,6)] = pm1m1;
	
	dAdU[conta_new(3,1,3,0,6,3,6,6)] = -2*m2_el/DTOT2 + pdgm2;
	dAdU[conta_new(3,1,3,1,6,3,6,6)] = pdsm2;
	dAdU[conta_new(3,1,3,3,6,3,6,6)] = 2.0/DTOT_el + pm2m2; 
	
	dAdU[conta_new(3,1,4,0,6,3,6,6)] = pdgm3;
	dAdU[conta_new(3,1,4,1,6,3,6,6)] = pdsm3;
	dAdU[conta_new(3,1,4,4,6,3,6,6)] = pm3m3;
	
	dAdU[conta_new(3,1,5,0,6,3,6,6)] = pdget;
	dAdU[conta_new(3,1,5,1,6,3,6,6)] = pdset;

//////////////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(3,2,0,0,6,3,6,6)] = 2*m2_el*m3_el/DTOT3;	 
	dAdU[conta_new(3,2,0,3,6,3,6,6)] = -m3_el/DTOT2;
	dAdU[conta_new(3,2,0,4,6,3,6,6)] = -m2_el/DTOT2;
	
	dAdU[conta_new(3,2,3,0,6,3,6,6)] = -m3_el/DTOT2;
	dAdU[conta_new(3,2,3,4,6,3,6,6)] =  1.0/DTOT_el; 
	
	dAdU[conta_new(3,2,4,0,6,3,6,6)] = -m2_el/DTOT2;
	dAdU[conta_new(3,2,4,3,6,3,6,6)] = 1.0/DTOT_el;

// 4 /////////////////////////////////////////////////////////////////////////////
			
	dAdU[conta_new(4,0,0,0,6,3,6,6)] = 2*m1_el*m3_el/DTOT3;	 
	dAdU[conta_new(4,0,0,2,6,3,6,6)] = -m3_el/DTOT2;
	dAdU[conta_new(4,0,0,4,6,3,6,6)] = -m1_el/DTOT2;
		
	dAdU[conta_new(4,0,2,0,6,3,6,6)] = -m3_el/DTOT2; 
	dAdU[conta_new(4,0,2,4,6,3,6,6)] = 1.0/DTOT_el; 
			
	dAdU[conta_new(4,0,4,0,6,3,6,6)] = -m1_el/DTOT2; 
	dAdU[conta_new(4,0,4,2,6,3,6,6)] = 1.0/DTOT_el;
	
//////////////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(4,1,0,0,6,3,6,6)] = 2*m2_el*m3_el/DTOT3;	 
	dAdU[conta_new(4,1,0,3,6,3,6,6)] = -m3_el/DTOT2;
	dAdU[conta_new(4,1,0,4,6,3,6,6)] = -m2_el/DTOT2;
	
	dAdU[conta_new(4,1,3,0,6,3,6,6)] = -m3_el/DTOT2;
	dAdU[conta_new(4,1,3,4,6,3,6,6)] = 1.0/DTOT_el;
	
	dAdU[conta_new(4,1,4,0,6,3,6,6)] = -m2_el/DTOT2;
	dAdU[conta_new(4,1,4,3,6,3,6,6)] = 1.0/DTOT_el;
	
//////////////////////////////////////////////////////////////////////// 
	 
	dAdU[conta_new(4,2,0,0,6,3,6,6)] = 2*m3_el*m3_el/DTOT3 + pdgdg;
	dAdU[conta_new(4,2,0,1,6,3,6,6)] = pdgds;
	dAdU[conta_new(4,2,0,2,6,3,6,6)] = pdgm1;
	dAdU[conta_new(4,2,0,3,6,3,6,6)] = pdgm2;
	dAdU[conta_new(4,2,0,4,6,3,6,6)] = -2*m3_el/DTOT2 + pdgm3;
	dAdU[conta_new(4,2,0,5,6,3,6,6)] = pdget;
	
	dAdU[conta_new(4,2,1,0,6,3,6,6)] = pdgds;
	dAdU[conta_new(4,2,1,1,6,3,6,6)] = pdsds;
	dAdU[conta_new(4,2,1,2,6,3,6,6)] = pdsm1;
	dAdU[conta_new(4,2,1,3,6,3,6,6)] = pdsm2;
	dAdU[conta_new(4,2,1,4,6,3,6,6)] = pdsm3;
	dAdU[conta_new(4,2,1,5,6,3,6,6)] = pdset;
	
	dAdU[conta_new(4,2,2,0,6,3,6,6)] = pdgm1;
	dAdU[conta_new(4,2,2,1,6,3,6,6)] = pdsm1;
	dAdU[conta_new(4,2,2,2,6,3,6,6)] = pm1m1;
	
	dAdU[conta_new(4,2,3,0,6,3,6,6)] = pdgm2;
	dAdU[conta_new(4,2,3,1,6,3,6,6)] = pdsm2;
	dAdU[conta_new(4,2,3,3,6,3,6,6)] = pm2m2; 
	
	dAdU[conta_new(4,2,4,0,6,3,6,6)] = -2*m3_el/DTOT2 + pdgm3;
	dAdU[conta_new(4,2,4,1,6,3,6,6)] = pdsm3;
	dAdU[conta_new(4,2,4,4,6,3,6,6)] = 2.0/DTOT_el + pm3m3;
	
	dAdU[conta_new(4,2,5,0,6,3,6,6)] = pdget;
	dAdU[conta_new(4,2,5,1,6,3,6,6)] = pdset;
	
// 5 //////////////////////////////////////////////////////////////////////
			
	dAdU[conta_new(5,0,0,0,6,3,6,6)] = (m1_el*(2*etot_el - 2*DTOT_el*pdg + DTOT2*pdgdg + 2*p_el))/DTOT3;
	dAdU[conta_new(5,0,0,1,6,3,6,6)] = (m1_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(5,0,0,2,6,3,6,6)] = -((etot_el - DTOT_el*(pdg + m1_el*pdgm1) + m1_el*pm1 + p_el)/DTOT2);
	dAdU[conta_new(5,0,0,3,6,3,6,6)] = (m1_el*(DTOT_el*pdgm2 - pm2))/DTOT2;
	dAdU[conta_new(5,0,0,4,6,3,6,6)] = (m1_el*(DTOT_el*pdgm3 - pm3))/DTOT2;
	dAdU[conta_new(5,0,0,5,6,3,6,6)] = (m1_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	
	dAdU[conta_new(5,0,1,0,6,3,6,6)] = (m1_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(5,0,1,1,6,3,6,6)] = (m1_el*pdsds)/DTOT_el;
	dAdU[conta_new(5,0,1,2,6,3,6,6)] = (pds + m1_el*pdsm1)/DTOT_el;
	dAdU[conta_new(5,0,1,3,6,3,6,6)] = (m1_el*pdsm2)/DTOT_el;
	dAdU[conta_new(5,0,1,4,6,3,6,6)] = (m1_el*pdsm3)/DTOT_el;
	dAdU[conta_new(5,0,1,5,6,3,6,6)] = (m1_el*pdset)/DTOT_el;
	 
	dAdU[conta_new(5,0,2,0,6,3,6,6)] = -((etot_el - DTOT_el*(pdg + m1_el*pdgm1) + m1_el*pm1 + p_el)/DTOT2);
	dAdU[conta_new(5,0,2,1,6,3,6,6)] = (pds + m1_el*pdsm1)/DTOT_el;
	dAdU[conta_new(5,0,2,2,6,3,6,6)] = (2*pm1 + m1_el*pm1m1)/DTOT_el;
	dAdU[conta_new(5,0,2,3,6,3,6,6)] = pm2/DTOT_el;
	dAdU[conta_new(5,0,2,4,6,3,6,6)] = pm3/DTOT_el;
	dAdU[conta_new(5,0,2,5,6,3,6,6)] = (1 + pet)/DTOT_el;
			
	dAdU[conta_new(5,0,3,0,6,3,6,6)] = (m1_el*(DTOT_el*pdgm2 - pm2))/DTOT2;
	dAdU[conta_new(5,0,3,1,6,3,6,6)] = (m1_el*pdsm2)/DTOT_el;
	dAdU[conta_new(5,0,3,2,6,3,6,6)] = pm2/DTOT_el;
	dAdU[conta_new(5,0,3,3,6,3,6,6)] = m1_el*pm2m2/DTOT_el;
	
	dAdU[conta_new(5,0,4,0,6,3,6,6)] = (m1_el*(DTOT_el*pdgm3 - pm3))/DTOT2;
	dAdU[conta_new(5,0,4,1,6,3,6,6)] = (m1_el*pdsm3)/DTOT_el;
	dAdU[conta_new(5,0,4,2,6,3,6,6)] = pm3/DTOT_el;
	dAdU[conta_new(5,0,4,4,6,3,6,6)] = m1_el*pm3m3/DTOT_el;
	
	dAdU[conta_new(5,0,5,0,6,3,6,6)] = (m1_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	dAdU[conta_new(5,0,5,1,6,3,6,6)] = (m1_el*pdset)/DTOT_el;
	dAdU[conta_new(5,0,5,2,6,3,6,6)] = (1 + pet)/DTOT_el;
	
	dAdU[conta_new(5,1,0,0,6,3,6,6)] = (m2_el*(2*etot_el - 2*DTOT_el*pdg + DTOT2*pdgdg + 2*p_el))/DTOT3;	 
	dAdU[conta_new(5,1,0,1,6,3,6,6)] = (m2_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(5,1,0,2,6,3,6,6)] = (m2_el*(DTOT_el*pdgm1 - pm1))/DTOT2; 
	dAdU[conta_new(5,1,0,3,6,3,6,6)] = -((etot_el - DTOT_el*(pdg + m2_el*pdgm2) + m2_el*pm2 + p_el)/DTOT2);
	dAdU[conta_new(5,1,0,4,6,3,6,6)] = (m2_el*(DTOT_el*pdgm3 - pm3))/DTOT2; 
	dAdU[conta_new(5,1,0,5,6,3,6,6)] = (m2_el*(-1 + DTOT_el*pdget - pet))/DTOT2;

	dAdU[conta_new(5,1,1,0,6,3,6,6)] = (m2_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(5,1,1,1,6,3,6,6)] = (m2_el*pdsds)/DTOT_el;
	dAdU[conta_new(5,1,1,2,6,3,6,6)] = (m2_el*pdsm1)/DTOT_el;
	dAdU[conta_new(5,1,1,3,6,3,6,6)] = (pds + m2_el*pdsm2)/DTOT_el;
	dAdU[conta_new(5,1,1,4,6,3,6,6)] = (m2_el*pdsm3)/DTOT_el;
	dAdU[conta_new(5,1,1,5,6,3,6,6)] = (m2_el*pdset)/DTOT_el;
	 
	dAdU[conta_new(5,1,2,0,6,3,6,6)] = (m2_el*(DTOT_el*pdgm1 - pm1))/DTOT2;
	dAdU[conta_new(5,1,2,1,6,3,6,6)] = (m2_el*pdsm1)/DTOT_el;
	dAdU[conta_new(5,1,2,2,6,3,6,6)] = (m2_el*pm1m1)/DTOT_el;
	dAdU[conta_new(5,1,2,3,6,3,6,6)] = pm1/DTOT_el;
			
	dAdU[conta_new(5,1,3,0,6,3,6,6)] = -((etot_el - DTOT_el*(pdg + m2_el*pdgm2) + m2_el*pm2 + p_el)/DTOT2);
	dAdU[conta_new(5,1,3,1,6,3,6,6)] = (pds + m2_el*pdsm2)/DTOT_el;  
	dAdU[conta_new(5,1,3,2,6,3,6,6)] = pm1/DTOT_el;
	dAdU[conta_new(5,1,3,3,6,3,6,6)] = (2*pm2 + m2_el*pm2m2)/DTOT_el;
	dAdU[conta_new(5,1,3,4,6,3,6,6)] = pm3/DTOT_el;
	dAdU[conta_new(5,1,3,5,6,3,6,6)] = (1 + pet)/DTOT_el;
	
	dAdU[conta_new(5,1,4,0,6,3,6,6)] = (m2_el*(DTOT_el*pdgm3 - pm3))/DTOT2;
	dAdU[conta_new(5,1,4,1,6,3,6,6)] = (m2_el*pdsm3)/DTOT_el;
	dAdU[conta_new(5,1,4,3,6,3,6,6)] = pm3/DTOT_el;
	dAdU[conta_new(5,1,4,4,6,3,6,6)] = (m2_el*pm3m3)/DTOT_el;
	
	dAdU[conta_new(5,1,5,0,6,3,6,6)] = (m2_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	dAdU[conta_new(5,1,5,1,6,3,6,6)] = (m2_el*pdset)/DTOT_el;
	dAdU[conta_new(5,1,5,3,6,3,6,6)] = (1 + pet)/DTOT_el;
	
	dAdU[conta_new(5,2,0,0,6,3,6,6)] = (m3_el*(2*etot_el - 2*DTOT_el*pdg + DTOT2*pdgdg + 2*p_el))/DTOT3;	 
	dAdU[conta_new(5,2,0,1,6,3,6,6)] = (m3_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(5,2,0,2,6,3,6,6)] = (m3_el*(DTOT_el*pdgm1 - pm1))/DTOT2; 
	dAdU[conta_new(5,2,0,3,6,3,6,6)] = (m3_el*(DTOT_el*pdgm2 - pm2))/DTOT2;
	dAdU[conta_new(5,2,0,4,6,3,6,6)] = -((etot_el - DTOT_el*(pdg + m3_el*pdgm3) + m3_el*pm3 + p_el)/DTOT2); 
	dAdU[conta_new(5,2,0,5,6,3,6,6)] = (m3_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	
	dAdU[conta_new(5,2,1,0,6,3,6,6)] = (m3_el*(DTOT_el*pdgds - pds))/DTOT2;
	dAdU[conta_new(5,2,1,1,6,3,6,6)] = (m3_el*pdsds)/DTOT_el;
	dAdU[conta_new(5,2,1,2,6,3,6,6)] = (m3_el*pdsm1)/DTOT_el;
	dAdU[conta_new(5,2,1,3,6,3,6,6)] = (m3_el*pdsm2)/DTOT_el;
	dAdU[conta_new(5,2,1,4,6,3,6,6)] = (pds + m3_el*pdsm3)/DTOT_el;
	dAdU[conta_new(5,2,1,5,6,3,6,6)] = (m3_el*pdset)/DTOT_el;
	
	dAdU[conta_new(5,2,2,0,6,3,6,6)] = (m3_el*(DTOT_el*pdgm1 - pm1))/DTOT2;
	dAdU[conta_new(5,2,2,1,6,3,6,6)] = (m3_el*pdsm1)/DTOT_el;
	dAdU[conta_new(5,2,2,2,6,3,6,6)] = (m3_el*pm1m1)/DTOT_el;
	dAdU[conta_new(5,2,2,4,6,3,6,6)] = pm1/DTOT_el;
	
	dAdU[conta_new(5,2,3,0,6,3,6,6)] = (m3_el*(DTOT_el*pdgm2 - pm2))/DTOT2;
	dAdU[conta_new(5,2,3,1,6,3,6,6)] = (m3_el*pdsm2)/DTOT_el;
	dAdU[conta_new(5,2,3,3,6,3,6,6)] = m3_el*pm2m2/DTOT_el;
	dAdU[conta_new(5,2,3,4,6,3,6,6)] = pm2/DTOT_el;
	
	dAdU[conta_new(5,2,4,0,6,3,6,6)] = -((etot_el - DTOT_el*(pdg + m3_el*pdgm3) + m3_el*pm3 + p_el)/DTOT2);
	dAdU[conta_new(5,2,4,1,6,3,6,6)] = (pds + m3_el*pdsm3)/DTOT_el;  
	dAdU[conta_new(5,2,4,2,6,3,6,6)] = pm1/DTOT_el;
	dAdU[conta_new(5,2,4,3,6,3,6,6)] = pm2/DTOT_el;
	dAdU[conta_new(5,2,4,4,6,3,6,6)] = (2*pm3 + m3_el*pm3m3)/DTOT_el;
	dAdU[conta_new(5,2,4,5,6,3,6,6)] = (1 + pet)/DTOT_el;
	
	dAdU[conta_new(5,2,5,0,6,3,6,6)] = (m3_el*(-1 + DTOT_el*pdget - pet))/DTOT2;
	dAdU[conta_new(5,2,5,1,6,3,6,6)] = (m3_el*pdset)/DTOT_el;
	dAdU[conta_new(5,2,5,4,6,3,6,6)] = (1 + pet)/DTOT_el;
/*
	for (i = 0; i < nScalarVariables; i++){
		for (j = 0; j < SpaceDimension; j++){
			for (l = 0; l < nScalarVariables; l++){
				for (m = 0; m < nScalarVariables; m++){
					if (dAdU[conta_new(i,j,l,m,6,3,6,6)] != dAdU[conta_new(i,j,m,l,6,3,6,6)])
						printf("Alert!!  Alert!!  Alert!!  \n\n\n");
				}
			}
		}
	}
*/


    for (i = 0; i < nScalarVariables*nScalarVariables; i++)   S(i) = 0.0;
																				// CONTROLLARE S					
    S(2*nScalarVariables + 0) = f_gauss(0);
    S(3*nScalarVariables + 0) = f_gauss(1);
    S(4*nScalarVariables + 0) = f_gauss(2); //OK
    S(5*nScalarVariables + 0) = r_gauss;
    S(5*nScalarVariables + 2) = f_gauss(0);
    S(5*nScalarVariables + 3) = f_gauss(1);
    S(5*nScalarVariables + 4) = f_gauss(2); //OK

    for (i = 0; i < nodesElement*nScalarVariables*nScalarVariables; i++)     Lstar[i] = 0.0;

    for (i = 0; i < nodesElement; i++){
			
		pp = i*nScalarVariables*nScalarVariables;
		
		Lstar[pp + 2*nScalarVariables + 0] = N[i]*S[2*nScalarVariables + 0];
		Lstar[pp + 3*nScalarVariables + 0] = N[i]*S[3*nScalarVariables + 0];
        Lstar[pp + 4*nScalarVariables + 0] = N[i]*S[4*nScalarVariables + 0];   // OK
				
		for (k = 0; k < nScalarVariables - 1; k++){
			Lstar[pp + 5*nScalarVariables + k] = N[i]*S[5*nScalarVariables + k];    //OK
		}
    }

    for (k = 0; k < nScalarVariables; k++){
		for ( s = 0; s < nNodalVariables; s++){
			
			p = s*nScalarVariables + k;
			
			for (i = 0; i < nScalarVariables; i++){
				
				pp = i*nNodalVariables + s;

				for (j = 0; j < SpaceDimension; j++){

					t = conta_new(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

					Lstar[p] += (-A[t]*gradV[(SpaceDimension*i+j)*nNodalVariables + s]);
                    
	
				}
			}
		}
	}

    for (i = 0; i < nScalarVariables*nScalarVariables; i++) B[i] = 0.0;

    for (i = 0; i < nScalarVariables; i++){
        for (j = 0; j < SpaceDimension; j++){
            for (k = 0; k < nScalarVariables; k++){
                pp = i * nScalarVariables + k;
                for (m = 0; m < nScalarVariables; m++){
                    tt = conta_new(i,j,k,m,nScalarVariables,SpaceDimension,nScalarVariables,nScalarVariables);

                    B[pp]  += dAdU[tt]*gradU[m*SpaceDimension + j];

                }
            }
        }
	}

    for (j = 0; j < nodesElement; j++){
        pp = j*nScalarVariables*nScalarVariables;
        for (i = 0; i < nScalarVariables; i++){
            p = i*nScalarVariables;
            for (k = 0; k < nScalarVariables; k++){
                Lstar[pp + p + k] -= N[j]*B[p + k];
            }
        }
    }


	invtauStab[0] =	stab_c2*(norm_u + SpeedSound)/h; 
    invtauStab[1] =	stab_c2*(norm_u + SpeedSound)/h;
	invtauStab[2] =	stab_c1*mu_mixture/(DTOT_el*h*h) + invtauStab[0];
	invtauStab[3] =	invtauStab[2];
    invtauStab[4] =	invtauStab[2];      // OK
	invtauStab[5] = stab_c1*lambda_mixture/(DTOT_el*Cp_mixture*h*h) + invtauStab[0];  

	
// controllare L
    L[0] = 0.0;
    L[1] = 0.0;
    L[2] = -S[2*nScalarVariables + 0]*U_gauss[0] - S[2*nScalarVariables + 1]*U_gauss[1];
    L[3] = -S[3*nScalarVariables + 0]*U_gauss[0] - S[3*nScalarVariables + 1]*U_gauss[1];
    L[4] = -S[4*nScalarVariables + 0]*U_gauss[0] - S[4*nScalarVariables + 1]*U_gauss[1];  // OK
    L[5] = 0.0;
    

    for (k = 0; k < nScalarVariables - 1 ; k++){
        L[5] -= S[5*nScalarVariables + k]*U_gauss[k];       // OK
    }

    for (i = 0; i < nScalarVariables; i++ ){
		for (k = 0; k < nScalarVariables; k++){
			for (j = 0; j < SpaceDimension; j++){

				s = conta_new(i,j,k,0,nScalarVariables,SpaceDimension,nScalarVariables,1);

				L[i] += A[s]*gradU[k*SpaceDimension + j];
			}
		}

		Residual[i] = -L[i];

		for (k = 0; k < nodesElement; k++){
            Residual[i] -= N[k]*UUp(k,i);
        }
	}

	for (i = 0; i < nScalarVariables; i++){
        for (k = 0; k < nodesElement; k++){
			FConv[i + k*nScalarVariables] = N[k]*L[i];
		}
    }

    // Build diffusive term: stress tensor and thermal diffusion

    
    // HERE FiND TAU AND q

    for (i = 0; i < sizeK;  i++)    K[i]  = 0.0;
    for (i = 0; i < sizeKT; i++)    KT[i] = 0.0;
			
	K[conta_new(0,0,0,0,3,3,6,3)] =  (-2.0 + ctau)*m1_el/DTOT2;
	K[conta_new(0,0,0,1,3,3,6,3)] =  ctau*m2_el/DTOT2;
	K[conta_new(0,0,0,2,3,3,6,3)] =  ctau*m3_el/DTOT2;
	
	K[conta_new(0,0,2,0,3,3,6,3)] =  (2.0 - ctau)/DTOT_el;
	K[conta_new(0,0,3,1,3,3,6,3)] = -ctau/DTOT_el;
	K[conta_new(0,0,4,2,3,3,6,3)] = -ctau/DTOT_el;
	
	K[conta_new(0,1,0,0,3,3,6,3)] = -m2_el/DTOT2;
	K[conta_new(0,1,0,1,3,3,6,3)] = -m1_el/DTOT2;
	K[conta_new(0,1,2,1,3,3,6,3)] =  1.0/DTOT_el; 
	K[conta_new(0,1,3,0,3,3,6,3)] =  1.0/DTOT_el;
	
	K[conta_new(0,2,0,0,3,3,6,3)] = -m3_el/DTOT2;
	K[conta_new(0,2,0,2,3,3,6,3)] = -m1_el/DTOT2;
	K[conta_new(0,2,2,2,3,3,6,3)] =  1.0/DTOT_el; 
	K[conta_new(0,2,4,0,3,3,6,3)] =  1.0/DTOT_el;
	
////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	K[conta_new(1,0,0,0,3,3,6,3)] = -m2_el/DTOT2;
	K[conta_new(1,0,0,1,3,3,6,3)] = -m1_el/DTOT2;
	K[conta_new(1,0,2,1,3,3,6,3)] =  1.0/DTOT_el;
	K[conta_new(1,0,3,0,3,3,6,3)] =  1.0/DTOT_el; 
	
	K[conta_new(1,1,0,0,3,3,6,3)] =  ctau*m1_el/DTOT2;
	K[conta_new(1,1,0,1,3,3,6,3)] =  (ctau - 2.0)*m2_el/DTOT2;
	K[conta_new(1,1,0,2,3,3,6,3)] =  ctau*m3_el/DTOT2;
	K[conta_new(1,1,2,0,3,3,6,3)] = -ctau/DTOT_el;
	K[conta_new(1,1,3,1,3,3,6,3)] =  (2.0 - ctau)/DTOT_el;
	K[conta_new(1,1,4,2,3,3,6,3)] = -ctau/DTOT_el;
	
	K[conta_new(1,2,0,1,3,3,6,3)] = -m3_el/DTOT2;
	K[conta_new(1,2,0,2,3,3,6,3)] = -m2_el/DTOT2;
	K[conta_new(1,2,3,2,3,3,6,3)] = 1.0/DTOT_el;
	K[conta_new(1,2,4,1,3,3,6,3)] = 1.0/DTOT_el;
	
	K[conta_new(2,0,0,0,3,3,6,3)] = -m3_el/DTOT2;
	K[conta_new(2,0,0,2,3,3,6,3)] = -m1_el/DTOT2;
	K[conta_new(2,0,2,2,3,3,6,3)] =  1.0/DTOT_el;
	K[conta_new(2,0,4,0,3,3,6,3)] =  1.0/DTOT_el; 
	
	K[conta_new(2,1,0,1,3,3,6,3)] = -m3_el/DTOT2;
	K[conta_new(2,1,0,2,3,3,6,3)] = -m2_el/DTOT2;
	K[conta_new(2,1,3,2,3,3,6,3)] = 1.0/DTOT_el;
	K[conta_new(2,1,4,1,3,3,6,3)] = 1.0/DTOT_el;
	
	K[conta_new(2,2,0,0,3,3,6,3)] =  ctau*m1_el/DTOT2;
	K[conta_new(2,2,0,1,3,3,6,3)] =  ctau*m2_el/DTOT2;
	K[conta_new(2,2,0,2,3,3,6,3)] =  (ctau - 2.0)*m3_el/DTOT2;
	K[conta_new(2,2,2,0,3,3,6,3)] = -ctau/DTOT_el;
	K[conta_new(2,2,3,1,3,3,6,3)] = -ctau/DTOT_el;
	K[conta_new(2,2,4,2,3,3,6,3)] = (2.0 - ctau)/DTOT_el;


	
	KT[0] =  pdg/(DG_el*R) - p_el/(DG_el*DG_el*R);
	KT[1] =  pds/(DG_el*R) + p_el/(DG_el*DG_el*R);
	KT[2] =  pm1/(DG_el*R);
	KT[3] =  pm2/(DG_el*R);
	KT[4] =  pm3/(DG_el*R);
	KT[5] =  pet/(DG_el*R);


    for (i = 0; i < SpaceDimension; i++){
        
        q(i) = 0.0;
		dt_diff(i) = 0.0;
		ds_diff(i) = 0.0;

        for (j = 0; j < SpaceDimension; j++){

            tau(i*SpaceDimension + j) = 0.0;
            
            for (l = 0; l < nScalarVariables; l++){
                for (m = 0; m < SpaceDimension; m++){

                    tau(i*SpaceDimension + j) += mu_mixture*K[conta_new(i,j,l,m,3,3,6,3)]*gradU[l*SpaceDimension + m];

                }
            }
        }
        for (l = 0; l < nScalarVariables; l++){
            q(i) -= lambda_mixture*KT(l)*gradU(l*SpaceDimension + i);
        }
    }

    ShockCapturing3d_new(mu_mixture, lambda_mixture, Cv_mixture, ros, h, dt_diff, ds_diff,tau, q, DTOT_el, gradU, Residual,U_gauss,a_el,norm_u,SpeedSound);

    // Build diffusive term: Diffusion tensor

	for ( i = 0; i < nScalarVariables*SpaceDimension; i++ )    G[i] = 0.0;
		
	for (i = 0; i < SpaceDimension; i++){

		for (j = 0; j < SpaceDimension; j++)
			G[(i + 2)*SpaceDimension + j] = -tau[i*SpaceDimension + j];
    }

	for (j = 0; j < SpaceDimension; j++){

		G[0*SpaceDimension + j] = -dt_diff[j];
		G[1*SpaceDimension + j] = -ds_diff[j];   // Numerical diffusivity for dust
		G[5*SpaceDimension + j] = q[j];

		for (i = 0; i < SpaceDimension; i++)
            G[5*SpaceDimension + j] += (-U_gauss[i + 2]/DTOT_el*tau[i*SpaceDimension + j]);

	}


    // Build diffusive term: Diffusion force

	for (s = 0; s < nNodalVariables; s++){

		FDiff[s] = 0.0;

		for (i = 0; i < nScalarVariables; i++){
			for ( j = 0; j < SpaceDimension; j++){
				FDiff[s] -= gradV[(SpaceDimension*i + j)*nNodalVariables + s]*G[i*SpaceDimension + j];
            }
        }
 	}

    // Stabilizing residual part

	for (s = 0; s < nNodalVariables; s++){

		FStab[s] = 0.0;

		for (k = 0; k < nScalarVariables; k++){
			FStab[s] += Lstar[s*nScalarVariables + k]*Residual[k]/invtauStab[k]*switchStab[k];
		}
	}

    // Force contribuution at the Gauss Point   

	int check = 1;

    for (i = 0; i < nNodalVariables; i++){

		F[i] = - sw_conv*FConv[i] - sw_diff*FDiff[i] - sw_stab*FStab[i];

		if (std::isnan(FConv[i]) == 1 || std::isnan(FDiff[i]) == 1 || std::isnan(FStab[i]) == 1){
			printf("%d %.3e %.3e %.3e\n", i, FConv[i], FDiff[i], FStab[i]);
			check = 0;
		}

	}
	if (check == 0)	{
		printf("%.3e %.3e %.3e %.3e %.3e %.3e \n", U_gauss(0), U_gauss(1), U_gauss(2), U_gauss(3), U_gauss(4), U_gauss(5));

		for (s = 0; s < nNodalVariables; s++){

			FStab[s] = 0.0;

			for (k = 0; k < nScalarVariables; k++){
				FStab[s] += Lstar[s*nScalarVariables + k]*Residual[k]/invtauStab[k]*switchStab[k];

				printf("%d %d %.3e %.3e %.3e\n", s, k, Lstar[s*nScalarVariables + k], Residual[k], invtauStab[k]);
			}
		}

		printf("stab_c2 = %.3e - norm_u = %.3e - SpeedSound = %.3e - stab_c1 = %.3e - mu_mixture = %.3e - lambda_mixture = %.3e"
			"Cp_mixture = %.3e \n", stab_c2, norm_u, SpeedSound, stab_c1, mu_mixture, lambda_mixture, Cp_mixture);

		double sps2 = DG_el*R*Cpmixed*(2*etot_el - DTOT_el*norm2u)/(2*DTOT_el*Cmixed2);
		
		printf("spsound2 = %.3e\n",sps2);

		printf("DG_el = %.3e - Cpmixed = %.3e - en-u2 = %.3e - DTOT_el = %.3e - Cmixed2 = %.3e\n",
		DG_el, Cpmixed, 2*etot_el - DTOT_el*norm2u, DTOT_el, Cmixed2);
		
		printf("%.3e %.3e \n",2*etot_el, DTOT_el*norm2u);

		printf("cos cos cos \n\n");
		


		abort();
	}

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector = F*data.volume /NGaussPoints; // Da vedere il fatto dei GP

    KRATOS_CATCH("")
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 2;
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 5;

    // Calculate the explicit residual vector
    BoundedVector<double, 15> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
        #pragma omp atomic
            r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        #pragma omp atomic
            r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY_SOLID) += rhs[aux +1];
            auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
            for (IndexType d = 0; d < dim; ++d) {
            #pragma omp atomic
                r_mom[d] += rhs[aux + 1 + (d + 1)];
            }
        #pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + dim + 2];
    }
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 3;
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 6;

    // Calculate the explicit residual vector
    BoundedVector<double, 24> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
        #pragma omp atomic
            r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        #pragma omp atomic
            r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY_SOLID) += rhs[aux +1];
            auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
            for (IndexType d = 0; d < dim; ++d) {
            #pragma omp atomic
                r_mom[d] += rhs[aux + 1 + (d + 1)];
            }
        #pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + dim + 2];
    }
}

template <>
void CompressibleNSBiphaseExplicit<2,3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 5;

    // Initialize and fill the mass matrix values
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six; rMassMatrix(0, 5) = one_twelve; rMassMatrix(0, 10) = one_twelve;
    rMassMatrix(1, 1) = one_six; rMassMatrix(1, 6) = one_twelve; rMassMatrix(1, 11) = one_twelve;
    rMassMatrix(2, 2) = one_six; rMassMatrix(2, 7) = one_twelve; rMassMatrix(2, 12) = one_twelve;
    rMassMatrix(3, 3) = one_six; rMassMatrix(3, 8) = one_twelve; rMassMatrix(3, 13) = one_twelve;
    rMassMatrix(4, 4) = one_six; rMassMatrix(4, 9) = one_twelve; rMassMatrix(4, 14) = one_twelve;

    rMassMatrix(5, 0) = one_twelve; rMassMatrix(5, 5) = one_six; rMassMatrix(5, 10) = one_twelve;
    rMassMatrix(6, 1) = one_twelve; rMassMatrix(6, 6) = one_six; rMassMatrix(6, 11) = one_twelve;
    rMassMatrix(7, 2) = one_twelve; rMassMatrix(7, 7) = one_six; rMassMatrix(7, 12) = one_twelve;
    rMassMatrix(8, 3) = one_twelve; rMassMatrix(8, 8) = one_six; rMassMatrix(8, 13) = one_twelve;
    rMassMatrix(9, 4) = one_twelve; rMassMatrix(9, 9) = one_six; rMassMatrix(9, 14) = one_twelve;

    rMassMatrix(10, 0) = one_twelve; rMassMatrix(10, 5) = one_twelve; rMassMatrix(10, 10) = one_six;
    rMassMatrix(11, 1) = one_twelve; rMassMatrix(11, 6) = one_twelve; rMassMatrix(11, 11) = one_six;
    rMassMatrix(12, 2) = one_twelve; rMassMatrix(12, 7) = one_twelve; rMassMatrix(12, 12) = one_six;
    rMassMatrix(13, 3) = one_twelve; rMassMatrix(13, 8) = one_twelve; rMassMatrix(13, 13) = one_six;
    rMassMatrix(14, 4) = one_twelve; rMassMatrix(14, 9) = one_twelve; rMassMatrix(14, 14) = one_six;
    

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

template <>
void CompressibleNSBiphaseExplicit<3,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 6;

    // Initialize and fill the mass matrix values
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    const unsigned int size = n_nodes * block_size;
    
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_ten; rMassMatrix(0, 6) = one_twenty; rMassMatrix(0, 12) =one_twenty; rMassMatrix(0,18) = one_twenty;
    rMassMatrix(1, 1) = one_ten; rMassMatrix(1, 7) = one_twenty; rMassMatrix(1, 13) = one_twenty; rMassMatrix(1,19) = one_twenty;
    rMassMatrix(2, 2) = one_ten; rMassMatrix(2, 8) = one_twenty; rMassMatrix(2, 14) = one_twenty; rMassMatrix(2,20) = one_twenty;
    rMassMatrix(3, 3) = one_ten; rMassMatrix(3, 9) = one_twenty; rMassMatrix(3, 15) = one_twenty; rMassMatrix(3,21) = one_twenty;
    rMassMatrix(4, 4) = one_ten; rMassMatrix(4, 10) = one_twenty; rMassMatrix(4, 16) = one_twenty; rMassMatrix(4,22) = one_twenty;
    rMassMatrix(5, 5) = one_ten; rMassMatrix(5, 11) = one_twenty; rMassMatrix(5, 17) = one_twenty; rMassMatrix(5,23) = one_twenty;

    rMassMatrix(6, 0) = one_twenty; rMassMatrix(6, 6) = one_ten; rMassMatrix(6, 12) =one_twenty; rMassMatrix(6,18) = one_twenty;
    rMassMatrix(7, 1) = one_twenty; rMassMatrix(7, 7) = one_ten; rMassMatrix(7, 13) = one_twenty; rMassMatrix(7,19) = one_twenty;
    rMassMatrix(8, 2) = one_twenty; rMassMatrix(8, 8) = one_ten; rMassMatrix(8, 14) = one_twenty; rMassMatrix(8,20) = one_twenty;
    rMassMatrix(9, 3) = one_twenty; rMassMatrix(9, 9) = one_ten; rMassMatrix(9, 15) = one_twenty; rMassMatrix(9,21) = one_twenty;
    rMassMatrix(10, 4) = one_twenty; rMassMatrix(10, 10) = one_ten; rMassMatrix(10, 16) = one_twenty; rMassMatrix(10,22) = one_twenty;
    rMassMatrix(11, 5) = one_twenty; rMassMatrix(11, 11) = one_ten; rMassMatrix(11, 17) = one_twenty; rMassMatrix(11,23) = one_twenty;

    rMassMatrix(12, 0) = one_twenty; rMassMatrix(12, 6) = one_twenty; rMassMatrix(12, 12) =one_ten; rMassMatrix(12,18) = one_twenty;
    rMassMatrix(13, 1) = one_twenty; rMassMatrix(13, 7) = one_twenty; rMassMatrix(13, 13) = one_ten; rMassMatrix(13,19) = one_twenty;
    rMassMatrix(14, 2) = one_twenty; rMassMatrix(14, 8) = one_twenty; rMassMatrix(14, 14) = one_ten; rMassMatrix(14,20) = one_twenty;
    rMassMatrix(15, 3) = one_twenty; rMassMatrix(15, 9) = one_twenty; rMassMatrix(15, 15) = one_ten; rMassMatrix(15,21) = one_twenty;
    rMassMatrix(16, 4) = one_twenty; rMassMatrix(16, 10) = one_twenty; rMassMatrix(16, 16) = one_ten; rMassMatrix(16,22) = one_twenty;
    rMassMatrix(17, 5) = one_twenty; rMassMatrix(17, 11) = one_twenty; rMassMatrix(17, 17) = one_ten; rMassMatrix(17,23) = one_twenty;

    rMassMatrix(18, 0) = one_twenty; rMassMatrix(18, 6) = one_twenty; rMassMatrix(18, 12) =one_twenty; rMassMatrix(18,18) = one_ten;
    rMassMatrix(19, 1) = one_twenty; rMassMatrix(18, 7) = one_twenty; rMassMatrix(19, 13) = one_twenty; rMassMatrix(19,19) = one_ten;
    rMassMatrix(20, 2) = one_twenty; rMassMatrix(20, 8) = one_twenty; rMassMatrix(20, 14) = one_twenty; rMassMatrix(20,20) = one_ten;
    rMassMatrix(21, 3) = one_twenty; rMassMatrix(21, 9) = one_twenty; rMassMatrix(21, 15) = one_twenty; rMassMatrix(21,21) = one_ten;
    rMassMatrix(22, 4) = one_twenty; rMassMatrix(22, 10) = one_twenty; rMassMatrix(22, 16) = one_twenty; rMassMatrix(22,22) = one_ten;
    rMassMatrix(23, 5) = one_twenty; rMassMatrix(23, 11) = one_twenty; rMassMatrix(23, 17) = one_twenty; rMassMatrix(23,23) = one_ten;
    

    // Here we assume that all the Gauss pt. have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Volume();
}

template <unsigned int TDim, unsigned int TNumNodes>
void CompressibleNSBiphaseExplicit<TDim, TNumNodes>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{   
    // Initialize the lumped mass vector
    constexpr IndexType size = TNumNodes * BlockSize;
    if (rLumpedMassVector.size() != BlockSize) {
        rLumpedMassVector.resize(size, false);
    }
    
    // Fill the lumped mass vector
    const double nodal_mass = GetGeometry().DomainSize() / TNumNodes;
    std::fill(rLumpedMassVector.begin(),rLumpedMassVector.end(),nodal_mass);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class CompressibleNSBiphaseExplicit<2,3>;
template class CompressibleNSBiphaseExplicit<3,4>;

}
