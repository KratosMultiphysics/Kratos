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

#if !defined(KRATOS_GAUSS_POINT_ERROR_PROCESS_H_INCLUDED )
#define  KRATOS_GAUSS_POINT_ERROR_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/gid_io.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/binbased_fast_point_locator.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"

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

/// Short class definition.
/** Detail class definition.
*/
class GaussPointErrorProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GaussPointErrorProcess
    KRATOS_CLASS_POINTER_DEFINITION(GaussPointErrorProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    // GaussPointErrorProcess(
    //     ModelPart& rModelPart,
    //     unsigned int DomainSize)
    //     : mrModelPart(rModelPart),
    //       mDomainSize(DomainSize)
    // {
    // }

    GaussPointErrorProcess(
        ModelPart& rModelPart,
        ModelPart& rReferenceModelPart,
        unsigned int DomainSize)
        : mrModelPart(rModelPart),
          mrReferenceModelPart(rReferenceModelPart),
          mDomainSize(DomainSize)
    {
    }

    /// Destructor.
    ~GaussPointErrorProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void ExecuteFinalize() override
    {
        KRATOS_TRY

        // Variables initialization
        double cut_area = 0.0;
        double total_area = 0.0;
        double cut_error_p_norm = 0.0;
        double cut_error_v_norm = 0.0;
        double total_error_p_norm = 0.0;
        double total_error_v_norm = 0.0;

        KRATOS_WATCH("gauss point error process");

        // Map reference values from reference mesh to analysis mesh
        BinBasedFastPointLocator<2> point_locator = BinBasedFastPointLocator<2>(mrReferenceModelPart);
        point_locator.UpdateSearchDatabase();

        // Loop the elements to compute the error error in each Gauss pt.
        if(mDomainSize == 2)
        {
            for(auto it_elem = mrModelPart.ElementsBegin(); it_elem != mrModelPart.ElementsEnd(); ++it_elem)
            {
                auto &r_geom = it_elem->GetGeometry();
                const auto elem_dist = it_elem->GetValue(ELEMENTAL_DISTANCES);

                const auto &r_p_0 = r_geom[0].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_1 = r_geom[1].FastGetSolutionStepValue(PRESSURE);
                const auto &r_p_2 = r_geom[2].FastGetSolutionStepValue(PRESSURE);
                const auto &r_vel_0 = r_geom[0].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_1 = r_geom[1].FastGetSolutionStepValue(VELOCITY);
                const auto &r_vel_2 = r_geom[2].FastGetSolutionStepValue(VELOCITY);

                // unsigned int aux_n = 0;
                // const double min_y = -1.0;
                // for (unsigned int i = 0; i < r_geom.PointsNumber(); ++i) {
                //     if (r_geom[i].Y() > min_y) {
                //         aux_n++;
                //     }
                // }

                // const bool is_valid = aux_n == r_geom.PointsNumber() ? true : false;
                const bool is_valid = true;

                // Split element case
                if (is_valid) {
                    if (this->IsCut(elem_dist)){
                        // Generate the splitting pattern
                        DivideTriangle2D3<Node>::Pointer p_divide_util = Kratos::make_shared<DivideTriangle2D3<Node>>(r_geom, elem_dist);
                        p_divide_util->GenerateDivision();

                        // Generate the modified shape functions util
                        auto p_geom = it_elem->pGetGeometry();
                        ModifiedShapeFunctions::Pointer p_ausas_modified_sh_func = Kratos::make_shared<Triangle2D3AusasModifiedShapeFunctions>(p_geom, elem_dist);

                        Matrix N_pos_side;
                        Vector w_gauss_pos_side;
                        GeometryType::ShapeFunctionsGradientsType DN_DX_pos_side;
                        p_ausas_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(N_pos_side, DN_DX_pos_side, w_gauss_pos_side, GeometryData::IntegrationMethod::GI_GAUSS_2);

                        // Save the coordinates of all the subdivision Gauss pts.
                        const unsigned int n_pos_gauss_pt = w_gauss_pos_side.size();
                        Matrix pos_gauss_pts_coords = ZeroMatrix(n_pos_gauss_pt, 2);
                        const unsigned int n_pos_sbdiv = (p_divide_util->GetPositiveSubdivisions()).size();

                        for (unsigned int i_sbdiv = 0; i_sbdiv < n_pos_sbdiv; ++i_sbdiv){
                            auto p_sbdiv = (p_divide_util->GetPositiveSubdivisions())[i_sbdiv];
                            auto N = p_sbdiv->ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
                            for (unsigned int i_gauss = 0;  i_gauss < 3; ++i_gauss){
                                pos_gauss_pts_coords(i_sbdiv*3+i_gauss,0) = ((*p_sbdiv)[0]).X()*N(i_gauss,0) + ((*p_sbdiv)[1]).X()*N(i_gauss,1) + ((*p_sbdiv)[2]).X()*N(i_gauss,2);
                                pos_gauss_pts_coords(i_sbdiv*3+i_gauss,1) = ((*p_sbdiv)[0]).Y()*N(i_gauss,0) + ((*p_sbdiv)[1]).Y()*N(i_gauss,1) + ((*p_sbdiv)[2]).Y()*N(i_gauss,2);
                            }
                        }

                        // Compute the error in each Gauss pt.
                        for (unsigned int i_gauss = 0; i_gauss < n_pos_gauss_pt; ++i_gauss){
                            // Locate the Gauss pt. in the reference mesh
                            array_1d<double,3> coords = ZeroVector(3);
                            coords[0] = pos_gauss_pts_coords(i_gauss,0);
                            coords[1] = pos_gauss_pts_coords(i_gauss,1);
                            Vector N_vect;
                            Element::Pointer p_elem;
                            const bool is_found = point_locator.FindPointOnMeshSimplified(coords, N_vect, p_elem);
                            KRATOS_WARNING_IF("GaussPointErrorProcess", is_found != true) << "Gauss pt. x: " << coords[0] << " y: " << coords[1] << " not found in element " << it_elem->Id() << std::endl;

                            // Compute the body-fitted values in the embedded mesh
                            double p_exact = 0.0;
                            array_1d<double,3> v_exact = ZeroVector(3);

                            if (is_found){
                                const auto &r_geom_elem = p_elem->GetGeometry();
                                for (unsigned int i_node = 0; i_node < r_geom_elem.PointsNumber(); ++i_node){
                                    v_exact += r_geom_elem[i_node].GetSolutionStepValue(VELOCITY) * N_vect[i_node];
                                    p_exact += r_geom_elem[i_node].GetSolutionStepValue(PRESSURE) * N_vect[i_node];
                                }
                            }

                            // Exact solution
                            // const double x_coord = pos_gauss_pts_coords(i_gauss,0);
                            // const double y_coord = pos_gauss_pts_coords(i_gauss,1);
                            // const double p_exact = this->ComputePressureExactSolution(x_coord, y_coord);
                            // const array_1d<double, 3> v_exact = this->ComputeVelocityExactSolution(x_coord, y_coord);

                            // Obtained solution
                            array_1d<double, 3> v_solu;
                            const double p_solu = r_p_0*N_pos_side(i_gauss,0) + r_p_1*N_pos_side(i_gauss,1) + r_p_2*N_pos_side(i_gauss,2);
                            v_solu[0] = r_vel_0[0]*N_pos_side(i_gauss,0) + r_vel_1[0]*N_pos_side(i_gauss,1) + r_vel_2[0]*N_pos_side(i_gauss,2);
                            v_solu[1] = r_vel_0[1]*N_pos_side(i_gauss,0) + r_vel_1[1]*N_pos_side(i_gauss,1) + r_vel_2[1]*N_pos_side(i_gauss,2);
                            v_solu[2] = r_vel_0[2]*N_pos_side(i_gauss,0) + r_vel_1[2]*N_pos_side(i_gauss,1) + r_vel_2[2]*N_pos_side(i_gauss,2);

                            // Velocity error norm
                            cut_error_v_norm += w_gauss_pos_side[i_gauss]*inner_prod(v_exact - v_solu, v_exact - v_solu);
                            total_error_v_norm += w_gauss_pos_side[i_gauss]*inner_prod(v_exact - v_solu, v_exact - v_solu);

                            // Pressure error norm
                            cut_error_p_norm += w_gauss_pos_side[i_gauss]*std::pow(p_exact - p_solu, 2);
                            total_error_p_norm += w_gauss_pos_side[i_gauss]*std::pow(p_exact - p_solu, 2);

                            // Add the local Gauss contribution to the areas
                            cut_area += w_gauss_pos_side[i_gauss];
                            total_area += w_gauss_pos_side[i_gauss];
                        }

                    // No split element case
                    } else {
                        // Check if the element is in the positive (fluid) region
                        if (this->IsPositive(elem_dist)){
                            // Get geometry data
                            Vector jac_vect;
                            jac_vect = r_geom.DeterminantOfJacobian(jac_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
                            auto N = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

                            // Compute the error in each Gauss pt.
                            for (unsigned int i_gauss = 0; i_gauss < 3; ++i_gauss){
                                // Locate the Gauss pt. in the reference mesh
                                array_1d<double,3> coords = ZeroVector(3);
                                coords[0] = r_geom[0].X()*N(i_gauss,0) + r_geom[1].X()*N(i_gauss,1) + r_geom[2].X()*N(i_gauss,2);
                                coords[1] = r_geom[0].Y()*N(i_gauss,0) + r_geom[1].Y()*N(i_gauss,1) + r_geom[2].Y()*N(i_gauss,2);
                                Vector N_vect;
                                Element::Pointer p_elem;
                                const bool is_found = point_locator.FindPointOnMeshSimplified(coords, N_vect, p_elem);
                                KRATOS_WARNING_IF("GaussPointErrorProcess", is_found != true) << "Gauss pt. x: " << coords[0] << " y: " << coords[1] << " not found!" << std::endl;

                                // Compute the body-fitted values in the embedded mesh
                                double p_exact = 0.0;
                                array_1d<double,3> v_exact = ZeroVector(3);

                                if (is_found){
                                    const auto &r_geom_elem = p_elem->GetGeometry();
                                    for (unsigned int i_node = 0; i_node < r_geom_elem.PointsNumber(); ++i_node){
                                        v_exact += r_geom_elem[i_node].GetSolutionStepValue(VELOCITY) * N_vect[i_node];
                                        p_exact += r_geom_elem[i_node].GetSolutionStepValue(PRESSURE) * N_vect[i_node];
                                    }
                                }

                                // Exact solution
                                // const double x_coord = r_geom[0].X()*N(i_gauss,0) + r_geom[1].X()*N(i_gauss,1) + r_geom[2].X()*N(i_gauss,2);
                                // const double y_coord = r_geom[0].Y()*N(i_gauss,0) + r_geom[1].Y()*N(i_gauss,1) + r_geom[2].Y()*N(i_gauss,2);
                                // const double p_exact = this->ComputePressureExactSolution(x_coord, y_coord);
                                // const array_1d<double, 3> v_exact = this->ComputeVelocityExactSolution(x_coord, y_coord);

                                // Obtained solution
                                array_1d<double, 3> v_solu;
                                const double p_solu = r_p_0*N(i_gauss,0) + r_p_1*N(i_gauss,1) + r_p_2*N(i_gauss,2);
                                v_solu[0] = r_vel_0[0]*N(i_gauss,0) + r_vel_1[0]*N(i_gauss,1) + r_vel_2[0]*N(i_gauss,2);
                                v_solu[1] = r_vel_0[1]*N(i_gauss,0) + r_vel_1[1]*N(i_gauss,1) + r_vel_2[1]*N(i_gauss,2);
                                v_solu[2] = r_vel_0[2]*N(i_gauss,0) + r_vel_1[2]*N(i_gauss,1) + r_vel_2[2]*N(i_gauss,2);

                                // Velocity error norm
                                total_error_v_norm += (jac_vect[i_gauss] / 6.0)*inner_prod(v_exact - v_solu, v_exact - v_solu);

                                // Pressure error norm
                                total_error_p_norm += (jac_vect[i_gauss] / 6.0)*std::pow(p_exact - p_solu, 2);

                                // Add the local Gauss contribution to the areas
                                total_area += (jac_vect[i_gauss] / 6.0);
                            }
                        }
                    }
                }
            }
        }
        else if(mDomainSize == 3)
        {
            KRATOS_ERROR << "3D case not implemented yet." << std::endl;
        }

        // Compute the square root of the norms
        cut_error_p_norm = std::sqrt(cut_error_p_norm);
        cut_error_v_norm = std::sqrt(cut_error_v_norm);
        total_error_p_norm = std::sqrt(total_error_p_norm);
        total_error_v_norm = std::sqrt(total_error_v_norm);

        // Print the obtained values
        std::cout << std::endl;
        std::cout << "Cut area: " << cut_area << std::endl;
        std::cout << "Total area: " << total_area << std::endl;
        std::cout << "Cut error p_norm: " << cut_error_p_norm << std::endl;
        std::cout << "Total error p_norm: " << total_error_p_norm << std::endl;
        std::cout << "Cut error v_norm: " << cut_error_v_norm << std::endl;
        std::cout << "Total error v_norm: " << total_error_v_norm << std::endl;
        std::cout << std::endl;

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
        return "GaussPointErrorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GaussPointErrorProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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


    ///@}
    ///@name Protected Operations
    ///@{

    bool IsCut(const Vector& rElementalDistances)
    {
        unsigned int n_pos = 0;
        unsigned int n_neg = 0;

        for (unsigned int i_node = 0; i_node < 3; ++i_node){
            if(rElementalDistances[i_node] > 0.0) {
                n_pos++;
            } else {
                n_neg++;
            }
        }

        if (n_pos != 0 && n_neg != 0){
            return true;
        } else {
            return false;
        }
    }

    bool IsPositive(const Vector& rElementalDistances)
    {
        unsigned int n_pos = 0;

        for (unsigned int i_node = 0; i_node < 3; ++i_node){
            if(rElementalDistances[i_node] > 0.0) {
                n_pos++;
            }
        }

        if (n_pos == 3){
            return true;
        } else {
            return false;
        }
    }

    array_1d<double,3> ComputeVelocityExactSolution(
        const double X,
        const double Y)
    {
        double tetha = atan2(Y,X);
        if (tetha < 0.0){
            tetha = 2.0*M_PI + tetha; // atan2 3rd and 4th quadrants correction
        }
        const double radius = std::sqrt(std::pow(X,2) + std::pow(Y,2));

        array_1d<double, 3> velocity(0.0);
        velocity[0] = radius*sin(tetha); // x-velocity component
        velocity[1] = -radius*cos(tetha);  // y_velocity component

        // set velocity
        return velocity;
    }

    double ComputePressureExactSolution(
        const double X,
        const double Y)
    {
        const double rho = 1.0;
        const double squared_radius = (std::pow(X,2) + std::pow(Y,2));
        return 0.5*rho*(squared_radius - 1.0);
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

    ModelPart& mrModelPart;
    ModelPart& mrReferenceModelPart;
    unsigned int mDomainSize;

    ///@}
    ///@name Private Operators
    ///@{


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
    GaussPointErrorProcess& operator=(GaussPointErrorProcess const& rOther);

    /// Copy constructor.
    //GaussPointErrorProcess(GaussPointErrorProcess const& rOther);


    ///@}

}; // Class GaussPointErrorProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GaussPointErrorProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GaussPointErrorProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GAUSS_POINT_ERROR_PROCESS_H_INCLUDED  defined
