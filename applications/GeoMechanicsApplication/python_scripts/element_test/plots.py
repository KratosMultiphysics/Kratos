import numpy as np
from ui_logger import log_message


def plot_principal_stresses_triaxial(ax, sigma_1, sigma_3):
    ax.plot(sigma_3, sigma_1, '-', color='blue', label='σ₁ vs σ₃')
    ax.set_title('σ₁ vs σ₃')
    ax.set_xlabel('σ₃ (Principal Stress 3) [kN/m²]')
    ax.set_ylabel('σ₁ (Principal Stress 1) [kN/m²]')
    ax.grid(True)
    ax.locator_params(nbins=8)

    min_val = 0
    max_val_x = max(sigma_3)
    max_val_y = min(sigma_1)
    padding_x = 0.1 * (max_val_x - min_val)
    padding_y = 0.1 * (max_val_y - min_val)
    ax.set_xlim(min_val, max_val_x + padding_x)
    ax.set_ylim(min_val, max_val_y + padding_y)
    ax.minorticks_on()

def plot_delta_sigma_triaxial(ax, vertical_strain, sigma_diff):
    ax.plot(vertical_strain, sigma_diff, '-', color='blue', label='|σ₁ - σ₃|')
    ax.set_title('|σ₁ - σ₃| vs εᵧᵧ')
    ax.set_xlabel('εᵧᵧ (Vertical Strain) [-]')
    ax.set_ylabel('|σ₁ - σ₃| [kN/m²]')
    ax.grid(True)
    ax.invert_xaxis()
    ax.locator_params(nbins=8)
    ax.minorticks_on()

def plot_volumetric_vertical_strain_triaxial(ax, vertical_strain, volumetric_strain):
    ax.plot(vertical_strain, volumetric_strain, '-', color='blue', label='Volumetric Strain')
    ax.set_title('εᵥ vs εᵧᵧ')
    ax.set_xlabel('εᵧᵧ (Vertical Strain) [-]')
    ax.set_ylabel('εᵥ (Volumetric Strain) [-]')
    ax.grid(True)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.locator_params(nbins=8)
    ax.minorticks_on()

def plot_mohr_coulomb_triaxial(ax, sigma_1, sigma_3, cohesion=None, friction_angle=None):
    if np.isclose(sigma_1, sigma_3):
        log_message("σ₁ is equal to σ₃. Mohr circle collapses to a point.", "warn")
    center = (sigma_1 + sigma_3) / 2
    radius = (sigma_1 - sigma_3) / 2
    theta = np.linspace(0, np.pi, 200)
    sigma = center + radius * np.cos(theta)
    tau = -radius * np.sin(theta)

    ax.plot(sigma, tau, label='Mohr-Coulomb', color='blue')

    if cohesion is not None and friction_angle is not None:
        phi_rad = np.radians(friction_angle)
        x_line = np.linspace(0, sigma_1, 200)
        y_line = x_line * np.tan(phi_rad) - cohesion
        ax.plot(x_line, -y_line, 'r--', label="Failure Criterion: τ = σ' tan(φ°) + c'")
        ax.legend(loc='upper left')

    ax.set_title("Mohr-Coulomb")
    ax.set_xlabel("σ' (Effective Stress) [kN/m²]")
    ax.set_ylabel("τ (Mobilized Shear Stress) [kN/m²]")
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_xlim(left=0, right= 1.2*np.max(sigma_1))
    ax.set_ylim(bottom=0, top = -0.6*np.max(sigma_1))
    ax.minorticks_on()

def plot_p_q_triaxial(ax, p_list, q_list):
    ax.plot(p_list, q_list, '-', color='blue', label="p' vs q")
    ax.set_title("p' vs q")
    ax.set_xlabel("p' (Mean Effective Stress) [kN/m²]")
    ax.set_ylabel("q (Deviatoric Stress) [kN/m²]")
    ax.grid(True)
    ax.invert_xaxis()
    ax.locator_params(nbins=8)
    ax.minorticks_on()

def plot_principal_stresses_direct_shear(ax, sigma_1, sigma_3):
    ax.plot(sigma_3, sigma_1, '-', color='blue', label='σ₁ vs σ₃')
    ax.set_title('σ₁ vs σ₃')
    ax.set_xlabel('σ₃ (Principal Stress 3) [kN/m²]')
    ax.set_ylabel('σ₁ (Principal Stress 1) [kN/m²]')
    ax.grid(True)
    ax.locator_params(nbins=8)

    min_x, max_x = min(sigma_3), max(sigma_3)
    min_y, max_y = min(sigma_1), max(sigma_1)
    x_padding = 0.1 * (max_x - min_x) if max_x > min_x else 1.0
    y_padding = 0.1 * (max_y - min_y) if max_y > min_y else 1.0
    ax.set_xlim(min_x - x_padding, max_x + x_padding)
    ax.set_ylim(min_y - y_padding, max_y + y_padding)

    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.legend(loc='best')
    ax.minorticks_on()

def plot_strain_stress_direct_shear(ax, shear_strain_xy, shear_stress_xy):
    gamma_xy = 2 * np.array(shear_strain_xy)
    ax.plot(np.abs(gamma_xy), np.abs(shear_stress_xy), '-', color='blue', label='τₓᵧ vs εₓᵧ')
    ax.set_title('τₓᵧ vs εₓᵧ')
    ax.set_xlabel('γₓᵧ (Shear Strain) [-]')
    ax.set_ylabel('τₓᵧ (Shear Stress) [kN/m²]')
    ax.grid(True)
    ax.locator_params(nbins=8)
    ax.minorticks_on()

def plot_mohr_coulomb_direct_shear(ax, sigma_1, sigma_3, cohesion=None, friction_angle=None):
    if np.isclose(sigma_1, sigma_3):
        log_message("σ₁ is equal to σ₃. Mohr circle collapses to a point.", "warn")
    center = (sigma_1 + sigma_3) / 2
    radius = (sigma_1 - sigma_3) / 2
    theta = np.linspace(0, np.pi, 400)
    sigma = center + radius * np.cos(theta)
    tau = -radius * np.sin(theta)

    ax.plot(sigma, tau, label='Mohr-Coulomb', color='blue')

    if cohesion is not None and friction_angle is not None:
        phi_rad = np.radians(friction_angle)
        max_sigma = center + radius
        x_line = np.linspace(0, max_sigma * 1.5, 400)
        y_line = x_line * np.tan(phi_rad) - cohesion
        ax.plot(x_line, -y_line, 'r--', label="Failure Criterion: τ = σ' tan(φ°) + c'")
        ax.legend(loc='upper left')

    ax.set_title("Mohr-Coulomb")
    ax.set_xlabel("σ' (Effective Stress) [kN/m²]")
    ax.set_ylabel("τ (Mobilized Shear Stress) [kN/m²]")
    ax.grid(True)
    ax.invert_xaxis()
    if sigma_1 > 0 or sigma_3 > 0:
        ax.set_xlim(left=1.2*np.max(sigma_3), right=1.2*np.max(sigma_1))
        ax.set_ylim(bottom=0, top = -0.9*(np.max(sigma_1) - np.max(sigma_3)))
    else:
        ax.set_xlim(left=0, right=1.2*np.max(sigma_1))
        ax.set_ylim(bottom=0, top = -0.9*np.max(sigma_1))
    ax.minorticks_on()

def plot_p_q_direct_shear(ax, p_list, q_list):
    ax.plot(p_list, q_list, '-', color='blue', label="p' vs q")
    ax.set_title("p' vs q")
    ax.set_xlabel("p' (Mean Effective Stress) [kN/m²]")
    ax.set_ylabel("q (Deviatoric Stress) [kN/m²]")
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_xlim(left=0, right=1.2*np.max(p_list))
    ax.locator_params(nbins=8)
    ax.minorticks_on()
