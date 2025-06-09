import numpy as np


def plot_sigma(ax, sigma_1, sigma_3):
    ax.plot(sigma_3, sigma_1, '-', color='blue', label='σ₁ vs σ₃')
    ax.set_title('σ₁ vs σ₃')
    ax.set_xlabel('σ₃ (Principal Stress 3) [kN/m²]')
    ax.set_ylabel('σ₁ (Principal Stress 1) [kN/m²]')
    ax.grid(True)

    min_val = 0
    max_val_x = max(sigma_3)
    max_val_y = min(sigma_1)
    padding_x = 0.1 * (max_val_x - min_val)
    padding_y = 0.1 * (max_val_y - min_val)

    ax.set_xlim(min_val, max_val_x + padding_x)
    ax.set_ylim(min_val, max_val_y + padding_y)

    ax.set_xticks(np.linspace(min_val, max_val_x, num=6))

def plot_delta_sigma(ax, vertical_strain, sigma_diff):
    ax.plot(vertical_strain, sigma_diff, '-', color='blue', label='|σ₁ - σ₃|')
    ax.set_title('|σ₁ - σ₃| vs εᵧᵧ')
    ax.set_xlabel('εᵧᵧ (Vertical Strain) [-]')
    ax.set_ylabel('|σ₁ - σ₃| [kN/m²]')
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_xticks(np.linspace(np.min(vertical_strain), np.max(vertical_strain), num=6))

def plot_volumetric_strain(ax, vertical_strain, volumetric_strain):
    ax.plot(vertical_strain, volumetric_strain, '-', color='blue', label='Volumetric Strain')
    ax.set_title('εᵥ vs εᵧᵧ')
    ax.set_xlabel('εᵧᵧ (Vertical Strain) [-]')
    ax.set_ylabel('εᵥ (Volumetric Strain) [-]')
    ax.grid(True)
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_xticks(np.linspace(np.min(vertical_strain), np.max(vertical_strain), num=6))

def plot_mohr_coulomb_circle(ax, sigma_1, sigma_3, cohesion=None, friction_angle=None):
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

def plot_p_q(ax, p_list, q_list):
    ax.plot(p_list, q_list, '-', color='blue', label="p' vs q")
    ax.set_title("p' vs q")
    ax.set_xlabel("p' (Mean Effective Stress) [kN/m²]")
    ax.set_ylabel("q (Deviatoric Stress) [kN/m²]")
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_xticks(np.linspace(np.min(p_list), np.max(p_list), num=6))
