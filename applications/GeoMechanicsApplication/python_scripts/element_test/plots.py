import numpy as np
from ui_logger import log_message
from ui_labels import (
    SIGMA1_LABEL, SIGMA3_LABEL, SIGMA1_SIGMA3_DIFF_LABEL, VERTICAL_STRAIN_LABEL,
    VOLUMETRIC_STRAIN_LABEL, SHEAR_STRAIN_LABEL, SHEAR_STRESS_LABEL, EFFECTIVE_STRESS_LABEL,
    MOBILIZED_SHEAR_STRESS_LABEL, P_STRESS_LABEL, Q_STRESS_LABEL, TITLE_SIGMA1_VS_SIGMA3,
    TITLE_DIFF_PRINCIPAL_SIGMA_VS_STRAIN, TITLE_VOL_VS_VERT_STRAIN, TITLE_MOHR, TITLE_P_VS_Q, TITLE_SHEAR_VS_STRAIN,
    LEGEND_MC, LEGEND_MC_FAILURE
)


def plot_principal_stresses_triaxial(ax, sigma_1, sigma_3):
    ax.plot(sigma_3, sigma_1, '-', color='blue', label=TITLE_SIGMA1_VS_SIGMA3)
    ax.set_title(TITLE_SIGMA1_VS_SIGMA3)
    ax.set_xlabel(SIGMA3_LABEL)
    ax.set_ylabel(SIGMA1_LABEL)
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
    ax.plot(vertical_strain, sigma_diff, '-', color='blue', label=SIGMA1_SIGMA3_DIFF_LABEL)
    ax.set_title(TITLE_DIFF_PRINCIPAL_SIGMA_VS_STRAIN)
    ax.set_xlabel(VERTICAL_STRAIN_LABEL)
    ax.set_ylabel(SIGMA1_SIGMA3_DIFF_LABEL)
    ax.grid(True)
    ax.invert_xaxis()
    ax.locator_params(nbins=8)
    ax.minorticks_on()

def plot_volumetric_vertical_strain_triaxial(ax, vertical_strain, volumetric_strain):
    ax.plot(vertical_strain, volumetric_strain, '-', color='blue', label=TITLE_VOL_VS_VERT_STRAIN)
    ax.set_title(TITLE_VOL_VS_VERT_STRAIN)
    ax.set_xlabel(VERTICAL_STRAIN_LABEL)
    ax.set_ylabel(VOLUMETRIC_STRAIN_LABEL)
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

    ax.plot(sigma, tau, label=LEGEND_MC, color='blue')

    if cohesion is not None and friction_angle is not None:
        phi_rad = np.radians(friction_angle)
        x_line = np.linspace(0, sigma_1, 200)
        y_line = x_line * np.tan(phi_rad) - cohesion
        ax.plot(x_line, -y_line, 'r--', label=LEGEND_MC_FAILURE)
        ax.legend(loc='upper left')

    ax.set_title(LEGEND_MC)
    ax.set_xlabel(EFFECTIVE_STRESS_LABEL)
    ax.set_ylabel(MOBILIZED_SHEAR_STRESS_LABEL)
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_xlim(left=0, right= 1.2*np.max(sigma_1))
    ax.set_ylim(bottom=0, top = -0.6*np.max(sigma_1))
    ax.minorticks_on()

def plot_p_q_triaxial(ax, p_list, q_list):
    ax.plot(p_list, q_list, '-', color='blue', label=TITLE_P_VS_Q)
    ax.set_title(TITLE_P_VS_Q)
    ax.set_xlabel(P_STRESS_LABEL)
    ax.set_ylabel(Q_STRESS_LABEL)
    ax.grid(True)
    ax.invert_xaxis()
    ax.locator_params(nbins=8)
    ax.minorticks_on()

def plot_principal_stresses_direct_shear(ax, sigma_1, sigma_3):
    ax.plot(sigma_3, sigma_1, '-', color='blue', label=TITLE_SIGMA1_VS_SIGMA3)
    ax.set_title(TITLE_SIGMA1_VS_SIGMA3)
    ax.set_xlabel(SIGMA3_LABEL)
    ax.set_ylabel(SIGMA1_LABEL)
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
    ax.minorticks_on()

def plot_strain_stress_direct_shear(ax, shear_strain_xy, shear_stress_xy):
    gamma_xy = 2 * np.array(shear_strain_xy)
    ax.plot(np.abs(gamma_xy), np.abs(shear_stress_xy), '-', color='blue', label=TITLE_SHEAR_VS_STRAIN)
    ax.set_title(TITLE_SHEAR_VS_STRAIN)
    ax.set_xlabel(SHEAR_STRAIN_LABEL)
    ax.set_ylabel(SHEAR_STRESS_LABEL)
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

    ax.plot(sigma, tau, label=LEGEND_MC, color='blue')

    if cohesion is not None and friction_angle is not None:
        phi_rad = np.radians(friction_angle)
        max_sigma = center + radius
        x_line = np.linspace(0, max_sigma * 1.5, 400)
        y_line = x_line * np.tan(phi_rad) - cohesion
        ax.plot(x_line, -y_line, 'r--', label=LEGEND_MC_FAILURE)
        ax.legend(loc='upper left')

    ax.set_title(LEGEND_MC)
    ax.set_xlabel(EFFECTIVE_STRESS_LABEL)
    ax.set_ylabel(MOBILIZED_SHEAR_STRESS_LABEL)
    ax.grid(True)
    ax.invert_xaxis()

    epsilon = 0.1
    relative_diff = np.abs(sigma_1 - sigma_3) / max(np.abs(sigma_1), 1e-6)

    if relative_diff < epsilon:
        ax.set_xlim(center - (1.2 * radius), center + (1.2 * radius))
        ax.set_ylim(bottom=0, top = -0.9*(np.max(sigma_1) - np.max(sigma_3)))

    else:
        if sigma_1 > 0 or sigma_3 > 0:
            ax.set_xlim(left=1.2*np.max(sigma_3), right=1.2*np.max(sigma_1))
            ax.set_ylim(bottom=0, top = -0.9*(np.max(sigma_1) - np.max(sigma_3)))
        else:
            ax.set_xlim(left=0, right=1.2*np.max(sigma_1))
            ax.set_ylim(bottom=0, top = -0.9*np.max(sigma_1))


    ax.minorticks_on()

def plot_p_q_direct_shear(ax, p_list, q_list):
    ax.plot(p_list, q_list, '-', color='blue', label=TITLE_P_VS_Q)
    ax.set_title(TITLE_P_VS_Q)
    ax.set_xlabel(P_STRESS_LABEL)
    ax.set_ylabel(Q_STRESS_LABEL)
    ax.grid(True)
    ax.invert_xaxis()
    ax.set_xlim(left=0, right=1.2*np.max(p_list))
    ax.locator_params(nbins=8)
    ax.minorticks_on()
