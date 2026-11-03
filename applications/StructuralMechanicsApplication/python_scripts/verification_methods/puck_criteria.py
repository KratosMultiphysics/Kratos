import numpy as np

def FibreFracture(sigma_1, R_pa_t, R_pa_c):
    strength = np.where(sigma_1 < 0.0, R_pa_c, R_pa_t)
    f_E = np.abs(sigma_1 / strength)
    mode = np.where(sigma_1 < 0.0, "FF Compression", "FF Tension")
    return f_E, mode

def InterFibreFracture(sigma_2, tau_21, R_tr_t, R_tr_c, R_trpa, p_trpa_t, p_trpa_c):
    root = np.sqrt(1.0 + 2.0 * p_trpa_c * R_tr_c / R_trpa)
    R_trtr_A = R_trpa / (2.0 * p_trpa_c) * (root - 1.0)
    p_trtr_c = 0.5 * (root - 1.0)
    tau_21_c = R_trpa * np.sqrt(1.0 + 2.0 * p_trtr_c)

    with np.errstate(divide="ignore", invalid="ignore"):
        mode_A = sigma_2 >= 0.0
        mode_B = (sigma_2 < 0.0) & (np.abs(sigma_2 / tau_21) <= R_trtr_A / np.abs(tau_21_c))
        mode_C = (sigma_2 < 0.0) & (np.abs(tau_21 / sigma_2) <= np.abs(tau_21_c) / R_trtr_A)

        f_E_A = (np.sqrt((tau_21 / R_trpa)**2 + (1.0 - p_trpa_t * R_tr_t / R_trpa)**2 * (sigma_2 / R_tr_t)**2)
                    + p_trpa_t * sigma_2 / R_trpa)
        f_E_B = (np.sqrt(tau_21**2 + (p_trpa_c * sigma_2)**2) + p_trpa_c * sigma_2) / R_trpa
        f_E_C = ((tau_21 / (2.0 * (1.0 + p_trtr_c) * R_trpa))**2 + (sigma_2 / R_tr_c)**2) * R_tr_c / -sigma_2

        theta_fp = np.degrees(np.arccos(np.sqrt(np.clip(
            (((R_trtr_A * tau_21) / (R_trpa * sigma_2))**2 + 1.0) / (2.0 * (1.0 + p_trtr_c)), 0.0, 1.0))))

        f_E = np.where(mode_A, f_E_A, np.where(mode_B, f_E_B, np.where(mode_C, f_E_C, np.nan)))
        mode = np.where(mode_A, "IFF Mode A", np.where(mode_B, "IFF Mode B", np.where(mode_C, "IFF Mode C", "NO IFF Mode")))
        theta_fp = np.where(mode_C, theta_fp, 0.0)
        return f_E, mode, theta_fp


def ReserveFactor(f_E):
    """RF = 1 / f_E; inf for an unloaded point (f_E = 0)."""
    with np.errstate(divide="ignore"):
        return 1.0 / f_E