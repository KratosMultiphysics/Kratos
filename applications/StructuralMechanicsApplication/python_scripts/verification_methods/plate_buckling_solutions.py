"""Buckling of simply supported (SSSS) rectangular plates with symmetric layup (B = 0).

Sign convention: compression positive (N0 = -N for N tension positive).
All functions return the critical load factor Lambda: buckling occurs at Lambda * N0.
"""
import numpy as np
import scipy.linalg


def NavierSpeciallyOrthotropic(D: np.ndarray, a: float, b: float, N0x: float, N0y: float,
                               max_half_waves: int = 30) -> tuple[float, int, int]:
    """Exact solution for D16 = D26 = 0 and no shear load. Returns (Lambda, m, n)."""
    m, n = np.meshgrid(np.arange(1, max_half_waves + 1), np.arange(1, max_half_waves + 1), indexing="ij")
    km2, kn2 = (m / a)**2, (n / b)**2
    stiffness = np.pi**2 * (D[0, 0] * km2**2 + 2.0 * (D[0, 1] + 2.0 * D[2, 2]) * km2 * kn2 + D[1, 1] * kn2**2)
    load = N0x * km2 + N0y * kn2
    lam = np.where(load > 0.0, stiffness / np.where(load > 0.0, load, 1.0), np.inf)
    i, j = np.unravel_index(np.argmin(lam), lam.shape)
    return float(lam[i, j]), int(m[i, j]), int(n[i, j])


def _SineIntegrals(size: int) -> dict:
    """Integrals over [0, 1] of d^i(sin(m pi x)) * d^j(sin(p pi x)), indices m, p = 1..size."""
    k = np.arange(1, size + 1, dtype=float)
    m, p = np.meshgrid(k, k, indexing="ij")
    odd = (m + p) % 2 == 1
    I01 = np.zeros((size, size))
    I01[odd] = 2.0 * m[odd] * p[odd] / (m[odd]**2 - p[odd]**2)
    Ic = np.zeros((size, size))
    Ic[odd] = -4.0 * m[odd] * p[odd] / (np.pi**2 * (m[odd]**2 - p[odd]**2)**2)
    return {(0, 0): np.eye(size) / 2.0,
            (1, 1): np.diag((k * np.pi)**2 / 2.0),
            (2, 0): np.diag(-(k * np.pi)**2 / 2.0),
            (0, 2): np.diag(-(k * np.pi)**2 / 2.0),
            (2, 2): np.diag((k * np.pi)**4 / 2.0),
            (0, 1): I01,
            (1, 0): I01.T,
            (2, 1): -(m * np.pi)**2 * I01, "centered": Ic}


def RitzSSSS(D: np.ndarray, a: float, b: float, N0x: float, N0y: float, N0xy: float,
             n_terms: int = 20, dN0x: float = 0.0, dN0y: float = 0.0) -> tuple[float, np.ndarray | None]:
    """Ritz method with n_terms x n_terms sine functions, including D16, D26, shear and linear load gradients.

    N0x(eta) = N0x + dN0x (eta - 1/2): N_xx varying across the width (in-plane bending)
    N0y(xi)  = N0y + dN0y (xi - 1/2):  N_yy varying along the length
    Returns (Lambda, W) with the mode coefficients W[m-1, n-1]; (inf, None) if nothing buckles.
    """
    D11, D12, D22, D16, D26, D66 = D[0, 0], D[0, 1], D[1, 1], D[0, 2], D[1, 2], D[2, 2]
    ref = np.sqrt(D11 * D22)
    alpha = (a / b) * (D22 / D11)**0.25
    beta = (D12 + 2.0 * D66) / ref
    zeta = D12 / (D12 + 2.0 * D66)
    gamma = D16 / (D11**0.75 * D22**0.25)
    delta = D26 / (D11**0.25 * D22**0.75)
    nxx = N0x * b**2 / (np.pi**2 * ref)
    nyy = N0y * a**2 / (np.pi**2 * ref)
    nxy = N0xy * a * b / (np.pi**2 * ref)
    nxx_grad = dN0x * b**2 / (np.pi**2 * ref)
    nyy_grad = dN0y * a**2 / (np.pi**2 * ref)

    I = _SineIntegrals(n_terms)
    K = (1.0 / (2.0 * alpha**2) * np.kron(I[2, 2], I[0, 0])
         + alpha**2 / 2.0       * np.kron(I[0, 0], I[2, 2])
         + beta * zeta          * np.kron(I[2, 0], I[0, 2])
         + (1.0 - zeta) * beta  * np.kron(I[1, 1], I[1, 1])
         + 2.0 * gamma / alpha  * np.kron(I[2, 1], I[0, 1])
         + 2.0 * delta * alpha  * np.kron(I[0, 1], I[2, 1]))
    G = (np.pi**2 / 2.0 * np.kron(I[1, 1], nxx * I[0, 0] + nxx_grad * I["centered"])
         + np.pi**2 / 2.0 * np.kron(nyy * I[0, 0] + nyy_grad * I["centered"], I[1, 1])
         + np.pi**2       * nxy * np.kron(I[1, 0], I[0, 1]))
    K = 0.5 * (K + K.T)     # only the symmetric part belongs to the strain energy
    G = 0.5 * (G + G.T)


    size = n_terms**2
    mu, v = scipy.linalg.eigh(G, K, subset_by_index=[size - 1, size - 1])
    if mu[0] <= 0.0:
        return np.inf, None
    return 1.0 / mu[0], v[:, 0].reshape(n_terms, n_terms)     # W[m-1, n-1]


def RitzModeShape(W: np.ndarray, xi: np.ndarray, eta: np.ndarray) -> np.ndarray:
    """Mode shape w(xi, eta) on a grid, normalized to max |w| = 1. xi = x/a, eta = y/b in [0, 1]."""
    m = np.arange(1, W.shape[0] + 1)
    n = np.arange(1, W.shape[1] + 1)
    w = np.sin(np.pi * np.outer(xi, m)) @ W @ np.sin(np.pi * np.outer(n, eta))
    return w / np.abs(w).max()
