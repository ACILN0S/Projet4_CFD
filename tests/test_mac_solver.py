"""Tests de la méthode 1 (grille MAC, projection de Chorin) et de ses diagnostics.

  1. Opérateurs sur un champ analytique → ordre 2 (laplacien, advection centrée), ordre 1 (upwind 1).
  2. Énergie cinétique K d'un champ analytique → ordre 2.
  3. Après projection, max |∇·u| est au niveau de l'erreur machine.
  4. Bilan d'énergie dK/dt = P − ε sur un run court à Re = 100 (advection centrée).

Champ test (type Taylor–Green, à divergence nulle, u = 0 en x = 0, 1 et v = 0 en y = 0, 1) :
    u =  sin(πx) cos(πy),   v = −cos(πx) sin(πy)
    Δu = −2π² u,  Δv = −2π² v
    A_u = ∂(uu)/∂x + ∂(vu)/∂y = (π/2) sin(2πx),   A_v = (π/2) sin(2πy)
    K = ½ ∫ (u² + v²) dS = 1/4
"""
import numpy as np

from cavity.diagnostics import dissipation, kinetic_energy, lid_power
from cavity.mac_solver import (advection, apply_bc, divergence, laplacian, new_fields,
                               project)
from cavity.params import Config, lid_velocity
from cavity.run import run

GRIDS = [16, 32, 64]


def observed_order(hs, errors):
    """Pente de log(erreur) en fonction de log(h)."""
    return np.polyfit(np.log(hs), np.log(errors), 1)[0]


def taylor_green(N):
    """Champ test placé aux positions MAC, fantômes compris (valeurs analytiques partout).

    u[i, jj] est en (x, y) = (i h, (jj − ½) h) ;  v[ii, j] est en ((ii − ½) h, j h).
    """
    h = 1.0 / N
    xu, yu = np.arange(N + 1) * h, (np.arange(N + 2) - 0.5) * h
    xv, yv = (np.arange(N + 2) - 0.5) * h, np.arange(N + 1) * h
    Xu, Yu = np.meshgrid(xu, yu, indexing="ij")
    Xv, Yv = np.meshgrid(xv, yv, indexing="ij")
    u = np.sin(np.pi * Xu) * np.cos(np.pi * Yu)
    v = -np.cos(np.pi * Xv) * np.sin(np.pi * Yv)
    return h, u, v, Xu, Yv


def test_laplacian_order2():
    errors, hs = [], []
    for N in GRIDS:
        h, u, v, _, _ = taylor_green(N)
        Lu, Lv = np.zeros_like(u), np.zeros_like(v)
        laplacian(u, v, h, Lu, Lv)
        # points intérieurs : u → i = 1..N-1, jj = 1..N ; v → ii = 1..N, j = 1..N-1
        err_u = np.max(np.abs(Lu[1:-1, 1:-1] + 2 * np.pi**2 * u[1:-1, 1:-1]))
        err_v = np.max(np.abs(Lv[1:-1, 1:-1] + 2 * np.pi**2 * v[1:-1, 1:-1]))
        errors.append(max(err_u, err_v))
        hs.append(h)
    assert abs(observed_order(hs, errors) - 2.0) < 0.1


def advection_errors(upwind):
    errors, hs = [], []
    for N in GRIDS:
        h, u, v, Xu, Yv = taylor_green(N)
        Au, Av = np.zeros_like(u), np.zeros_like(v)
        advection(u, v, h, upwind, Au, Av)
        err_u = np.max(np.abs(Au[1:-1, 1:-1] - 0.5 * np.pi * np.sin(2 * np.pi * Xu[1:-1, 1:-1])))
        err_v = np.max(np.abs(Av[1:-1, 1:-1] - 0.5 * np.pi * np.sin(2 * np.pi * Yv[1:-1, 1:-1])))
        errors.append(max(err_u, err_v))
        hs.append(h)
    return hs, errors


def test_advection_centered_order2():
    assert abs(observed_order(*advection_errors(upwind=False)) - 2.0) < 0.15


def test_advection_upwind1_order1():
    assert abs(observed_order(*advection_errors(upwind=True)) - 1.0) < 0.15


def test_kinetic_energy_order2():
    """Champ polynomial u = x(1−x) y², v = 0 : K = ½ · (1/30) · (1/5) = 1/300.

    (Sur le champ de Taylor–Green, la quadrature est EXACTE — erreur nulle — car l'intégrande est
    trigonométrique : on ne verrait pas l'ordre. D'où ce champ polynomial.)
    """
    errors, hs = [], []
    for N in GRIDS:
        h, u, v, Xu, _ = taylor_green(N)
        Yu = (np.arange(N + 2) - 0.5)[None, :] * h
        u = Xu * (1.0 - Xu) * Yu**2
        v = np.zeros_like(v)
        errors.append(abs(kinetic_energy(u, v, h) - 1.0 / 300.0))
        hs.append(h)
    assert abs(observed_order(hs, errors) - 2.0) < 0.1


def test_projection_gives_zero_divergence():
    N = 32
    h = 1.0 / N
    rng = np.random.default_rng(0)
    u, v = new_fields(N)
    u[1:-1, 1:-1] = rng.standard_normal((N - 1, N))   # faces intérieures seulement
    v[1:-1, 1:-1] = rng.standard_normal((N, N - 1))
    apply_bc(u, v, 1.0)
    div = np.zeros((N, N))
    divergence(u, v, h, div)
    div_before = np.max(np.abs(div))
    project(u, v, h)
    divergence(u, v, h, div)
    assert np.max(np.abs(div)) < 1e-12 * div_before


def test_energy_balance_re100():
    """dK/dt (différences centrées sur la série) doit égaler P − ε à O(Δt²) près."""
    res = run(Config(Re=100.0, N=32, t_final=2.0, n_diag=1))
    t, K, P, eps = res["t"], res["K"], res["P"], res["eps"]
    dKdt = (K[2:] - K[:-2]) / (t[2:] - t[:-2])
    residual = np.max(np.abs(dKdt - (P[1:-1] - eps[1:-1])))
    assert residual < 1e-3 * np.max(P)
    assert np.max(res["div_max"]) < 1e-10
