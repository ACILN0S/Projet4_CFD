"""Paramètres d'une simulation de cavité entraînée (commun aux méthodes M1 et M2).

Adimensionnement (CLAUDE.md §1) : L = U = ρ = 1, donc
    ν = 1 / Re,   temps en unités de L/U,   vitesse du couvercle U_t(t) = 1 − exp(−t).

Tous les nombres « réglables » sont regroupés ici, dans un objet `Config` immuable
(on en crée un nouveau pour chaque run, on ne le modifie jamais).
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Config:
    """Configuration d'un run.

    Re         : nombre de Reynolds U L / ν
    N          : nombre de cellules par direction (pas h = 1/N)
    scheme     : schéma d'advection, "centered" (ordre 2) ou "upwind1" (ordre 1, contre-exemple)
    cfl        : nombre de Courant visé, Δt ≤ cfl · h / U  (U = 1 : vitesse max = celle du couvercle)
    fo         : nombre de Fourier visé, Δt ≤ fo · h² / ν
    t_final    : temps final (unités L/U)
    n_diag     : on enregistre K, P, ε… tous les n_diag pas de temps
    tol_steady : arrêt anticipé si ‖u^{n+1} − u^n‖∞ / Δt < tol_steady (régime stationnaire)
    """
    Re: float = 100.0
    N: int = 32
    scheme: str = "centered"
    cfl: float = 0.5
    fo: float = 0.25
    t_final: float = 10.0
    n_diag: int = 10
    tol_steady: float = 0.0          # 0 ⇒ pas d'arrêt anticipé : on va jusqu'à t_final

    @property
    def nu(self):
        """Viscosité cinématique adimensionnée ν = 1/Re."""
        return 1.0 / self.Re

    @property
    def h(self):
        """Pas d'espace (grille uniforme, carré unité)."""
        return 1.0 / self.N

    @property
    def dt(self):
        """Pas de temps Δt = min(CFL · h / U, Fo · h² / ν).

        - CFL = 0.5 : RK3-SSP est stable pour l'advection centrée jusqu'à CFL ≈ √3 en 1D ;
          en 2D les deux directions s'ajoutent, 0.5 garde une marge confortable.
        - Fo = 0.25 : limite classique de la diffusion explicite 2D (Euler) ; RK3 la tolère
          jusqu'à ≈ 0.6, on reste donc prudent. Ce terme ne limite qu'à bas Re.
        """
        return min(self.cfl * self.h / 1.0, self.fo * self.h**2 / self.nu)


def lid_velocity(t):
    """Vitesse du couvercle U_t(t) = 1 − exp(−t)  (rampe de l'énoncé, τ = L/U = 1)."""
    return 1.0 - np.exp(-t)
