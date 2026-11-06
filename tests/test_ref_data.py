"""Contrôles de cohérence physique des données de référence saisies (ref_data/).

But : détecter une erreur de saisie/transcription AVANT de s'en servir pour valider nos solveurs.
On n'a pas l'article de Ghia sous la main : ces tests sont notre garde-fou.
"""

import numpy as np
import pytest

from cavity.ref_data import GHIA_SUSPECT, load_profile, load_vortices

GHIA_RE = [100, 400, 1000, 3200, 5000, 7500, 10000]
ERTURK_RE = [1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 21000]

# Débit à travers une ligne médiane : ∫ u dy = 0 (et ∫ v dx = 0) par conservation de la masse.
# Avec 17 points de Ghia, la règle des trapèzes fait une erreur de quadrature de quelques 1e-3 ;
# une coquille d'un chiffre (ex. -0.86636 au lieu de -0.08636) donne un écart ~ 0.03.
FLUX_TOL_GHIA = 0.02
FLUX_TOL_ERTURK = 0.02


def _flux(coord, values):
    ok = ~np.isnan(values)
    return np.trapezoid(values[ok], coord[ok])


@pytest.mark.parametrize("name, re_list", [
    ("ghia_u", GHIA_RE), ("ghia_v", GHIA_RE), ("erturk_u", ERTURK_RE), ("erturk_v", ERTURK_RE),
])
def test_profiles_have_expected_reynolds(name, re_list):
    assert sorted(load_profile(name)) == re_list


@pytest.mark.parametrize("name", ["ghia_u", "erturk_u"])
def test_u_profile_wall_values(name):
    """u = 0 sur la paroi du bas (y = 0) et u = 1 sous le couvercle (y = 1)."""
    for coord, u in load_profile(name).values():
        assert coord[0] == 0.0 and coord[-1] == 1.0
        assert u[0] == 0.0 and u[-1] == 1.0


@pytest.mark.parametrize("name", ["ghia_v", "erturk_v"])
def test_v_profile_wall_values(name):
    """v = 0 sur les parois latérales x = 0 et x = 1."""
    for coord, v in load_profile(name).values():
        assert coord[0] == 0.0 and coord[-1] == 1.0
        assert v[0] == 0.0 and v[-1] == 0.0


@pytest.mark.parametrize("name, tol", [
    ("ghia_u", FLUX_TOL_GHIA), ("ghia_v", FLUX_TOL_GHIA),
    ("erturk_u", FLUX_TOL_ERTURK), ("erturk_v", FLUX_TOL_ERTURK),
])
def test_mass_conservation_through_centerline(name, tol):
    for reynolds, (coord, vel) in load_profile(name).items():
        assert abs(_flux(coord, vel)) < tol, f"{name} Re={reynolds}"


def test_suspect_ghia_values_break_mass_conservation():
    """Contrôle négatif : la coquille Re=3200 non masquée doit être détectée par le test de débit."""
    coord, u = load_profile("ghia_u", mask_suspect=False)[3200]
    assert abs(_flux(coord, u)) > FLUX_TOL_GHIA


def test_suspect_values_are_masked():
    for name, reynolds, c_bad in GHIA_SUSPECT:
        coord, vals = load_profile(name)[reynolds]
        assert np.isnan(vals[np.isclose(coord, c_bad)]).all()


def test_unknown_profile_raises():
    with pytest.raises(KeyError):
        load_profile("inconnu")


def test_vortices_primary_signs_and_location():
    """Tourbillon principal : psi < 0, omega < 0 (rotation horaire), centre proche du milieu."""
    primaries = [r for r in load_vortices() if r["vortex"] == "primary"]
    assert len(primaries) >= 10
    for r in primaries:
        assert r["psi"] < 0 and r["omega"] < 0
        if r["x"] is not None:
            assert 0.45 < r["x"] < 0.6 and 0.5 < r["y"] < 0.6


def test_vortices_re1000_sources_agree():
    """À Re = 1000, Erturk (601², extrapolé) et Botella & Peyret (spectral) s'accordent à 1e-4 près."""
    rec = {r["source"]: r for r in load_vortices() if r["Re"] == 1000 and r["vortex"] == "primary"}
    assert abs(rec["Erturk2005_extrap"]["psi"] - rec["BotellaPeyret1998"]["psi"]) < 1e-4
    assert abs(rec["Erturk2005"]["psi"] - rec["BotellaPeyret1998"]["psi"]) < 5e-4
