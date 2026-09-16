"""
ZF22 (Zhang & Frigaard 2022, JFM 947 A32) validation cases.

Table 1 = dimensional inputs, Table 2 = the dimensionless groups ZF22 derives
from them, Table 3 = breakthrough time and displacement efficiency.

Apparatus (ZF22 Section 2.3, ZF23 Section 2.1):
    r_o = 22.23 mm, r_i = 17.46 mm  ->  d_hat = 2.385 mm, r_a_hat = 19.845 mm
    annulus length 4.8 m, vertical, two Newtonian fluids.

KNOWN DEFECT IN TABLE 1, case 9: Q_hat is printed as 2.38e-4 m^3/s, but
Q/A = 2.38e-4 / 5.947e-4 = 0.400 m/s, contradicting the printed w0 = 0.04 m/s.
Case 10 prints the same Q with w0 = 0.4 m/s, which IS consistent.  Case 9's Re
of 100 (vs case 10's 1000) also requires w0 = 0.04.  So the case-9 flow rate is
a typo for 2.38e-5.  We drive every case from w0, which is unambiguous.
"""

from dataclasses import dataclass

R_O_HAT = 22.23e-3
R_I_HAT = 17.46e-3
D_HAT = 0.5 * (R_O_HAT - R_I_HAT)        # 2.385 mm
R_A_HAT = 0.5 * (R_O_HAT + R_I_HAT)      # 19.845 mm
DELTA_STAR = D_HAT / R_A_HAT
ANNULUS_LENGTH = 4.8


@dataclass(frozen=True)
class ZF22Case:
    case: int
    # ---- Table 1, dimensional ----
    w0: float
    q0: float
    rho1: float
    mu1: float
    rho2: float
    mu2: float
    # ---- Table 2, dimensionless (as printed) ----
    e: float
    Re: float
    m: float
    b: float
    # ---- Table 3 ----
    t_br_d2dga: float
    t_br_3d: float
    eta_e_d2dga: float
    eta_e_3d: float


CASES = [
    #        w0     q0        rho1  mu1    rho2      mu2    e    Re    m    b     tbr2D tbr3D etaE2D etaE3D
    ZF22Case(1, 0.032, 1.90e-5, 1000, 0.004, 885.31, 0.02, 0.8, 20, 0.2, -50, 0.44, 0.33, 0.66, 0.61),
    ZF22Case(2, 0.032, 1.90e-5, 1000, 0.004, 1229.38, 0.02, 0.8, 20, 0.2, 100, 0.95, 0.95, 0.95, 0.95),
    ZF22Case(3, 0.032, 1.90e-5, 1000, 0.004, 1022.94, 0.02, 0.6, 20, 0.2, 10, 0.79, 0.67, 0.92, 0.91),
    ZF22Case(4, 0.008, 4.76e-6, 1000, 0.001, 1143.37, 0.005, 0.6, 20, 0.2, 1000, 0.99, 0.98, 1.00, 0.98),
    ZF22Case(5, 0.032, 1.90e-5, 1000, 0.004, 1229.38, 0.02, 0.4, 20, 0.2, 100, 0.95, 0.93, 0.97, 0.95),
    ZF22Case(6, 0.04, 2.38e-5, 1000, 0.005, 1358.41, 0.001, 0.4, 20, 5.0, 100, 0.93, 0.94, 0.93, 0.96),
    ZF22Case(7, 0.032, 1.90e-5, 1000, 0.004, 1022.94, 0.008, 0.2, 20, 0.5, 10, 0.78, 0.70, 0.90, 0.90),
    ZF22Case(8, 0.08, 4.76e-5, 1000, 0.010, 1143.37, 0.005, 0.2, 20, 2.0, 10, 0.78, 0.70, 0.84, 0.86),
    ZF22Case(9, 0.04, 2.38e-5, 1000, 0.001, 1071.68, 0.005, 0.1, 100, 0.2, 100, 0.97, 0.96, 0.97, 0.95),
    ZF22Case(10, 0.4, 2.38e-4, 1000, 0.001, 1716.83, 0.005, 0.1, 1000, 0.2, 100, 0.97, 0.96, 0.97, 0.95),
]


def build_scaling(case: ZF22Case):
    """ZF22 case -> BF25-convention Scaling object."""
    from d2dga.scaling import HerschelBulkleyFluid, Scaling
    return Scaling(
        HerschelBulkleyFluid.newtonian("displaced", case.rho1, case.mu1),
        HerschelBulkleyFluid.newtonian("displacing", case.rho2, case.mu2),
        r_a_hat_star=R_A_HAT,
        delta_star=DELTA_STAR,
        mean_velocity=case.w0,
    )
