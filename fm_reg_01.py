"""FM-REG-01: mirror-symmetry decomposition of the crack-tip J-integral.

Linear elastic, plane strain + antiplane, exact Williams near-tip fields.
Mirror operator (Ma)(x) = R a(Rx); projections P_pm = (I pm M)/2.
Prediction (Ishikawa-type decomposition), with E' = E/(1-nu^2), mu = E/(2(1+nu)):
    J_total = K_I^2/E' + K_II^2/E' + K_III^2/(2 mu)
    J(P+ u) = K_I^2/E'                     (mode I  : mirror-even)
    J(P- u) = K_II^2/E' + K_III^2/(2 mu)   (mode II, III : mirror-odd)
    cross term J - J+ - J- = 0
These hold ONLY for reflection across the crack plane, R = diag(1,-1,1).
"""
import numpy as np

E, NU = 200e9, 0.3                         # steel-like, Pa
MU = E / (2 * (1 + NU)); KAPPA = 3 - 4 * NU; LAM = E * NU / ((1 + NU) * (1 - 2 * NU)); EP = E / (1 - NU**2)

def u_field(x1, x2, KI, KII, KIII):
    """Williams near-tip displacements; crack on x1<0, x2=0; branch cut at theta=+-pi."""
    r = np.hypot(x1, x2); t = np.arctan2(x2, x1); s = np.sqrt(r / (2 * np.pi)); c2, s2 = np.cos(t / 2), np.sin(t / 2)
    u1 = KI / (2 * MU) * s * c2 * (KAPPA - 1 + 2 * s2**2) + KII / (2 * MU) * s * s2 * (KAPPA + 1 + 2 * c2**2)
    u2 = KI / (2 * MU) * s * s2 * (KAPPA + 1 - 2 * c2**2) - KII / (2 * MU) * s * c2 * (KAPPA - 1 - 2 * s2**2)
    u3 = 2 * KIII / MU * s * s2
    return np.array([u1, u2, u3])

def mirror(field, R):
    """(M a)(x) = R a(R x), acting in the x1-x2 plane (x3 is translation-invariant)."""
    return lambda x1, x2: R[:, None] * field(R[0] * x1, R[1] * x2)

def J_integral(field, radius, n=20000):
    """J = int_Gamma (W n1 - sigma_ij n_j du_i/dx1) ds on a circle, theta in (-pi, pi)."""
    th = -np.pi + (np.arange(n) + 0.5) * 2 * np.pi / n
    x1, x2 = radius * np.cos(th), radius * np.sin(th); h = 1e-7 * radius
    d1 = (field(x1 + h, x2) - field(x1 - h, x2)) / (2 * h)   # du_i/dx1
    d2 = (field(x1, x2 + h) - field(x1, x2 - h)) / (2 * h)   # du_i/dx2
    e11, e22, e12 = d1[0], d2[1], 0.5 * (d1[1] + d2[0])
    tr = e11 + e22
    s11, s22, s12 = LAM * tr + 2 * MU * e11, LAM * tr + 2 * MU * e22, 2 * MU * e12
    s13, s23 = MU * d1[2], MU * d2[2]
    W = 0.5 * (s11 * e11 + s22 * e22 + 2 * s12 * e12) + 0.5 * (s13 * d1[2] + s23 * d2[2])
    n1, n2 = np.cos(th), np.sin(th)
    t1, t2, t3 = s11 * n1 + s12 * n2, s12 * n1 + s22 * n2, s13 * n1 + s23 * n2
    integrand = W * n1 - (t1 * d1[0] + t2 * d1[1] + t3 * d1[2])
    return float(np.sum(integrand) * radius * 2 * np.pi / n)

def decompose(KI, KII, KIII, R, radius=1e-3):
    f = lambda a, b: u_field(a, b, KI, KII, KIII); m = mirror(f, R)
    plus = lambda a, b: 0.5 * (f(a, b) + m(a, b)); minus = lambda a, b: 0.5 * (f(a, b) - m(a, b))
    J, Jp, Jm = J_integral(f, radius), J_integral(plus, radius), J_integral(minus, radius)
    return dict(J=J, J_plus=Jp, J_minus=Jm, cross=J - Jp - Jm)

def theory(KI, KII, KIII):
    return dict(J=KI**2 / EP + KII**2 / EP + KIII**2 / (2 * MU), J_plus=KI**2 / EP, J_minus=KII**2 / EP + KIII**2 / (2 * MU))
