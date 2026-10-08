import json, numpy as np
from fm_reg_01 import decompose, theory, u_field, mirror
rel = lambda a, b: abs(a - b) / max(abs(b), 1e-300)
RC, RW = np.array([1., -1., 1.]), np.array([-1., 1., 1.])      # correct axis, spec's axis
cases = {"pure_I": (30e6, 0, 0), "pure_II": (0, 20e6, 0), "pure_III": (0, 0, 15e6), "mixed": (30e6, 20e6, 15e6)}
out, checks = {}, {}
for name, K in cases.items():
    th = theory(*K); c = decompose(*K, RC); w = decompose(*K, RW)
    out[name] = {"theory": th, "R_crack_plane": c, "R_spec_diag(-1,1,1)": w}
    checks[f"{name}: J total matches"] = rel(c["J"], th["J"]) < 1e-4
    checks[f"{name}: J+ = K_I^2/E' (correct R)"] = abs(c["J_plus"] - th["J_plus"]) < 1e-4 * th["J"]
    checks[f"{name}: J- = K_II^2/E'+K_III^2/2mu (correct R)"] = abs(c["J_minus"] - th["J_minus"]) < 1e-4 * th["J"]
    checks[f"{name}: cross term = 0 (correct R)"] = abs(c["cross"]) < 1e-4 * th["J"]
# path independence
K = cases["mixed"]; Js = [decompose(*K, RC, radius=r)["J"] for r in (1e-5, 1e-3, 1e-1)]
checks["path independence (3 radii)"] = max(Js) - min(Js) < 1e-4 * Js[1]
# parity rules from the spec: mode I even, II and III odd under the correct R
for name, K, sign in [("mode I even", (1, 0, 0), 1), ("mode II odd", (0, 1, 0), -1), ("mode III odd", (0, 0, 1), -1)]:
    f = lambda a, b, K=K: u_field(a, b, *K); m = mirror(f, RC)
    x1, x2 = np.array([0.3, -0.2, 0.05]), np.array([0.4, 0.1, -0.7])
    checks[f"parity: {name}"] = np.allclose(m(x1, x2), sign * f(x1, x2))
# the spec's axis fails to separate modes
w = out["mixed"]["R_spec_diag(-1,1,1)"]; th = out["mixed"]["theory"]
checks["spec axis diag(-1,1,1) does NOT give J+ = K_I^2/E'"] = abs(w["J_plus"] - th["J_plus"]) > 1e-2 * th["J"]
for k, v in checks.items(): print(("PASS " if v else "FAIL ") + k)
print(f"{sum(checks.values())}/{len(checks)} checks pass")
m = out["mixed"]; print(f"\nmixed mode (K_I=30, K_II=20, K_III=15 MPa*sqrt(m)), J in N/m:")
for lab in ("theory", "R_crack_plane", "R_spec_diag(-1,1,1)"):
    d = m[lab]; print(f"  {lab:22s} J={d['J']:.4f}  J+={d['J_plus']:.4f}  J-={d['J_minus']:.4f}" + (f"  cross={d['cross']:.2e}" if 'cross' in d else ""))
json.dump({"checks": checks, "results": out}, open("fm_reg_01_results.json", "w"), indent=1)
assert all(checks.values())
