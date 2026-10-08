# FM-REG-01: mirror-symmetry decomposition of the crack-tip J-integral

Run `python verify_fm_reg_01.py`: **21/21 checks pass.** Needs NumPy only.

The model uses exact Williams near-tip fields (linear elastic, plane strain plus antiplane), the mirror operator
(Ma)(x) = R a(Rx), the projections P± = (I ± M)/2, and the J-integral evaluated numerically on circular contours.

| Mixed mode (K_I = 30, K_II = 20, K_III = 15 MPa·√m), J in N/m | J | J₊ | J₋ | cross |
|---|---|---|---|---|
| Theory | 7377.50 | 4095.00 | 3282.50 | 0 |
| R = diag(1, −1, 1), crack plane | 7377.50 | 4095.00 | 3282.50 | −6e-8 |
| R = diag(−1, 1, 1), original spec | 7377.50 | 0.00 | 0.00 | 7377.50 |

**Confirmed from the original spec:**
- the parity rules (mode I even, modes II and III odd);
- the J₊/J₋ decomposition with zero cross-coupling;
- path independence.

**Corrected:** the reflection must be taken across the crack plane. The original axis sends all of J into the
cross term, so the decomposition gives nothing.

**Removed:** the periodic wrap and the "Mark 27" bound, which have no defined role here. The receipt hash in the
original was a placeholder pattern. The receipt in `FM-REG-01_spec_v2.json` is a real SHA-256 over the code and
results.

**Scope:** this covers the analytic singular fields only. Applying it to real parts needs computed K values
(FEM or a handbook) plus fracture toughness K_IC from testing, and it does not predict exact failure times or
atoms.
