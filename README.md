# FM-REG-01: mirror-symmetry decomposition of the crack-tip J-integral

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23244957.svg)](https://doi.org/10.5281/zenodo.23244957)

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


## Prior work
The mirror (symmetric/antisymmetric) decomposition of the J-integral into mode contributions is established
in the fracture-mechanics literature. This repository is an independent, open, test-verified implementation
of it, not a claim of a new method.

- H. Ishikawa, H. Kitagawa, H. Okamura (1979). *J integral of a mixed mode crack and its application.*
  Proc. 3rd Int. Conf. on Mechanical Behaviour of Materials (ICM3), Cambridge.
- O. Huber, J. Nickel, G. Kuhn (1993). *On the decomposition of the J-integral for 3D crack problems.*
  International Journal of Fracture 64, 339–348. https://doi.org/10.1007/BF00017849

**What this repository adds:** exact Williams-field verification of J₊ = K_I²/E′ and J₋ = K_II²/E′ + K_III²/2μ,
zero cross-coupling, path independence, and a negative test showing that reflecting along the crack direction
instead of across the crack plane sends all of J into the cross term (21/21 checks).

## How to cite
Politzer, B. S. (2026). *FM-REG-01: Mirror-symmetry decomposition of the crack-tip J-integral* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.23244957
