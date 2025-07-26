from typing import Tuple

# A projective point is a 3-tuple (X, Y, Z)
Point = Tuple[int, int, int]

def proj_add(P: Point, Q: Point, a: int) -> Point:
    """
    Projective addition on a short Weierstrass curve
        Y²Z = X³ + aXZ² + bZ³
    using the complete addition law given in the question.
    Parameters
    ----------
    P, Q : (X, Y, Z)
        Projective points over an arbitrary field (int, gmpy2.mpz,
        Fp, Fp2 element, etc.).  All arithmetic is delegated to the
        underlying type.
    a : int (or same type as field elements)
        Curve coefficient a in the Weierstrass equation.
    Returns
    -------
    (X3, Y3, Z3) : Point
        Projective representative of P ⊕ Q.
    """

    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q

    # ------------------------------------------------------------
    # 1.  Handle the identity (point at infinity)
    # ------------------------------------------------------------
    if Z1 == 0:                       # [0:1:0] is encoded as (0,1,0)
        return Q
    if Z2 == 0:
        return P

    # ------------------------------------------------------------
    # 2.  Prepare helper variables
    # ------------------------------------------------------------
    U1 = Y2 * Z1
    U2 = Y1 * Z2
    V1 = X2 * Z1
    V2 = X1 * Z2

    # ------------------------------------------------------------
    # 3.  Case 1: same x-coordinate  (V1 == V2)
    # ------------------------------------------------------------
    if V1 == V2:
        if U1 != U2:                  # P = -Q
            return (0, 1, 0)
        else:                         # P == Q  (doubling)
            if Y1 == 0:               # 2-torsion
                return (0, 1, 0)

            W = a * Z1**2 + 3 * X1**2
            S = Y1 * Z1
            B = X1 * Y1 * S
            H = W**2 - 8 * B
            X3 = 2 * H * S
            Y3 = W * (4 * B - H) - 8 * Y1**2 * S**2
            Z3 = 8 * S**3
            return (X3, Y3, Z3)

    # ------------------------------------------------------------
    # 4.  Case 2: distinct x-coordinates  (general addition)
    # ------------------------------------------------------------
    U = U1 - U2
    V = V1 - V2
    W = Z1 * Z2
    A = U**2 * W - V**3 - 2 * V**2 * V2
    X3 = V * A
    Y3 = U * (V**2 * V2 - A) - V**3 * U2
    Z3 = V**3 * W
    return (X3, Y3, Z3)

# Example with a = 4, b = 20, curve y² = x³ + 4x + 20 over F29
p = 29
def mod(x): return x % p

a = 4
P = (mod(3), mod(10), 1)      # affine (3,10)
Q = (mod(9), mod(10), 1)      # affine (9,10)

R = proj_add(P, Q, a)
print("P+Q =", R)             # (X:Y:Z) in projective form

# Converting back to affine (if Z ≠ 0):
if R[2] != 0:
    inv = pow(R[2], -1, p)
    x = (R[0] * inv) % p
    y = (R[1] * inv) % p
    print("affine coordinates:", (x, y))