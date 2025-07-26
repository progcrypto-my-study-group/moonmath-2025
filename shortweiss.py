from typing import Tuple, Optional, List
import itertools

class ShortWeierstrassCurve:
    """
    Affine Short Weierstrass curve  y^2 = x^3 + a·x + b  (mod p)
    over a prime field F_p with p > 3.
    """

    def __init__(self, p: int, a: int, b: int):
        assert p > 3 and (4 * a ** 3 + 27 * b ** 2) % p != 0, "Not an elliptic curve"
        self.p = p
        self.a = a % p
        self.b = b % p
        self.O = None            # point at infinity

    # ------------------------------------------------------------------
    # Field helpers
    # ------------------------------------------------------------------
    def _add(self, x: int, y: int) -> int: return (x + y) % self.p
    def _sub(self, x: int, y: int) -> int: return (x - y) % self.p
    def _mul(self, x: int, y: int) -> int: return (x * y) % self.p
    def _inv(self, n: int) -> int:
        """Fermat inverse"""
        return pow(n, -1, self.p)

    # ------------------------------------------------------------------
    # Point validation
    # ------------------------------------------------------------------
    def is_valid(self, P: Optional[Tuple[int, int]]) -> bool:
        if P is self.O:                  # point at infinity
            return True
        x, y = P
        lhs = (y * y) % self.p
        rhs = (x ** 3 + self.a * x + self.b) % self.p
        return lhs == rhs

    # ------------------------------------------------------------------
    # Group law: P + Q
    # ------------------------------------------------------------------
    def add(self, P: Optional[Tuple[int, int]],
            Q: Optional[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        if P is self.O: return Q
        if Q is self.O: return P
        x1, y1 = P
        x2, y2 = Q

        if x1 == x2:
            if y1 == y2:                # doubling
                if y1 == 0:
                    return self.O
                s = (3 * x1 * x1 + self.a) * self._inv(2 * y1) % self.p
            elif (y1 + y2) % self.p == 0:
                return self.O           # P = -Q
            else:
                raise ValueError("Invalid points")
        else:
            s = self._sub(y2, y1) * self._inv(self._sub(x2, x1)) % self.p

        x3 = self._sub(self._sub(self._mul(s, s), x1), x2)
        y3 = self._sub(self._mul(s, self._sub(x1, x3)), y1)
        return (x3, y3)

    # ------------------------------------------------------------------
    # k·P  (double-and-add)
    # ------------------------------------------------------------------
    def mul(self, k: int, P: Optional[Tuple[int, int]]) -> Optional[Tuple[int, int]]:
        if P is self.O or k % self.p == 0:
            return self.O
        k = k % (self.order() if hasattr(self, '_order') else self.p)
        R = self.O
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R

    # ------------------------------------------------------------------
    # Enumerate all affine points (only for small p)
    # ------------------------------------------------------------------
    def affine_points(self) -> List[Tuple[int, int]]:
        pts = []
        for x, y in itertools.product(range(self.p), repeat=2):
            if self.is_valid((x, y)):
                pts.append((x, y))
        return pts

    def __repr__(self):
        return f"E: y^2 = x^3 + {self.a}x + {self.b} over F_{self.p}"


# ----------------------------------------------------------------------
# Demo with the exercise curve
# ----------------------------------------------------------------------
if __name__ == "__main__":
    E = ShortWeierstrassCurve(p=5, a=1, b=1)
    print(E)
    pts = E.affine_points()
    print("Affine points:", sorted(pts))
    print("Cardinality  :", len(pts) + 1)   # +1 for O

    # quick sanity checks
    P = (2, 1)
    Q = (3, 1)
    R = E.add(P, Q)
    print("P + Q =", R)
    print("3*P   =", E.mul(3, P))