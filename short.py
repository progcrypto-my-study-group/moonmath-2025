from typing import Tuple, Optional, List

class F5:
    """prime field F5"""
    p = 5
    def __init__(self, n: int):
        self.n = n % 5
    def __add__(self, other): return F5(self.n + other.n)
    def __sub__(self, other): return F5(self.n - other.n)
    def __mul__(self, other): return F5(self.n * other.n)
    def __pow__(self, k):     return F5(pow(self.n, k, 5))
    def __neg__(self):        return F5(-self.n)
    def inv(self):
        return F5(pow(self.n, -1, 5))
    def __eq__(self, other):  return self.n == other.n
    def __repr__(self):       return str(self.n)

class ShortWeierstrass:
    """Affine Short Weierstrass curve y^2 = x^3 + a x + b over F5"""
    def __init__(self, a: F5, b: F5):
        self.a, self.b = a, b
        # discriminant check
        delta = F5(4) * a**3 + F5(27) * b**2
        assert delta != F5(0), "Singular curve"

    # ------------------------------------------------------------------
    # point validation
    def on_curve(self, x: F5, y: F5) -> bool:
        return y*y == x**3 + self.a * x + self.b

    # ------------------------------------------------------------------
    # list all affine points
    def affine_points(self) -> List[Tuple[int, int]]:
        pts = []
        for x in range(5):
            for y in range(5):
                if self.on_curve(F5(x), F5(y)):
                    pts.append((x, y))
        return pts

    # ------------------------------------------------------------------
    # pretty string
    def __repr__(self):
        return f"E: y^2 = x^3 + {self.a}x + {self.b} over F5"



F = F5
a = F(2)
b = F(4)

# discriminant check
delta = F(4)*a**3 + F(27)*b**2
print("Δ =", delta, "≠ 0 → curve non-singular")

E = ShortWeierstrass(a, b)
print(E)

# point P = (0,2)
x0, y0 = F(0), F(2)
assert E.on_curve(x0, y0), "Point not on curve"
print("P =", (x0.n, y0.n))

# all affine points
pts = E.affine_points()
print("Affine points:", sorted(pts))
print("Cardinality  :", len(pts) + 1)   # +1 for O