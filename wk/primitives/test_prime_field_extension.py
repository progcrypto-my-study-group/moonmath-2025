# import pytest
# from itertools import product
# from .fields import Field

# # ------------------------------------------------------------------
# # 1.  Fixtures (tiny reusable fields)
# # ------------------------------------------------------------------
# @pytest.fixture
# def F3():
#     return Field(3)

# @pytest.fixture
# def F5():
#     return Field(5)

# # ------------------------------------------------------------------
# # 2.  Correct number of elements |𝔽ₚᵐ| = pᵐ
# # ------------------------------------------------------------------
# @pytest.mark.parametrize("p,m", [(3, 2), (5, 2), (7, 3)])
# def test_field_size(p, m):
#     P = Poly.irreducible(p, m)
#     assert len(list(Fpm.elements(P))) == p ** m

# # ------------------------------------------------------------------
# # 3.  Field axioms
# # ------------------------------------------------------------------
# @pytest.mark.parametrize("p,m", [(3, 2), (5, 2)])
# def test_field_axioms(p, m):
#     from myfields import Fpm, Poly
#     P = Poly.irreducible(p, m)
#     zero, one = Fpm.zero(P), Fpm.one(P)

#     # additive group
#     for a in Fpm.elements(P):
#         assert a + zero == a
#         assert a + (-a) == zero

#     # multiplicative group
#     for a in Fpm.elements(P):
#         if a == zero:
#             continue
#         assert a * one == a
#         assert a * a.inv() == one

# # ------------------------------------------------------------------
# # 4.  Polynomial evaluation vs direct construction
# # ------------------------------------------------------------------
# def test_polynomial_evaluation(F3):
#     P = Poly([1, 0, 1], F3)  # t² + 1  (irreducible mod 3)
#     t = Fpm(Poly([0, 1], F3), P)  # indeterminate 't'

#     # (t + 2)(2t + 2) = 2t  (see the MoonMath example)
#     lhs = (t + Fpm.from_int(2, P)) * (Fpm.from_int(2, P) * t + Fpm.from_int(2, P))
#     rhs = Fpm.from_int(0, P) * t + Fpm.from_int(2, P)
#     assert lhs == rhs

# # ------------------------------------------------------------------
# # 5.  Inverse via Extended Euclidean Algorithm
# # ------------------------------------------------------------------
# @pytest.mark.parametrize("p,m", [(3, 2), (5, 3)])
# def test_inverse_euclidean(p, m):
#     P = Poly.irreducible(p, m)
#     zero = Fpm.zero(P)
#     for _ in range(50):          # random spot checks
#         a = Fpm.random(P)
#         if a == zero:
#             continue
#         assert a * a.inv() == Fpm.one(P)

# # ------------------------------------------------------------------
# # 6.  Tower sub-field inclusion
# # ------------------------------------------------------------------
# def test_subfield_embedding(F5):
#     P = Poly.irreducible(5, 2)  # 𝔽₅²
#     # 𝔽₅ embeds as degree-0 polynomials
#     for k in range(5):
#         embedded = Fpm.from_int(k, P)
#         assert embedded.coef.degree() == 0
#         assert int(embedded.coef[0]) == k
