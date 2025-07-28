from .fields import Field


# TODO: Implementing this properly
class Polynomial:
    def __init__(self, coefficients, field_order):
        self.field_order = field_order
        self.field = Field(field_order)
        self.coefficients = [self.field(co) for co in coefficients]
        self.degree = max(len(coefficients) - 1, 0)


class PrimeFieldExtension:
    def __init__(
        self,
        order,
    ):
        pass


def poly_xgcd():
    pass


def poly_mod(poly, modulus):
    pass
