# A field is built from a set and 2 ops, add and multiply
# Here we are only interested in
class FieldElement:
    def __init__(self, value, order):
        if order < 0:
            raise ValueError("Order must be a positive integer.")
        self.order = order
        self.value = value
        if not (0 <= value < order):
            self.value = self.value % order

    def __add__(self, other: "FieldElement"):
        if self.order != other.order:
            raise TypeError("Cannot add elements from different fields.")
        new_value = (self.value + other.value) % self.order
        return self.__class__(new_value, self.order)

    def __mul__(self, other: "FieldElement"):
        if self.order != other.order:
            raise TypeError("Cannot multiply elements from different fields.")
        new_value = (self.value * other.value) % self.order
        return self.__class__(new_value, self.order)

    def __rmul__(self, other: "FieldElement"):
        return self * other

    def __sub__(self, other: "FieldElement"):
        if self.order != other.order:
            raise TypeError("Cannot subtract elements from different fields.")
        new_value = (self.value - other.value) % self.order
        return self.__class__(new_value, self.order)

    def __pow__(self, exponent: "FieldElement"):
        # Using python's pow(base, exp, mod)
        new_value = pow(self.value, exponent, self.order)
        return self.__class__(new_value, self.order)

    def __truediv__(self, other: "FieldElement"):
        if self.order != other.order:
            raise TypeError("Cannot divide elements from different fields.")
        if other.value == 0:
            raise ValueError("Zero division error.")
        # a / b = a * b^(p-2) (mod p) by Fermat's Little Theorem
        # b^(p-2) is the modular multiplicative inverse of b.
        inverse = pow(other.value, self.order - 2, self.order)
        new_value = (self.value * inverse) % self.order
        return self.__class__(new_value, self.order)

    def __eq__(self, other: "FieldElement"):
        # Short circuit check: only check the conditions if its the same type
        return (
            type(self) is type(other)
            and self.value == other.value
            and self.order == other.order
        )

    def inv(self):
        # r^(p-2) (mod p) by Fermat's Little Theorem
        inverse = pow(self.value, self.order - 2, self.order)
        return self.__class__(inverse, self.order)

    def __str__(self):
        return f"FieldElement_{self.order}({self.value})"


# Returns a function that generates the field elements of the field.
# Convenience function for ease of use inspired by sage.
def Field(order: int):
    if order < 0:
        raise ValueError("Order must be a positive integer.")

    def create_element(value: int):
        return FieldElement(value, order)

    return create_element
