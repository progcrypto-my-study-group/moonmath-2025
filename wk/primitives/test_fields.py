import pytest
from .fields import Field, FieldElement


class TestFieldElement:
    """Tests for Field primitives"""

    @pytest.fixture
    def prime(self):
        """A small prime for fast tests."""
        return 19

    @pytest.fixture
    def zero(self, prime):
        return FieldElement(0, prime)

    @pytest.fixture
    def one(self, prime):
        return FieldElement(1, prime)

    @pytest.fixture
    def five(self, prime):
        return FieldElement(5, prime)

    @pytest.fixture
    def eight(self, prime):
        return FieldElement(8, prime)

    def test_field_creation_valid(self):
        Fp = Field(5)
        value_of_2 = Fp(2)
        assert value_of_2.value == 2
        assert value_of_2.order == 5

    def test_value_bigger_than_order(self):
        Fp = Field(5)
        value_of_10 = Fp(10)
        assert value_of_10.value == 10 % 5
        assert value_of_10.order == 5

    def test_negative_value(self):
        Fp = Field(5)
        value_of_negative_one = Fp(-1)
        assert value_of_negative_one.value == 4
        assert value_of_negative_one.order == 5

    def test_field_order_cannot_be_negative(self):
        with pytest.raises(ValueError, match="Order must be a positive integer"):
            Field(-4)

    def test_field_element_order_cannot_be_negative(self):
        with pytest.raises(ValueError, match="Order must be a positive integer"):
            FieldElement(1, -5)

    # ---------- Addition ----------
    def test_add_same_field(self, five, eight):
        result = five + eight
        expected = FieldElement(13, 19)
        assert result == expected

    def test_add_wraps(self, eight):
        result = eight + FieldElement(15, 19)
        expected = FieldElement(4, 19)  # (8+15) % 19
        assert result == expected

    def test_add_different_fields(self):
        a = FieldElement(1, 19)
        b = FieldElement(1, 23)
        with pytest.raises(
            TypeError, match="Cannot add elements from different fields"
        ):
            a + b

    # ---------- Multiplication ----------
    def test_mul_same_field(self, five, eight):
        result = five * eight
        expected = FieldElement(2, 19)  # 5*8 % 19
        assert result == expected

    def test_mul_different_fields(self):
        a = FieldElement(3, 19)
        b = FieldElement(4, 23)
        with pytest.raises(
            TypeError, match="Cannot multiply elements from different fields"
        ):
            a * b

    # ---------- Subtraction ----------
    def test_sub_same_field(self, eight):
        result = eight - FieldElement(3, 19)
        expected = FieldElement(5, 19)
        assert result == expected

    def test_sub_wraps_negative(self, five):
        result = five - FieldElement(10, 19)
        expected = FieldElement(14, 19)  # (5-10) % 19 = 14
        assert result == expected

    def test_sub_different_fields(self):
        a = FieldElement(1, 19)
        b = FieldElement(2, 23)
        with pytest.raises(
            TypeError, match="Cannot subtract elements from different fields"
        ):
            a - b

    # ---------- Division ----------
    def test_div_same_field(self, five, eight):
        # 8 / 5 = 8 * 5^(19-2) mod 19
        result = eight / five
        expected = FieldElement(13, 19)  # 5^(-1) mod 19 = 4, 8*4 % 19 = 14
        assert result == expected

    def test_div_by_one(self, eight):
        assert eight / FieldElement(1, 19) == eight

    def test_div_by_zero_raises(self, eight):
        # pow(0, 19-2, 19) fails
        with pytest.raises(ValueError):
            eight / FieldElement(0, 19)

    def test_div_different_fields(self):
        a = FieldElement(1, 19)
        b = FieldElement(1, 23)
        with pytest.raises(
            TypeError, match="Cannot divide elements from different fields"
        ):
            a / b

    # ---------- Equality ----------
    def test_equality(self, five):
        assert five == FieldElement(5, 19)
        assert five != FieldElement(6, 19)
        assert five != FieldElement(5, 23)  # different prime
        assert five != 5  # different type

    # ---------- Inverse ----------
    def test_inverse(self, five):
        assert five.inv() == FieldElement(4, 19)
