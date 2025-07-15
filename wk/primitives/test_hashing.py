import pytest
from hashlib import sha256
from .fields import Field, FieldElement
from .hashing import naive_hash_to_field, pedersen_style_hash_to_field


class TestNaiveHashingIntoFields:
    @pytest.fixture
    def tiny_field(self):
        """A miniature field of prime order 101 for fast testing."""
        return Field(1031)

    @pytest.fixture
    def generator(self, tiny_field):
        return tiny_field(2)

    def test_output_is_field_element(self, generator):
        """The function must return an element of the same field as the generator."""
        out = naive_hash_to_field("hello", generator)
        assert isinstance(out, FieldElement)
        assert out.order == generator.order

    def test_deterministic(self, generator):
        """Same input → same output."""
        a = naive_hash_to_field("foo", generator)
        b = naive_hash_to_field("foo", generator)
        assert a == b

    def test_different_inputs_different_outputs(self, generator):
        """Different strings should almost certainly map to different elements."""
        seen = set()
        for s in ("a", "b", "c", "x"):
            hashed = naive_hash_to_field(s, generator)
            seen.add(hashed.value)
        assert len(seen) == 4

    def test_empty_string(self, generator):
        """Edge case: empty plaintext."""
        out = naive_hash_to_field("", generator)
        assert isinstance(out, FieldElement)
        # quick sanity check: exponent should be sha256(b'').hex converted to int
        expected_exp = int(sha256(b"").hexdigest(), 16)
        assert out == (generator**expected_exp)


class TestPedersenStyleHashing:
    def test_output_is_same_order(self):
        hash_target = "this is the hash message"
        out = pedersen_style_hash_to_field(hash_target)
        assert out.order == 1031

    def test_deterministic(self):
        hash_target = "some message"
        out_1 = pedersen_style_hash_to_field(hash_target)
        out_2 = pedersen_style_hash_to_field(hash_target)
        assert out_1 == out_2

    def test_diff_input_should_result_in_diff_output(self):
        out_1 = pedersen_style_hash_to_field("message 1")
        out_2 = pedersen_style_hash_to_field("message 2")
        assert out_1 != out_2
