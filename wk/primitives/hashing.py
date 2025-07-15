import hashlib
import random
from .fields import Field, FieldElement


# Hash to field should be similar to hash to cyclic group.
# H : {0, 1}∗ -> G
def naive_hash_to_field(plaintext: str, generator: FieldElement):
    hasher = hashlib.sha256(plaintext.encode())
    ciphertext = hasher.hexdigest()
    ciphertext_value = int(ciphertext, 16)
    # Calculate the exponential map
    return generator**ciphertext_value


# This does not fully follow the actual pederson hash impl
# which is usually based on a curve. This is simply taking
# random values with a hash for a group.
def make_random_generators(field_order: int):
    F = Field(field_order)
    generators = []
    for i in range(32):
        hasher = hashlib.sha256(f"pederson_style_{i}".encode())
        ciphertext = hasher.hexdigest()
        ciphertext_val = int(ciphertext, 16)
        generators.append(F(ciphertext_val))
    return generators


PEDERSEN_STYLE_GENERATORS = make_random_generators(1031)


def pedersen_style_hash_to_field(plaintext: str):
    hasher = hashlib.sha256(plaintext.encode())
    ciphertext = hasher.hexdigest()
    ciphertext_value = int(ciphertext, 16)
    # Convert into byte array.
    ciphertext_bytes = ciphertext_value.to_bytes(32)
    # We are doing something different than what the book,
    # encoding chunks(bytes) instead of bits.
    # This is also what modern pederson hashes do.
    F = Field(1031)

    # Calculate the hash
    result = F(1)
    for i in range(len(PEDERSEN_STYLE_GENERATORS)):
        result *= PEDERSEN_STYLE_GENERATORS[i] ** ciphertext_bytes[i]
    return result
