import pytest
from typing import Tuple
from .fields import Field
from .projective_planes import ProjPoint


@pytest.fixture
def F2():
    return 2


@pytest.fixture
def F3():
    return 3


@pytest.fixture
def F5():
    return 5


@pytest.fixture
def F7():
    return 7


@pytest.fixture
def F13():
    return 13


# TEST 1: correct number of points
@pytest.mark.parametrize("q,expected", [(2, 7), (3, 13), (5, 31), (7, 57)])
def test_point_count(q: int, expected: int):
    assert len(ProjPoint.all_points(q)) == expected


# TEST 2: scaling invariance
@pytest.mark.parametrize("coords", [(1, 2, 3), (0, 1, 1), (1, 0, 0)])
@pytest.mark.parametrize("q", [5, 7])
def test_scaling(coords: Tuple[int, int, int], q: int):
    x, y, z = coords
    p = ProjPoint(x, y, z, q)
    for k in range(1, q):
        assert p == ProjPoint(k * x % q, k * y % q, k * z % q)


# TEST 3: affine vs infinity split
@pytest.mark.parametrize("q", [3, 5])
def test_affine_infinity_split(q: int):
    pts = ProjPoint.all_points(q)
    affine = [p for p in pts if p.Z.value != 0]
    infinity = [p for p in pts if p.Z.value == 0]
    assert len(affine) == q * q
    assert len(infinity) == q + 1


# TEST 4: line membership (Fano plane q=2)
def test_fano_lines():
    pts = ProjPoint.all_points(2)
    lines = [
        {ProjPoint(1, 0, 0), ProjPoint(0, 1, 0), ProjPoint(1, 1, 0)},  # z=0
        {ProjPoint(0, 0, 1), ProjPoint(0, 1, 0), ProjPoint(0, 1, 1)},  # x=0
        {ProjPoint(0, 0, 1), ProjPoint(1, 0, 0), ProjPoint(1, 0, 1)},  # y=0
        {ProjPoint(0, 0, 1), ProjPoint(1, 1, 0), ProjPoint(1, 1, 1)},
    ]
    for line in lines:
        assert len(line) == 3
        assert all(pt in pts for pt in line)


# TEST 5: no zero point
@pytest.mark.parametrize("q", [2, 3, 5, 7])
def test_no_zero_point(q: int):
    F = Field(q)
    assert ProjPoint(F(0), F(0), F(0)) not in ProjPoint.all_points(q)
