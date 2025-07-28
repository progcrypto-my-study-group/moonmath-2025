import pytest
from typing import Tuple
import logging

from .fields import Field
from .projective_planes import ProjPoint

LOGGER = logging.getLogger(__name__)


# TEST 1: correct number of points
@pytest.mark.parametrize("q,expected", [(2, 7), (3, 13), (5, 31), (7, 57)])
def test_point_count(q: int, expected: int):
    LOGGER.info(ProjPoint.all_points(q))
    assert len(ProjPoint.all_points(q)) == expected


# TEST 2: scaling invariance
@pytest.mark.parametrize("coords", [(1, 2, 3), (0, 1, 1), (1, 0, 0)])
@pytest.mark.parametrize("q", [5, 7])
def test_scaling(coords: Tuple[int, int, int], q: int):
    F = Field(q)
    x, y, z = coords
    p = ProjPoint(F(x), F(y), F(z))
    for k in range(1, q):
        x_coord = F(k * x)
        y_coord = F(k * y)
        z_coord = F(k * z)
        assert p == ProjPoint(x_coord, y_coord, z_coord)


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
    F = Field(2)
    lines = [
        {
            ProjPoint(F(1), F(0), F(0)),
            ProjPoint(F(0), F(1), F(0)),
            ProjPoint(F(1), F(1), F(0)),
        },  # z=0
        {
            ProjPoint(F(0), F(0), F(1)),
            ProjPoint(F(0), F(1), F(0)),
            ProjPoint(F(0), F(1), F(1)),
        },  # x=0
        {
            ProjPoint(F(0), F(0), F(1)),
            ProjPoint(F(1), F(0), F(0)),
            ProjPoint(F(1), F(0), F(1)),
        },  # y=0
        {
            ProjPoint(F(0), F(0), F(1)),
            ProjPoint(F(1), F(1), F(0)),
            ProjPoint(F(1), F(1), F(1)),
        },
    ]
    for line in lines:
        assert len(line) == 3
        assert all(pt in pts for pt in line)


# TEST 5: no zero point
@pytest.mark.parametrize("q", [2, 3, 5, 7])
def test_no_zero_point(q: int):
    F = Field(q)
    assert ProjPoint(F(0), F(0), F(0)) not in ProjPoint.all_points(q)


@pytest.mark.parametrize(
    "q,point,expected",
    [
        (3, (1, 2, 2), (1, 2, 2)),
        (5, (1, 0, 2), (1, 0, 2)),
        (5, (2, 3, 4), (1, 4, 2)),
        (5, (0, 2, 3), (0, 1, 4)),
        (5, (1, 2, 0), (1, 2, 0)),
        (5, (1, 0, 0), (1, 0, 0)),
        (5, (0, 1, 0), (0, 1, 0)),
        (5, (0, 0, 1), (0, 0, 1)),
    ],
)
def test_base_point(
    q: int, point: Tuple[int, int, int], expected: Tuple[int, int, int]
):
    F = Field(q)
    p = ProjPoint(F(point[0]), F(point[1]), F(point[2]))
    base_p = ProjPoint.base_point(p)
    expected_p = ProjPoint(F(expected[0]), F(expected[1]), F(expected[2]))
    assert base_p == expected_p
