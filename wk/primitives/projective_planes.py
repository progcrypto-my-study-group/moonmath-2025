from .fields import Field, FieldElement


class ProjPoint:
    def __init__(self, x: FieldElement, y: FieldElement, z: FieldElement):
        if not (x.order == y.order == z.order):
            raise ValueError("Coordinates must be of the same field.")
        self.X = x
        self.Y = y
        self.Z = z

    def __eq__(self, other: "ProjPoint"):
        k = self.Z / other.Z
        return self.X == k * other.X and self.Y == k * other.Y

    # Given an order of a field, produce all possible projective planes points
    @staticmethod
    def all_points(order: int):
        coord_list = []
        F = Field(order)
        for x_coord in range(order):
            for y_coord in range(order):
                for z_coord in range(order):
                    coord_list.append(ProjPoint(F(x_coord), F(y_coord), F(z_coord)))
        # TODO: remove points that are repeating of [k*X; k*Y ;k*Z]
        coord_list = coord_list[1:]
        return coord_list
