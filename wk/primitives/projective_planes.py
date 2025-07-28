from .fields import Field, FieldElement


class ProjPoint:
    def __init__(self, x: FieldElement, y: FieldElement, z: FieldElement):
        if not (x.order == y.order == z.order):
            raise ValueError("Coordinates must be of the same field.")
        self.X = x
        self.Y = y
        self.Z = z

    def __hash__(self):
        return hash(f"Order{self.X.order}_Coord{self.X}{self.Y}{self.Z}".encode())

    def __eq__(self, other: "ProjPoint"):
        base_self = ProjPoint.base_point(self)
        base_other = ProjPoint.base_point(other)

        return (
            base_self.X == base_other.X
            and base_self.Y == base_other.Y
            and base_self.Z == base_other.Z
        )

    def __repr__(self):
        return f"[{self.X.value}; {self.Y.value}; {self.Z.value}]"

    def __str__(self):
        return self.__repr__()

    # Calculate base point [X; Y; Z] from [k*X; k*Y; k*Z].
    @staticmethod
    def base_point(point: "ProjPoint"):
        F = Field(point.X.order)
        ZERO_IN_F = F(0)
        x_coord = point.X
        y_coord = point.Y
        z_coord = point.Z

        # Attempting to calculate scaling, by looking at every coordinate from
        # the right.
        if x_coord != ZERO_IN_F:
            # If Z is not zero, we handle scaling base on Z.
            if y_coord != ZERO_IN_F:
                y_coord = y_coord / x_coord
            if z_coord != ZERO_IN_F:
                z_coord = z_coord / x_coord
            x_coord = x_coord / x_coord  # should be F(1)
        elif y_coord != ZERO_IN_F:
            # If Z is zero, and Y is not, we handle scaling base on Y.
            if z_coord != ZERO_IN_F:
                z_coord = z_coord / y_coord
            y_coord = y_coord / y_coord  # should be F(1)
        elif z_coord != ZERO_IN_F:
            # If both Z and Y are zero, we simply just get the X value divide itself.
            z_coord = z_coord / z_coord  # should be F(1)
        return ProjPoint(x_coord, y_coord, z_coord)

    # Given an order of a field, produce all possible projective planes points,
    # excluding repeated points of the same base point.
    @staticmethod
    def all_points(order: int):
        coord_list = set()
        F = Field(order)
        for x_coord in range(order):
            for y_coord in range(order):
                for z_coord in range(order):
                    new_point = ProjPoint(F(x_coord), F(y_coord), F(z_coord))
                    new_point = ProjPoint.base_point(new_point)

                    # Remove any 0,0,0
                    if new_point != ProjPoint(F(0), F(0), F(0)):
                        coord_list.add(new_point)
        coord_list = list(coord_list)
        return coord_list
