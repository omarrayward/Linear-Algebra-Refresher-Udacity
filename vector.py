from math import acos, sqrt, pi
from decimal import Decimal, getcontext

getcontext().prec = 30


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(c) for c in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __iter__(self):
        self.current = 0
        return self

    def next(self):
        if self.current >= len(self.coordinates):
            raise StopIteration
        else:
            current_value = self.coordinates[self.current]
            self.current += 1
            return current_value

    def __len__(self):
        return len(self.coordinates)

    def __getitem__(self, i):
        return self.coordinates[i]

    def __str__(self):
        return 'Vector: {}'.format([round(coord, 3)
                                    for coord in self.coordinates])

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def is_zero(self):
        return set(self.coordinates) == set([Decimal(0)])

    def plus(self, other):
        return Vector(map(sum, zip(self.coordinates, other.coordinates)))

    def minus(self, other):
        return Vector([coords[0] - coords[1]
                       for coords in zip(self.coordinates, other.coordinates)])

    def times_scalar(self, factor):
        return Vector([Decimal(factor) * coord for coord in self.coordinates])

    def magnitude(self):
        return Decimal(sqrt(sum([coord * coord
                                 for coord in self.coordinates])))

    def normalize(self):
        try:
            return self.times_scalar(Decimal('1.0') / self.magnitude())
        except ZeroDivisionError:
            raise Exception('Cannot normalize the zero vector')

    def dot_product(self, other):
        return sum(x * y for x, y in zip(self.coordinates, other.coordinates))

    def get_angle_rad(self, other):
        dot_prod = round(self.normalize().dot_product(other.normalize()), 3)
        return acos(dot_prod)

    def get_angle_deg(self, other):
        degrees_per_rad = 180. / pi
        return degrees_per_rad * self.get_angle_rad(other)

    def is_parallel(self, other):
        return (self.is_zero() or other.is_zero() or
                self.get_angle_rad(other) in [0, pi])

    def is_orthogonal(self, other):
        return round(self.dot_product(other), 3) == 0

    def get_projected_vector(self, other):
        """
        Gets projection of vector v in b
        """
        b_normalized = other.normalize()
        return b_normalized.times_scalar(self.dot_product(b_normalized))

    def get_orthogonal_vector(self, other):
        return self.minus(self.get_projected_vector(other))

    def cross_product(self, other):
        [x1, y1, z1] = self.coordinates
        [x2, y2, z2] = other.coordinates
        x = (y1 * z2) - (y2 * z1)
        y = -((x1 * z2) - (x2 * z1))
        z = (x1 * y2) - (x2 * y1)
        return Vector([x, y, z])

    def area_parallelogram(self, other):
        return self.cross_product(other).magnitude()

    def area_triangle(self, other):
        return self.cross_product(other).magnitude() / 2

# A point represents a location, defined with x, y, z (2, -1)
# A vector represents a change in direction and it is also defined with x, y, z
# This time x, y and z represent the change in each of the coordinates.
# A vectors doesn't have a specific location.


if __name__ == '__main__':
    v = Vector([8.218, -9.341])
    w = Vector([-1.129, 2.111])
    addition = v.plus(w)
    print 'addition: {}'.format(addition)

    v = Vector([7.119, 8.215])
    w = Vector([-8.223, 0.878])
    subtraction = v.minus(w)
    print 'subtraction: {}'.format(subtraction)

    v = Vector([1.671, -1.012, -0.318])
    multiplication = v.times_scalar(7.41)
    print 'multiplication: {}'.format(multiplication)


# Vectors have magnitude and direction. The magnitude of a Vector is calculated
# using the pythagorian formula adapted to as many dimensions as the vector has
# e.g. 3 dimensions = square root of x^^2 + y^^2 + z^^2

# The zero vector is the one that has all its coordinates equal to 0. The zero
# vector has no direction and no normalization
# A unit vector is any vector in which it's magnitude is 1

# The process to find a Unit vector from a given vector is called normalization
# A. Find magnitude
# B. Multiply vector by 1/magnitude of vector

    v = Vector([-0.221, 7.437])
    first_magintude = v.magnitude()
    print 'first_magintude: {}'.format(round(first_magintude, 3))

    v = Vector([8.813, -1.331, -6.247])
    second_magintude = v.magnitude()
    print 'second_magintude: {}'.format(round(second_magintude, 3))

    v = Vector([5.581, -2.136])
    first_normalization = v.normalize()
    print 'first_normailization: {}'.format(first_normalization)

    v = Vector([1.996, 3.108, -4.554])
    second_normalization = v.normalize()
    print 'second_normailization: {}'.format(second_normalization)

# The inner multiplication of vectors is the multiplication of their magnitudes
# times the multiplications of the cos of their angle.
# v1 * v2 = magnitude(v1) * magnitude(v2) * cos(v1-v2).
# cos is bounded by -1 and +1 so:
# v1 * v2 <= magnitude(v1) * magnitude(v2) // Inequality of Cauchy-Schwartz
# if v1 * v2 == magnitude(v1) * magnitude(v2) => Angle is 0deg (same angle)
# if v1 * v2 == - (magnitude(v1) * magnitude(v2)) => Angle is 180deg
# (opposite angle)
# if v1 * v2 == 0 and v1 != 0 and v2 != 0 => Angle is 90deg
# v * v => magnitude(v) * magnitude(v) * cos(0) => magnitude(v) = sqrt(v * v)


# It can also be expressed as:
# v1 * v2 = sum of multiplication of each of its coordinates
# So the angle of 2 vectors:
# arc-cosin of (v1*v2)/(magnitude(v1) * magnitude(2)) => radians
# arc-cosin of (normalize(v1) * normalize(v2))

# http://betterexplained.com/articles/vector-calculus-understanding-the-dot-product/  # NOQA
# Basically the dot product gives a number that applies the directional growth
# of a vector into another one

    v = Vector([7.887, 4.138])
    w = Vector([-8.802, 6.776])
    dot_product = v.dot_product(w)
    print 'first_dot_product: {}'.format(round(dot_product, 3))

    v = Vector([-5.955, -4.904, -1.874])
    w = Vector([-4.496, -8.755, 7.103])
    dot_product = v.dot_product(w)
    print 'second_dot_product: {}'.format(round(dot_product, 3))

    v = Vector([3.183, -7.627])
    w = Vector([-2.668, 5.319])
    angle_rads = v.get_angle_rad(w)
    print 'first_angle_rads: {}'.format(angle_rads)

    v = Vector([7.35, 0.221, 5.188])
    w = Vector([2.751, 8.259, 3.985])
    angle_degrees = v.get_angle_deg(w)
    print 'first_angle_rads: {}'.format(angle_degrees)


# 2 vectors are parallel if one is a scalar multiple of another one
# e.g: v and 2v and 1/2v and Zero vector

# 2 vectors are orthogonal if v1 * v2 is = 0. There are 2 possibilities.
# A vector that is 90 deg with another vector or the Zero vector.

# Zero vector is both parallel and orthogonal to all other vectors
# Zero vector is the only one that is orthogonal to itself

    v = Vector([-7.579, -7.88])
    w = Vector([22.737, 23.64])
    is_parallel = v.is_parallel(w)
    is_orthogonal = v.is_orthogonal(w)

    print '1 parallel: {}, orthogonal: {}'.format(is_parallel, is_orthogonal)

    v = Vector([-2.029, 9.97, 4.172])
    w = Vector([-9.231, -6.639, -7.245])
    is_parallel = v.is_parallel(w)
    is_orthogonal = v.is_orthogonal(w)

    print '2 parallel: {}, orthogonal: {}'.format(is_parallel, is_orthogonal)

    v = Vector([-2.328, -7.284, -1.214])
    w = Vector([-1.821, 1.072, -2.94])
    is_parallel = v.is_parallel(w)
    is_orthogonal = v.is_orthogonal(w)
    print '3 parallel: {}, orthogonal: {}'.format(is_parallel, is_orthogonal)

    v = Vector([2.118, 4.827])
    w = Vector([0, 0])
    is_parallel = v.is_parallel(w)
    is_orthogonal = v.is_orthogonal(w)

    print '4 parallel: {}, orthogonal: {}'.format(is_parallel, is_orthogonal)

# Orthogonality is a tool for decomposing objects into combinations of simpler
# objects in a structured way

# Proyecting a vector into another vector:
# https://www.udacity.com/course/viewer#!/c-ud953/l-4374471116/m-4583493277
# v = v (parallel to base vector) + v (orthogonal to base vector)
# v parallel is a cathetus, v orthogonal is the other cathetus and v is the
# hipotenuse
# cos angle = magnitude(v parallel) / magnitude(v) =>
# magnitude(v parallel) = cos angle * magnitude(v)

# (in this case v2 is the base vector)
# from before (v1 * v2) / (magnitude(v1) * magnitude(v2)) = cos(v1-v2)
# Substituting cos in both functions :

# magnitude(v parallel) = magnitude(v) * (v1 * v2) / (mag(v1) * mag(v2)) =>
# magnitude(v1 parallel) = (v1 * v2) /  magnitude(v2)) =>
# magnitude(v1 parallel) = v1 * normalize(v2)

# v1 parallel has the same direction as v2
# magnitude(v1 parallel) * normalize(v2) = v1 parallel =>
# (v1 * normalize(v2)) * normalize(v2) = v1 parallel

    v = Vector([3.039, 1.879])
    w = Vector([0.825, 2.036])
    projected_vector = v.get_projected_vector(w)

    print 'projected vector is: {}'.format(projected_vector)

    v = Vector([-9.88, -3.264, -8.159])
    w = Vector([-2.155, -9.353, -9.473])
    orthogonal_vector = v.get_orthogonal_vector(w)

    print 'orthogonal vector is: {}'.format(orthogonal_vector)

    v = Vector([3.009, -6.172, 3.692, -2.51])
    w = Vector([6.404, -9.144, 2.759, 8.718])
    projected_vector = v.get_projected_vector(w)
    orthogonal_vector = v.get_orthogonal_vector(w)

    print 'second projected vector is: {}'.format(projected_vector)

    print 'second orthogonal vector is: {}'.format(orthogonal_vector)


# Cross product: It's only possible in 3-dimension vectors and it gives
# another vector.

# Cross product of v and W is a vector that is orthogonal to both and its
# magintute = magnitude(v) * magnitude(w) * sin(angle v-w)

# Cross product of 2 parallel vectors is the 0 vector
# If v or w are the 0 vector the cross product will be the 0 vector
# Cross product is anticommutative, which means that the order matters in this
# case the vectors have opposite direction

# cross product(v, w) = - cross product(v, w)

# the area of the parallelogram created with v and w is the magnitude of the
# cross product of v and w

# And the area of the triangle created by those 2 vectors v and w is half of
# the area of the parallelogram


# http://betterexplained.com/articles/cross-product/

    v1 = Vector([8.462, 7.893, -8.187])
    w1 = Vector([6.984, -5.975, 4.778])

    v2 = Vector([-8.987, -9.838, 5.031])
    w2 = Vector([-4.268, -1.861, -8.866])

    v3 = Vector([1.5, 9.547, 3.691])
    w3 = Vector([-6.007, 0.124, 5.772])

    first_cross_product = v1.cross_product(w1)
    print 'cross product is: {}'.format(first_cross_product)

    area_parallelogram = v2.area_parallelogram(w2)
    print 'area parallelogram is: {}'.format(round(area_parallelogram, 3))

    area_triangle = v3.area_triangle(w3)
    print 'area triangle is: {}'.format(round(area_triangle, 3))
