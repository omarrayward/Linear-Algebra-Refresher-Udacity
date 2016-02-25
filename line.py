from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps
########################################################################

# Intersections.

# Flat objects such as lines and planes are defined by linear equations:
# - Can add and subtract vars and constants
# - Can multiply a var by a constant
# e.g x + 2y = 1 // y/2 - 2z = x   => Linear equations
# x**2 - 1 = y  // y/x = 3 => Non linear equations

# Finding intersections is basically solving systems of equations.
# Give a set of variables with a set of contraints a system of equations
# is created, and to solve it we need to find the intersection of these
# equations.
# We are basicaly converting a real world problem into a geometrical problem.

# The most basic problem is finding intersection within lines of 2 dimensions
# A line is defined with its basepoint (x0) and the direction of the vector(v)

# x(t) = x0 + times_scalar(v, t) => x(t) = x0 + t * v (parametrization of line)
# Parametrization of a line is the act of conveting a system of equations into
# a set of one dimension equations that all take the same parameters in this
# case "t" (parameter)
# we are defining the line with the parameter t

# Note that x0 is a point, but it's treated as a vector in order to computer
# other points (x(t))

# A given line has infinetely many base points (basically each of its points)
# A given line has infinetely many directions (we can multiply the direction
# vector by any postive scalar and it will give the same direction)

# e.g x(t) = [1, -2] + t*[3, 1]
# for t = 0 => x = [1, -2] => Base point
# for t = 1 => x = [1, -2] + [3, 1] = [3, -1] => Another point on the line

# The line y = mx + b => Base point (0, b)  // m is slope // Direction vector [1, m] # NOQA
# This representation of lines is incomplete since we don't have the posibility
# to represent vertical lines

# A more generic way of representing a line would be Ax + By = k (A, B not both 0) # NOQA
# k => constant term
# Base point => (0, k/B) and (k/A, 0)
# if k = 0
# [A, B] (dot product) [x, y] = 0 => [x, y] is orthogonal to [A, B] =>
# [B, -A] is orthogonal to [A, B] (that's how we figure out the direction vector) # NOQA
# Direction vector => [B, -A]
# Ax + By = k => Gives us an implicit representation of the line such that for
# a set of points (x, y)
# [A, B] (dot product) [x, y] = k
# A line is defined by its normal vector [A, B] and its constant term k
# An orthogonal or normal vector to the line is the one represented by [A, B]

# 2 lines are parallel if their normal vectors are parallel
# if 2 lines are not parallel they'll have an intersection in a single point
# If 2 lines are parallel there are 2 options:
# - They never cross
# - They are the same line, which means that they are conincidental and
# intersect in infinite points

# So 2 lines could have:
# - 1 point of intersection
# - 0 points of intersection
# - infinite points of intersection

# - First check if 2 lines are parallel by checking their normal vecotrs
# - If normal vectors are parallel, then check if they are the same vector
# - Make a line connecting a point of line 1 and line 2 and if that third
# line is orthogonal to either normal vectors for line 1 and 2 => line 1 and 2
# are coincident or the same

# - If we have 2 lines that are not parallel, to know how they intersect:
# Ax + By = k1
# Cx + Dy = k2
# Either A or C need to not be 0, b/c if not the lines would both be horizontal
# and these lines intersect with each other
# if A = 0 => swap A by C, B by D and k1 by k2 => A will never be 0
# (AD - BC) is never 0 b/c that would mean that both lines are parallel
# y = (-Ck1 + Ak2) / (AD - BC)
# x = (Dk1 - Bk2) / (AD - BC)


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0'] * self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]
            basepoint_coords[initial_index] = c / initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Line.first_nonzero_index(n)
            terms = [write_coefficient(n[i],
                                       is_initial_term=(i == initial_index)) +
                     'x_{}'.format(i + 1)
                     for i in range(self.dimension)
                     if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)

    def is_parallel(self, line2):
        return self.normal_vector.is_parallel(line2.normal_vector)

    def __eq__(self, line2):
        if self.normal_vector.is_zero():
            if not line2.normal_vector.is_zero():
                return False

            diff = self.constant_term - line2.constant_term
            return MyDecimal(diff).is_near_zero()

        elif line2.normal_vector.is_zero():
            return False

        if not self.is_parallel(line2):
            return False

        basepoint_difference = self.basepoint.minus(line2.basepoint)
        return basepoint_difference.is_orthogonal(self.normal_vector)

    def intersection(self, line2):

        a, b = self.normal_vector.coordinates
        c, d = line2.normal_vector.coordinates
        k1 = self.constant_term
        k2 = line2.constant_term
        denom = ((a * d) - (b * c))

        if MyDecimal(denom).is_near_zero():
            if self == line2:
                return self
            else:
                return None

        one_over_denom = Decimal('1') / ((a * d) - (b * c))
        x_num = (d * k1 - b * k2)
        y_num = (-c * k1 + a * k2)

        return Vector([x_num, y_num]).times_scalar(one_over_denom)


# first system
# 4.046x + 2.836y = 1.21
# 10.115x + 7.09y = 3.025

line1 = Line(Vector([4.046, 2.836]), 1.21)
line2 = Line(Vector([10.115, 7.09]), 3.025)

print 'first system instersects in: {}'.format(line1.intersection(line2))


# second system
# 7.204x + 3.182y = 8.68
# 8.172x + 4.114y = 9.883

line3 = Line(Vector([7.204, 3.182]), 8.68)
line4 = Line(Vector([8.172, 4.114]), 9.883)

print 'second system instersects in: {}'.format(line3.intersection(line4))

# third system
# 1.182x + 5.562y = 6.744
# 1.773x + 8.343y = 9.525

line5 = Line(Vector([1.182, 5.562]), 6.744)
line6 = Line(Vector([1.773, 8.343]), 9.525)

print 'third system instersects in: {}'.format(line5.intersection(line6))
