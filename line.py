from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


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
