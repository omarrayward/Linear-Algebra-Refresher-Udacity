from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


class Plane(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

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

            initial_index = Plane.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c / initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
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
            initial_index = Plane.first_nonzero_index(n)
            terms = [write_coefficient(
                     n[i], is_initial_term=(i == initial_index)) +
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

    def is_parallel(self, plane2):
        return self.normal_vector.is_parallel(plane2.normal_vector)

    def __eq__(self, plane2):
        if self.normal_vector.is_zero():
            if not plane2.normal_vector.is_zero():
                return False

            diff = self.constant_term - plane2.constant_term
            return MyDecimal(diff).is_near_zero()

        elif plane2.normal_vector.is_zero():
            return False

        if not self.is_parallel(plane2):
            return False

        basepoint_difference = self.basepoint.minus(plane2.basepoint)
        return basepoint_difference.is_orthogonal(self.normal_vector)

    def __iter__(self):
        self.current = 0
        return self

    def next(self):
        if self.current >= len(self.normal_vector):
            raise StopIteration
        else:
            current_value = self.normal_vector[self.current]
            self.current += 1
            return current_value

    def __len__(self):
        return len(self.normal_vector)

    def __getitem__(self, i):
        return self.normal_vector[i]

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)


if __name__ == '__main__':
    # first system of planes:
    # -0.412x + 3.806y + 0.728z = -3.46
    # 1.03x - 9.515y - 1.82z = 8.65

    plane1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
    plane2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)

    print '1 is parallel: {}'.format(plane1.is_parallel(plane2))
    print '1 is equal: {}'.format(plane1 == plane2)

    # second system of planes:
    # 2.611x + 5.528y + 0.283z = 4.6
    # 7.715x + 8.306y + 5.342z = 3.76

    plane3 = Plane(Vector([2.611, 5.528, 0.283]), 4.6)
    plane4 = Plane(Vector([7.715, 8.306, 5.342]), 3.76)

    print '2 is parallel: {}'.format(plane3.is_parallel(plane4))
    print '2 is equal: {}'.format(plane3 == plane4)

    # third system of planes:
    # -7.926x + 8.625y - 7.212z = -7.952
    # -2.642x + 2.875y - 2.404z = -2.443

    plane5 = Plane(Vector([-7.926, 8.625, -7.212]), -7.952)
    plane6 = Plane(Vector([-2.642, 2.875, -2.404]), -2.443)

    print '3 is parallel: {}'.format(plane5.is_parallel(plane6))
    print '3 is equal: {}'.format(plane5 == plane6)
