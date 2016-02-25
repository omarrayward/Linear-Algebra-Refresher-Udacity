from copy import deepcopy
from decimal import Decimal, getcontext
from hyperplane import Hyperplane
from plane import Plane
from vector import Vector

getcontext().prec = 30


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


# System with 3 equations. 3 lines.
# We are looking for intersections of these 3 lines.
# If there are no intersections in 3 or more lines the system is "Inconsistent"

# Three dimensional spaces are defined by linear equations of 3 variables
# Ax + By + Cz = k => [A B C] (dot product) [x y z] = k
# [A B C] => Normal vector
# Changing k shifts the plane but it doesn't change the direction of the plane
# 2 parallels planes do not intersect
# Same rules as with lines apply here.
# If 2 planes are parallel:
# if a vector connection one point on each plane is orthogonal to the planes
# normal vectors => Same plane

# Planes that are not parallel intersect on a line not a point.

# With 2 planes (a system of 2, 3-dimenstion equations) could have:
# - No intersections (parallel planes)
# - Intesect in a line
# - Intersect in a whole plane

# With 3 planes (a system of 3, 3-dimenstion equations), the planes could
# intersect:
# - In a line
# - No solution
# - A plane
# - A single point

# The coefficients of a linear equation are the coordinates of a normal vector
# to the object the equation defines

# We need at least 2 lines in 2 variables to obtain a unique intersection
# We need at least 3 planes in 3 variables to obtain a unique intersection

# In order to solve a system wtih 3 linear equations and 3 variables we need
# a set of "tools" to manipulate its data (they need to be reversible):
# 1. Swap equations
# 2. Multiply each of the sides of the planes that is not 0
# 3. Add a multiple of one equation to another equation

# Gaussanian elimination:
# Clear a variable in successive equations until there is only one equation
# with one variable.
# At that point 3 eq with 3 variables have been converted to a system of 3
# equations, one with 3 variables, one with 2 variables and another one with
# one variable (triangular form)
# The triangular form doesn't represent the same set of planes as the original
# system, but it represents the same result.

# If we are resolving the system and at any point we get the one of the sides
# to be 0 and the side to be not 0 => Inconsistent system, the system has no
# solution

# If the system gets to a point in which we have 0 = 0, that means that the
# last equation that we were soloving was redundant and gives no additional
# information.
# The system in this case has been converted from a 3 plane system to a 2 plane
# system, and now we are not looking for a single point of intersection but
# for a line of intersection

# At this time we need to parametrize the solution set. We need to identify the
# pivot variables (lead variables), which are variables in the leading term,
# when the system is in a triangular form. The other variable is called the
# free variable.
# This free variables will become a parameter in the parametrization. So in
# order to find the values of of the pivotal variables for any given point in
# the line we'll use the free variable (parameter) to compute it

# Then we have to make the coefficient of pivot variable to be 1
# Then we clear the pivot variable of the second equation in the first
# equation leaving the system:
# x + Az = B
# y + Cz = D
# The leading variables have coefficient one and they only appear on one
# equation each. The system is simplified as much as possible
# "reduced row-echelon form" or rref
# [x y z ] = [(B - AZ) (D - Cz) z ] =>
# [B D 0] + z* [-A -C 1]
# Base point [B D 0]
# Direction vector [-A -C 1]

# In order for a system to intersect in only one point it needs at least
# 3 equations, 3 planes

# A system is inconsistent if we find 0 = k during Gaussian elimination
# A consistent solution has a unique solution and each variable is a pivotal
# or lead variable
# Number of free variables is equal to the dimension of the solution set


def _get_new_plane(coefficient, plane):
    new_normal_vector = plane.normal_vector.times_scalar(coefficient)
    return Hyperplane(normal_vector=new_normal_vector,
                      constant_term=coefficient * plane.constant_term)


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = ('All planes in the system should '
                                          'live in the same dimension')
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i + 1, p)
                for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        self[row] = _get_new_plane(coefficient, self[row])

    def add_multiple_times_row_to_row(self, coefficient, row_to_add,
                                      row_to_be_added_to):
        recipient_plane = self[row_to_be_added_to]
        new_plane = _get_new_plane(coefficient, self[row_to_add])
        new_normal_vector = \
            recipient_plane.normal_vector.plus(new_plane.normal_vector)
        constant_term = new_plane.constant_term + recipient_plane.constant_term
        self[row_to_be_added_to] = Hyperplane(normal_vector=new_normal_vector,
                                              constant_term=constant_term)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                p.first_nonzero_index(p.normal_vector)
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = system.dimension

        col = 0
        for row in range(num_equations):
            while col < num_variables:
                c = MyDecimal(system[row].normal_vector[col])
                if c.is_near_zero():
                    swap_succeeded = system.did_swap_with_row_below(row, col)
                    if not swap_succeeded:
                        col += 1
                        continue

                system.clear_coefficients_bellow(row, col)
                col += 1
                break

        return system

    def did_swap_with_row_below(self, row, col):
        num_equations = len(self)

        for k in range(row + 1, num_equations):
            coefficient = MyDecimal(self[k].normal_vector[col])
            if not coefficient.is_near_zero():
                self.swap_rows(row, k)
                return True

        return False

    def clear_coefficients_bellow(self, row, col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector[col])

        for row_to_be_added_to in range(row + 1, num_equations):
            n = self[row_to_be_added_to].normal_vector
            gamma = n[col]
            alpha = -gamma / beta
            self.add_multiple_times_row_to_row(alpha, row, row_to_be_added_to)

    def clear_coefficients_above(self, row, col):
        for row_to_be_added_to in range(row)[::-1]:
            n = self[row_to_be_added_to].normal_vector
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, row_to_be_added_to)

    def compute_rref(self):
        tf = self.compute_triangular_form()

        num_equations = len(tf)
        pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()

        for row in range(num_equations)[::-1]:
            pivot_var = pivot_indices[row]
            if pivot_var < 0:
                continue
            tf.scale_row_to_make_coefficient_equal_one(row, pivot_var)
            tf.clear_coefficients_above(row, pivot_var)

        return tf

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)

    def do_gaussian_elimination(self):
        rref = self.compute_rref()

        try:
            rref.raise_excepion_if_contradictory_equation()
            rref.raise_excepion_if_too_few_pivots()
        except Exception as e:
            return e.message

        num_variables = rref.dimension
        solution_coordinates = [rref.planes[i].constant_term
                                for i in range(num_variables)]

        return Vector(solution_coordinates)

    def raise_excepion_if_contradictory_equation(self):
        for plane in self.planes:
            try:
                plane.first_nonzero_index(plane.normal_vector)

            except Exception as e:
                if str(e) == 'No nonzero elements found':
                    constant_term = MyDecimal(plane.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)

                else:
                    raise e

    def raise_excepion_if_too_few_pivots(self):
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indices])
        num_variables = self.dimension

        if num_pivots < num_variables:
            raise Exception(self.INF_SOLUTIONS_MSG)

    def compute_solution(self):
        try:
            return self.do_gaussian_elimination_and_parametrization()

        except Exception as e:
            if str(e) == self.NO_SOLUTIONS_MSG:
                return str(e)
            else:
                raise e

    def do_gaussian_elimination_and_parametrization(self):
        rref = self.compute_rref()
        rref.raise_excepion_if_contradictory_equation()

        direction_vectors = rref.extract_direction_vectors_for_parametrization()  # NOQA
        basepoint = rref.extract_basepoint_for_parametrization()

        return Parametrization(basepoint, direction_vectors)

    def extract_direction_vectors_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()
        free_variable_indices = set(range(num_variables)) - set(pivot_indices)

        direction_vectors = []

        for free_var in free_variable_indices:
            vector_coords = [0] * num_variables
            vector_coords[free_var] = 1
            for index, plane in enumerate(self.planes):
                pivot_var = pivot_indices[index]
                if pivot_var < 0:
                    break
                vector_coords[pivot_var] = -plane.normal_vector[free_var]

            direction_vectors.append(Vector(vector_coords))

        return direction_vectors

    def extract_basepoint_for_parametrization(self):
        num_variables = self.dimension
        pivot_indices = self.indices_of_first_nonzero_terms_in_each_row()

        basepoint_coords = [0] * num_variables

        for index, plane in enumerate(self.planes):
            pivot_var = pivot_indices[index]
            if pivot_var < 0:
                break
            basepoint_coords[pivot_var] = plane.constant_term

        return Vector(basepoint_coords)


p0 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')

s = LinearSystem([p0, p1, p2, p3])
# Print initial system
# print s

p0 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')

s = LinearSystem([p0, p1, p2, p3])
s.swap_rows(0, 1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print 'test case 1 failed'

s.swap_rows(1, 3)
if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
    print 'test case 2 failed'

s.swap_rows(3, 1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print 'test case 3 failed'

s.multiply_coefficient_and_row(1, 0)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print 'test case 4 failed'

s.multiply_coefficient_and_row(-1, 2)
new_s2 = Plane(normal_vector=Vector(['-1', '-1', '1']), constant_term='-3')
if not (s[0] == p1 and
        s[1] == p0 and
        s[2] == new_s2 and
        s[3] == p3):
    print 'test case 5 failed'

s.multiply_coefficient_and_row(10, 1)
new_s1 = Plane(normal_vector=Vector(['10', '10', '10']), constant_term='10')
if not (s[0] == p1 and
        s[1] == new_s1 and
        s[2] == new_s2 and
        s[3] == p3):
    print 'test case 6 failed'

s.add_multiple_times_row_to_row(0, 0, 1)
if not (s[0] == p1 and
        s[1] == new_s1 and
        s[2] == new_s2 and
        s[3] == p3):
    print 'test case 7 failed'

s.add_multiple_times_row_to_row(1, 0, 1)
added_s1 = Plane(normal_vector=Vector(['10', '11', '10']), constant_term='12')
if not (s[0] == p1 and
        s[1] == added_s1 and
        s[2] == new_s2 and
        s[3] == p3):
    print 'test case 8 failed'

s.add_multiple_times_row_to_row(-1, 1, 0)
new_s0 = Plane(normal_vector=Vector(['-10', '-10', '-10']),
               constant_term='-10')
if not (s[0] == new_s0 and
        s[1] == added_s1 and
        s[2] == new_s2 and
        s[3] == p3):
    print 'test case 9 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='2')
s = LinearSystem([p1, p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2):
    print 'test case 1 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='2')
s = LinearSystem([p1, p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == Plane(constant_term='1')):
    print 'test case 2 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')
s = LinearSystem([p1, p2, p3, p4])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2 and
        t[2] == Plane(normal_vector=Vector(['0', '0', '-2']),
                      constant_term='2') and
        t[3] == Plane()):
    print 'test case 3 failed'

p1 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1', '2', '-5']), constant_term='3')
s = LinearSystem([p1, p2, p3])
t = s.compute_triangular_form()
if not (t[0] == Plane(normal_vector=Vector(['1', '-1', '1']),
                      constant_term='2') and
        t[1] == Plane(normal_vector=Vector(['0', '1', '1']),
                      constant_term='1') and
        t[2] == Plane(normal_vector=Vector(['0', '0', '-9']),
                      constant_term='-2')):
    print 'test case 4 failed'


# ***************

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='2')
s = LinearSystem([p1, p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']),
                      constant_term='-1') and
        r[1] == p2):
    print 'test case 1 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='2')
s = LinearSystem([p1, p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print 'test case 2 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')
s = LinearSystem([p1, p2, p3, p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']),
                      constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0', '0', '-2']),
                      constant_term='2') and
        r[3] == Plane()):
    print 'test case 3 failed'

p1 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1', '2', '-5']), constant_term='3')
s = LinearSystem([p1, p2, p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']),
                      constant_term=Decimal('23') / Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0', '1', '0']),
                      constant_term=Decimal('7') / Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0', '0', '1']),
                      constant_term=Decimal('2') / Decimal('9'))):
    print 'test case 4 failed'


# first system
p1 = Plane(Vector([5.862, 1.178, -10.366]), -8.15)
p2 = Plane(Vector([-2.931, -0.589, 5.183]), -4.075)
system1 = LinearSystem([p1, p2])
print 'first system: {}'.format(system1.do_gaussian_elimination())


# # second system
p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)
system2 = LinearSystem([p1, p2, p3])
print 'second system: {}'.format(system2.do_gaussian_elimination())


# third system
p1 = Plane(Vector([5.262, 2.739, -9.878]), -3.441)
p2 = Plane(Vector([5.111, 6.358, 7.638]), -2.152)
p3 = Plane(Vector([2.016, -9.924, -1.367]), -9.278)
p4 = Plane(Vector([2.167, -13.543, -18.883]), -10.567)
system3 = LinearSystem([p1, p2, p3, p4])
print 'thrid system: {} '.format(system3.do_gaussian_elimination())


class Parametrization(object):

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = (
        'The basepoint and direction vectors should all live in the same '
        'dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)

    def __str__(self):

        output = ''
        for coord in range(self.dimension):
            output += 'x_{} = {} '.format(coord + 1,
                                          round(self.basepoint[coord], 3))
            for free_var, vector in enumerate(self.direction_vectors):
                output += '+ {} t_{}'.format(round(vector[coord], 3),
                                             free_var + 1)
            output += '\n'
        return output


p1 = Plane(normal_vector=Vector([0.786, 0.786, 0.588]), constant_term=-0.714)
p2 = Plane(normal_vector=Vector([-0.131, -0.131, 0.244]), constant_term=0.319)

system = LinearSystem([p1, p2])
print system.compute_solution()


p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)

system = LinearSystem([p1, p2, p3])
print system.compute_solution()

p1 = Plane(Vector([0.935, 1.76, -9.365]), -9.955)
p2 = Plane(Vector([0.187, 0.352, -1.873]), -1.991)
p3 = Plane(Vector([0.374, 0.704, -3.746]), -3.982)
p4 = Plane(Vector([-0.561, -1.056, 5.619]), 5.973)

print system.compute_solution()

# Notes
# It says that parametrization code is in the instructor notes, but it isn't
# https://www.udacity.com/course/viewer#!/c-ud953/l-4624329808/e-4897268655/m-4897268656  # NOQA
# This is were it redirects
# https://storage.googleapis.com/supplemental_media/udacityu/4897268656/linsys.py

# Problem with system with 4 lines
# Document with all the Vectors and planes to copy paste and don't have
# formatting bugs

# https://www.udacity.com/course/viewer#!/c-ud953/l-4624329808/e-4897268655/m-4897268656  # NOQA
# wrong answer
# the second equation of first system should be
# -0.131x -0.131y + 0.244z = 0.319
# To be seen here
# https://www.udacity.com/course/viewer#!/c-ud953/l-4624329808/e-4897268655/m-4897268657

# hyperplane in n dimensions is defined as a single linear equation in
# variables
# A hyperplane in n dimensions is an (n-1) dimensional object. // I don't
# understand this concept yet
# A 2 dimenstional plane in 3 dimensions. Single linear equation with 3
# dimensions
# Gaussian elimination works as well


p1 = Hyperplane(normal_vector=Vector([0.786, 0.786]), constant_term=0.786)
p2 = Hyperplane(normal_vector=Vector([-0.131, -0.131]), constant_term=-0.131)

system = LinearSystem([p1, p2])
print system.compute_solution()


p1 = Hyperplane(normal_vector=Vector([2.102, 7.489, -0.786]),
                constant_term=-5.113)
p2 = Hyperplane(normal_vector=Vector([-1.131, 8.318, -1.209]),
                constant_term=-6.775)
p3 = Hyperplane(normal_vector=Vector([9.015, 5.873, -1.105]),
                constant_term=-0.831)

system = LinearSystem([p1, p2, p3])
print system.compute_solution()

p1 = Hyperplane(normal_vector=Vector([0.786, 0.786, 8.123, 1.111, -8.363]),
                constant_term=-9.955)
p2 = Hyperplane(normal_vector=Vector([0.131, -0.131, 7.05, -2.813, 1.19]),
                constant_term=-1.991)
p3 = Hyperplane(normal_vector=Vector([9.015, -5.873, -1.105, 2.013, -2.802]),
                constant_term=-3.982)

system = LinearSystem([p1, p2, p3])
print system.compute_solution()
