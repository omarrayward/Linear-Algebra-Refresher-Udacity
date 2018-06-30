# Purpose: Combining the information learned from the course linear algebra into an interactive code which can be
# used to test equations easily.

from vector import Vector
from line import Line
from linear_system import LinearSystem
from plane import Plane


def error_check_input(text, lower=-99999999999, upper=99999999999):
    while True:
        try:
            operation = float(input(text))
            if operation >= lower and operation <= upper:
                break
            else:
                print("Only enter one of the specified options")
        except ValueError:
            print("You may only enter a number")
    return operation


# Function creates a vector using user input
def enter_vector(times=3):
    coordinates = []
    for x in range(times):
        coordinates.append(float(input("Enter coordinate: ")))
    return Vector(coordinates)


def enter_vectors(vectors=2, choice=True, dimension=None):
    if choice:
        dimension = int(error_check_input("How many dimensions would you like your vectors to include(2-5): ", 2, 5))
    my_vectors = []
    for i in range(vectors):
        print("Entering Vector #%s:" % str(i + 1))
        v = enter_vector(dimension)
        my_vectors.append(v)
    return my_vectors


def vector():
    print("Vector is a quantity with both magnitude and direction\n"
          "Using vectors, we may perform mathematical operations including: \n"
          "\t1) Vector addition\n"
          "\t2) Vector Subtraction\n "
          "\t3) Multiply the vector by a scalar\n"
          "\t4) Find one vector's magnitude\n"
          "\t5) Normalize a vector\n"
          "\t6) Calculate dot product \n"
          "\t7) Compute angle\n"
          "\t8) Determine whether parallel\n"
          "\t9) Determine whether orthogonal\n"
          "\t10) Determine the projected vector\n"
          "\t11) Determine the orthogonal vector\n"
          "\t12) Determine the cross product\n"
          "\t13) Find the area of the parallelogram formed\n"
          "\t14) Find the area of the triangle formed\n")
    option_to_execute = error_check_input("Pick one of the above operations to perform on vectors: ", 1, 14)
    print("===================================================")
    if option_to_execute == 1:
        v, w = enter_vectors(2)
        addition = v.plus(w)
        print("Together the vectors form: %s" % addition)
    elif option_to_execute == 2:
        v, w = enter_vectors(2)
        subtraction = v.minus(w)
        print("Subtracting the vectors forms: %s" % subtraction)
    elif option_to_execute == 3:
        v = enter_vectors(1)[0]
        scalar = error_check_input("Enter scalar: ")
        multiplied = v.times_scalar(scalar)
        print("Multiplying the vector by scalar returns: %s" % multiplied)
    elif option_to_execute == 4:
        v = enter_vectors(1)[0]
        magnitude = v.magnitude()
        print("The magnitude of the vector is: %s" % magnitude)
    elif option_to_execute == 5:
        v = enter_vectors(1)[0]
        normalized = v.normalize()
        print("The vector normalized is: %s" % normalized)
    elif option_to_execute == 6:
        v, w = enter_vectors(2)
        dot_product = v.dot_product(w)
        print('dot_product: {}'.format(round(dot_product, 3)))
    elif option_to_execute == 7:
        v, w = enter_vectors(2)
        angle_degrees = v.get_angle_deg(w)
        angle_radiants = v.get_angle_rad(w)
        print("The first angle is:")
        print("%s in radiant and %s in degrees " % (angle_degrees, angle_radiants))
    elif option_to_execute == 8:
        v, w = enter_vectors(2)
        is_parallel = v.is_parallel(w)
        if is_parallel:
            print("The vectors are parallel")
        else:
            print("The vectors aren't parallel")
    elif option_to_execute == 9:
        v, w = enter_vectors(2)
        is_orthogonal = v.is_orthogonal(w)
        if is_orthogonal:
            print("The vectors are orthogonal")
        else:
            print("The vectors aren't orthogonal")
    elif option_to_execute == 10:
        v, w = enter_vectors(2)
        projected_vector = v.get_projected_vector(w)
        print("The projected vector is: %s" % projected_vector)
    elif option_to_execute == 11:
        v, w = enter_vectors(2)
        orthogonal_vector = v.get_orthogonal_vector(w)
        print("The orthogonal vector: %s" % orthogonal_vector)
    elif option_to_execute == 12:
        v, w = enter_vectors(2, False, 3)
        cross_product = v.cross_product(w)
        print("The cross product is: %s" % cross_product)
    elif option_to_execute == 13:
        v, w = enter_vectors(2, False, 3)
        area_parallelogram = v.area_parallelogram(w)
        print("The parallelogram formed area is: %s" % area_parallelogram)
    elif option_to_execute == 14:
        v, w = enter_vectors(2, False, 3)
        area_triangle = v.area_triangle(w)
        print("The triangle formed area is: %s" % area_triangle)


def enter_lines(times):
    current_lines = []
    for i in range(times):
        print("Enter the two-dimensional vector below:")
        vec = enter_vector(2)
        base_point = error_check_input("Enter the base point: ")
        current_lines.append(Line(vec, base_point))
    return current_lines


def lines():
    print("A line can be defined by a basepoint and a direction vector. ")
    line1, line2 = enter_lines(2)
    print("The lines intersect in %s" % line1.intersection(line2))


def linear_system():
    print("Linear systems use gaussian elimination to\n"
          "\t1) Determine amount of possible solutions for planes\n"
          "\t2) Solve a system of hyperplanes \n")
    option_to_execute = error_check_input("Pick one of the above: ", 1, 2)
    if option_to_execute == 1:
        total_planes = int(error_check_input("Enter amount of planes(2/5): ", 2, 5))
        myPlanes = enter_plane(total_planes)
        system1 = LinearSystem([i for i in myPlanes])
        print("System intersection: ", system1.do_gaussian_elimination())
    elif option_to_execute == 2:
        total_planes = int(error_check_input("Enter amount of planes(2/5): ", 2, 5))
        myPlanes = enter_plane(total_planes)
        system1 = LinearSystem([i for i in myPlanes])
        print(system1.compute_solution())


def enter_plane(times):
    planes = []
    for i in range(times):
        print("Plane #%s" % str(i + 1))
        print("First, enter a 3 dimensional vector: ")
        vector = enter_vector(3)
        base_point = error_check_input("Now, enter base point: ")
        planes.append(Plane(vector, base_point))
    return planes


def plane():
    print("A plane is similar to line, but in three dimensions\nLet's check whether two planes are parallel/equal: ")
    p1, p2 = enter_plane(2)
    if p1.is_parallel(p2):
        print("The planes are parallel")
    else:
        print("The planes are not parallel")
    if p1 == p2:
        print("The planes are equal")
    else:
        print("The planes are not equal")


print("Welcome to the interactive math problem solving library. In this library you will find solutions to well-known"
      "\nmath problems such as vectors, intersection, linear systems and planes.")
while True:
    print("\n"
          "\t1) Vectors\n"
          "\t2) Intersection of lines\n"
          "\t3) Linear Systems\n"
          "\t4) Planes\n"
          "\t5) Exit Program")
    section_view = error_check_input("Choose the operation you would like to execute from the above: ", 1, 4)
    if section_view == 1:
        vector()
    elif section_view == 2:
        lines()
    elif section_view == 3:
        linear_system()
    elif section_view == 4:
        plane()
    elif section_view == 5:
        break
