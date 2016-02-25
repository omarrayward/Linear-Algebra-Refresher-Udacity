### Vectors


A point represents a location, defined with x, y, z (2, -1, 0).
A vector represents a change in direction and it is also defined with x, y, z.
This time x, y and z represent the change in each of the coordinates.
A vectors doesn't have a specific location.


Vectors have magnitude and direction. The magnitude of a vector is calculated
using the pythagorian formula adapted to as many dimensions as the vector has.
e.g. 3 dimensions = `square-root(x**2 + y**2 + z**2)`


The Zero vector is the one that has all its coordinates equal to 0. The Zero
vector has no direction and no normalization.

A unit vector is any vector in which it's magnitude is 1.

The process of finding a unit vector from a given vector is called normalization.
There are 2 steps involved when normilizing a vector:
- A. Find magnitude.
- B. Multiply vector by 1/magnitude of vector.

#### Dot product
The inner multiplication or dot product multiplication of vectors is a way of
multiplying 2 vectors and getting a scalar (a couple of references [here](http://betterexplained.com/articles/vector-calculus-understanding-the-dot-product/) and [here](http://mathworld.wolfram.com/InnerProduct.html)).
The dot product represents a scalar a number that applies the directional growth
of a vector into another one.

It is calculated by multiplying their magnitudes times the multiplications of the
cos of their angle.

`v1 * v2 = magnitude(v1) * magnitude(v2) * cos(v1-v2)`

cos is bounded by -1 and +1 so:

`v1 * v2 <= magnitude(v1) * magnitude(v2)` // Inequality of Cauchy-Schwartz.

if `v1 * v2 == magnitude(v1) * magnitude(v2)` => Angle is 0deg (same angle).

if `v1 * v2 == - (magnitude(v1) * magnitude(v2))` => Angle is 180deg
(opposite angle).

if `v1 * v2 == 0 and v1 != 0 and v2 != 0` => Angle is 90deg

`v * v => magnitude(v) * magnitude(v) * cos(0)` => `magnitude(v) = sqrt(v * v)`

It can also be expressed as:

`v1 * v2 = sum of multiplication of each of its coordinates`

To find out the angle of 2 vectors:

`arc-cosin ((v1*v2)/(magnitude(v1) * magnitude(2)))` => radians.

`arc-cosin((normalize(v1) * normalize(v2)))`.

#### Parallel and orthogonal vectors
2 vectors are parallel if one is a scalar multiple of another one.
e.g: v and 2v and 1/2v and Zero vector

2 vectors are orthogonal if v1 * v2 is = 0. There are 2 possibilities.
A vector that is 90 deg with another vector or the Zero vector.

Zero vector is both parallel and orthogonal to all other vectors.
Zero vector is the only one that is orthogonal to itself.

Orthogonality is a tool for decomposing objects into combinations of simpler
objects in a structured way.

#### Proyecting a vector into another vector
[Proyecting a vector into another vector](https://www.udacity.com/course/viewer#!/c-ud953/l-4374471116/m-4583493277):

`v1 = v1-parallel (to base vector) + v1-orthogonal (to base vector)`

`v1-parallel` is a cathetus, `v2-orthogonal` is the other cathetus and `v1` is the
hipotenuse.


`(*)` `cos angle = magnitude(v1-parallel) / magnitude(v1)` =>

`magnitude(v1-parallel) = cos angle * magnitude(v1)`.


Being v2 the base vector:

`(v1 * v2) / (magnitude(v1) * magnitude(v2)) = cos(v1-v2)`.

Substituting cos from `(*)`:

`magnitude(v1-parallel) = magnitude(v1) * (v1 * v2) / (magnitude(v1) * magnitude(v2))` =>

(eliminating `magnitude(v1)` in right equation) =>

`magnitude(v1-parallel) = (v1 * v2) /  magnitude(v2))` =>

`magnitude(v1-parallel) = v1 * normalize(v2)`

v1 parallel has the same direction as v2.

`magnitude(v1-parallel) * normalize(v2) = v1-parallel` =>

`(v1 * normalize(v2)) * normalize(v2) = v1-parallel`


#### Cross product
Cross product: It's only possible in 3-dimension vectors and it gives another
vector.

Cross product of v and w is a vector that is orthogonal to both and its
`magintute = magnitude(v) * magnitude(w) * sin(angle v-w)`.

Cross product of 2 parallel vectors is the 0 vector.
If v or w are the 0 vector the cross product will be the 0 vector.
Cross product is anticommutative, which means that the order matters, in this
case changing the order results on vectors having opposite directions.

`cross product(v, w) = - cross product(v, w)`

The area of the parallelogram created with v and w is the magnitude of the
cross product of v and w.

The area of the triangle created by those 2 vectors v and w is half of
the area of the parallelogram

For more information in cross products go [here](http://betterexplained.com/articles/cross-product/)
