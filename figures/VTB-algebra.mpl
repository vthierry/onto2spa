#
# Checks VTB algebra operations
#

with(LinearAlgebra):

# Defines the dimension as square
d_:= 3: d := d_^2:

# Defines a literal vector
v := s -> Vector(d, (i) -> cat(s, '_', i)):

# Defines a binding matrix
B := v -> Matrix(d, d, (j, i) -> if 1 <= i - d_ * iquo(j-1, d_) and i - d_ * iquo(j-1, d_) <= d_ then sqrt(d_) * v[1 + irem(i-1, d_) + d_ * irem(j-1, d_)] else 0 fi):
B_ := v -> Matrix(d, d, (j, i) -> if 1 <= j - d_ * iquo(i-1, d_) and j - d_ * iquo(i-1, d_) <= d_ then sqrt(d_) * v[1 + irem(i-1, d_) + d_ * irem(j-1, d_)] else 0 fi):

# Defines the explicit binding formula
b := (y, x) -> Vector(d, (i) -> sqrt(d_) * add(y[k + d_ * irem(i-1, d_)] * x[k + d_ * iquo(i-1, d_)], k = 1 .. d_)):

##b_ := (y, x) -> Vector(d, (i) -> sqrt(d_) * add(y[1 + irem(i-1, d_) + d_ * irem((k + d_ * iquo(i-1, d_))-1, d_)] * x[k + d_ * iquo(i-1, d_)], k = 1 .. d_)):
##simplify(b(y, x) - b_(y, x));


# Defines the explicit transpose formula
sigma := (i) -> 1 + d_ * irem(i-1, d_) + iquo(i-1, d_):
t := (y) -> Vector(d, (i) -> y[sigma(i)]):

# Defines the identity vector
u := Vector(d, (i) -> if i = sigma(i) then 1/sqrt(d_) else 0 fi):

# Define the left mirror matrix
M := Matrix(d, d, (i, j) -> if j = sigma(i) then 1 else 0 fi):

# Define the explicit composition formula
o := (y, x) -> Vector(d, (i) -> sqrt(d_) * add(y[k + d_ * iquo(i-1, d_)] * x[1 + d_ * (k - 1) + irem(i-1, d_)], k = 1 .. d_)):

# Verifies the formulas
x := v('x'): y := v('y'): z := v('z'):

# The explicit binding formula test
ok_binding := evalb(0 = Norm(simplify(b(y, x) - b_(y, x))));
ok_binding := evalb(0 = Norm(simplify(b(y, x) - Multiply(B(y), x))));
ok_binding := evalb(0 = Norm(simplify(B(y) - B_(y))));

# The explicit transpose and identity formula test
ok_transpose := evalb(0 = Norm(simplify(Vector(d, (i) -> sigma(sigma(i)) - i))));
ok_transpose := evalb(0 = Norm(simplify(B(u) - IdentityMatrix(d, d))));
ok_transpose := evalb(0 = Norm(simplify(B(t(y))- Transpose(B(y)))));

# The mirror matrix test
ok_mirror := evalb(0 = Norm(simplify(Multiply(M, Multiply(B(x), y)) - Multiply(B(y), x))));
ok_mirror := evalb(0 = Norm(simplify(Multiply(M, Multiply(M, B(y))) - B(y))));

# The explicit composition formula test
ok_composition := evalb(0 = Norm(simplify(B(o(y, x)) - Multiply(B(y), B(x)))));
ok_composition := evalb(0 = Norm(simplify(t(o(y, x)) - o(t(x), t(y)))));
ok_composition := evalb(0 = Norm(simplify(o(z, o(y, x)) - o(o(z, y), x))));
ok_composition := evalb(0 = Norm(simplify(o(l * x + y, z) - (l * o(x, z) + o(y, z)))));

# Then ...
X := Matrix(d, (i, j) -> cat('x', '_', i, j)):

Multiply(M, y);
