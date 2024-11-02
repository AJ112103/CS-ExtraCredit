class Complex:
    def __init__(self, real, img):
        self.real = real
        self.img = img

    def __add__(self, y):
        return Complex(self.real + y.real, self.img + y.img)

    def __sub__(self, y):
        return Complex(self.real - y.real, self.img - y.img)

    def __mul__(self, y):
        return Complex(self.real * y.real - self.img * y.img, self.real * y.img + self.img * y.real)

    def __truediv__(self, y):
        denominator = y.real**2 + y.img**2
        if denominator == 0:
            raise ValueError("Cannot divide by 0")
        return Complex((self.real * y.real + self.img * y.img) / denominator, (self.img * y.real - self.real * y.img) / denominator)

    def __abs__(self):
        return (self.real**2 + self.img**2)**0.5

    def conjugate(self):
        return Complex(self.real, -self.img)

class Vector:
    def __init__(self, field, n, *coordinates):
        if field not in ["real", "complex"]:
            raise ValueError("Field must be real or complex")
        
        self.field = field
        self.n = n

        if len(coordinates) != n:
            raise ValueError(f"There should be {n} coordinates")

        if field == "real":
            for coord in coordinates:
                if not isinstance(coord, (int, float)):
                    raise TypeError("coordinates must be real numbers i.e. integers or floats")
        elif field == "complex":
            for coord in coordinates:
                if not isinstance(coord, Complex):
                    raise TypeError("coordinates must be complex numbers i.e. instances of the class Complex)")

        self.coordinates = list(coordinates)

class Matrix:
    def __init__(self, field, n=None, m=None, *values_or_vectors):
        if field not in ["real", "complex"]:
            raise ValueError("Field must be real or complex")
        
        self.field = field

        if all(isinstance(v, Vector) for v in values_or_vectors):
            vectors = values_or_vectors
            m = len(vectors)
            if m == 0:
                raise ValueError("At least one vector is required to form a matrix")
            n = len(vectors[0].coordinates)
            
            for vec in vectors:
                if vec.field != field:
                    raise TypeError("All vectors must match the specified field")
                if len(vec.coordinates) != n:
                    raise ValueError("All vectors must have the same length")

            self.entries = [[vec.coordinates[i] for vec in vectors] for i in range(n)]
            self.n = n
            self.m = m

        else:
            if n is None or m is None:
                raise ValueError("n and m must be specified")
            if len(values_or_vectors) != n * m:
                raise ValueError(f"There should be {n * m} values")

            if field == "real":
                for val in values_or_vectors:
                    if not isinstance(val, (int, float)):
                        raise TypeError("values must be real numbers")
            elif field == "complex":
                for val in values_or_vectors:
                    if not isinstance(val, Complex):
                        raise TypeError("values must be complex numbers")

            self.entries = []
            for i in range(n):
                row = list(values_or_vectors[i * m: (i + 1) * m])
                self.entries.append(row)
            self.n = n
            self.m = m

    def __add__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError("Add with Matrix only")
        if self.field != other.field:
            raise TypeError("Field mismatch")
        if self.n != other.n or self.m != other.m:
            raise ValueError("Dimensions must match")

        result_entries = []
        for i in range(self.n):
            row = []
            for j in range(self.m):
                row.append(self.entries[i][j] + other.entries[i][j])
            result_entries.extend(row)
        return Matrix(self.field, self.n, self.m, *result_entries)

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError("Multiply with Matrix only")
        if self.field != other.field:
            raise TypeError("Fields do not match")
        if self.m != other.n:
            raise ValueError("Dimensions not compatible")

        result_entries = []
        for i in range(self.n):
            row = []
            for j in range(other.m):
                sum_product = self.entries[i][0] * other.entries[0][j]
                for k in range(1, self.m):
                    sum_product += self.entries[i][k] * other.entries[k][j]
                row.append(sum_product)
            result_entries.extend(row)
        return Matrix(self.field, self.n, other.m, *result_entries)

    def get_row(self, i):
        if i < 0 or i >= self.n:
            raise IndexError("Row out of range")
        return Matrix(self.field, 1, self.m, *self.entries[i])

    def get_column(self, j):
        if j < 0 or j >= self.m:
            raise IndexError("Column out of range")
        return Matrix(self.field, self.n, 1, *[self.entries[i][j] for i in range(self.n)])

    def transpose(self):
        transposed_entries = [self.get_column(j).entries for j in range(self.m)]
        transposed_flat = [elem for sublist in transposed_entries for elem in sublist]
        return Matrix(self.field, self.m, self.n, *transposed_flat)

    def conjugate(self):
        conjugated_entries = []
        for row in self.entries:
            conjugated_entries.extend([elem.conjugate() if isinstance(elem, Complex) else elem for elem in row])
        return Matrix(self.field, self.n, self.m, *conjugated_entries)

    def transpose_conjugate(self):
        return self.transpose().conjugate()
