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
          return ((self.real**2 + self.img**2)**0.5)
    
    def conjugate(self):
          return Complex(self.real, -(self.img))
    
class Vector:
    def __init__(self, field, n, *coordinates):
        if field not in ["real", "complex"]:
            raise ValueError("Field must be either 'real' or 'complex'")
        
        self.field = field
        self.n = n

        if len(coordinates) != n:
            raise ValueError(f"Expected {n} coordinates, but got {len(coordinates)}")

        if field == "real":
            for coord in coordinates:
                if not isinstance(coord, (int, float)):
                    raise TypeError("All coordinates must be real numbers (int or float)")
        elif field == "complex":
            for coord in coordinates:
                if not isinstance(coord, Complex):
                    raise TypeError("All coordinates must be complex numbers (instances of Complex)")

        self.coordinates = list(coordinates)

