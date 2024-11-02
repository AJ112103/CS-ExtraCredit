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
