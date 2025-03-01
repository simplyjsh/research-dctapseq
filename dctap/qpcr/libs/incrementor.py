class Incrementor:
    def __init__(self, value=0):
        self.value = value

    def __int__(self):
        return self.value

    def __index__(self):
        return self.value

    def __str__(self):
        return str(self.value)

    def __add__(self, other):
        if isinstance(other, Incrementor):
            return Incrementor(self.value + other.value)
        elif isinstance(other, int):
            return Incrementor(self.value + other)
        return NotImplemented

    def __radd__(self, other):
        # This covers the case int + Incrementor
        if isinstance(other, int):
            return Incrementor(other + self.value)
        return NotImplemented

    def step(self):
        self.value += 1

    def skip(self, steps: int):
        self.value += steps
