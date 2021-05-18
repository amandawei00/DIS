import numpy as np
import matplotlib.pyplot as plt

def x(l, n):

    x = []
    for i in range(n):
        a = float(l/n) * (i + 1 - 0.5)
        x.append(a)

    print(x)
    return x

a = [1.0881, 6.4935]
b = [1546.875, 4374.375]

plt.plot(a,b)
plt.show()