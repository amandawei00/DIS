import numpy as np
import matplotlib.pyplot as plt

"""
# qsq2 = 10 GeV^2
x = np.logspace(-4, -2, 5)
y = [4.461676, 3.39909, 2.656326, 1.8648, 1.327180]
y = [y[i] / 4 for i in range(len(y))]


#qsq2 = 0.5 GeV^2
x = np.logspace(-5, -3, 5)
y = [0.6341769607426025, 0.5345270131919134, 0.45033344746914367, 0.3794213100572658, 0.32015807115119566]
y = [y[i] / 4 for i in range(len(y))]


#qsq2 = 1.5 GeV^2
x = np.logspace(-5, -2, 9)
y = [1.7994275397059067, 1.5704603213640764, 1.3684580622326927, 1.1906973459438115, 1.0347522528646285, 0.8983555849367808, 0.7795312073752181, 0.676448278203514, 0.5875345516031919]
y = [y[i] / 4 for i in range(len(y))]




#qsq2 = 2.5 GeV^2
x = np.logspace(-4, -2, 9)
y = []
y = [y[i]/4 for i in range(len(y))]



# qsq2 = 5.0 GeV^2
x = np.logspace(-4, -2, 7)
y = []
y = [y[i]/4 for i in range(len(y))] """


#qsq2 = 800.0 GeV^2
x = [0.0185, 0.021, 0.0242, 0.0322]
y = [0.002141647, 0.0022550336137, 0.002387356, 0.002672106]
y = [y[i]/4 for i in range(len(y))]

plt.plot(x, y)
plt.xscale('log')
plt.xlim(10.e-6, 10.e-2)
plt.ylim(0, 1.8)
plt.show()