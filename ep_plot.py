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


""" #qsq2 = 800.0 GeV^2
x = [0.001, 0.00133352, 0.00177828, 0.00237137, 0.00316228]
y1 = [1.3759847033681756, 1.361175146291811, 1.34940189281754, 1.341227964841948, 1.3359151697300244]
y2 = [2.54532613667, 2.3602336871219314, 2.1854479088331846, 2.0204814977947048, 1.8648601005260643]
# y = [ 4.461676029943761, 3.7281722840899203, 3.0931746386507255, 2.545326136672553, 2.0744102498269568, 1.6710809420193744, 1.3271801894209667]
y2 = [y2[i]/2 for i in range(len(y2))] """

# run 20, qsq2 = 10. GeV^2
"""x = np.logspace(-4, -2.5, 5)
y = [4.461676029943761, 3.64365533535095, 2.948456918411565, 2.3602336871219314, 1.8648601005260643]
y = [(y[i]/max(y)) * 1.4 for i in range(len(y))]

# plt.xlim(1.e-6, 1.e-2)
"""

# run 21, qsq2 = 8000 GeV^2
x = [0.18, 0.25, 0.4, 0.65]
y_calc = [4.904070283, 4.61380915, 4.347081647, 4.189467577]
# y_calc = [y_calc[i] - 4.0 for i in range(len(y_calc))]
# y_calc = [(y_calc[i]/max(y_calc))* 0.2781 for i in range(len(y_calc))]
y_exp = [0.2781, 0.2022, 0.1006, 0.0103]
# plt.xlim(0.1, 1.)


plt.plot(x, y_calc, color='r')
plt.plot(x, y_exp)
plt.xscale('log')
plt.xlim(0.1, 1.)
plt.ylim(0, 5)
plt.show()
