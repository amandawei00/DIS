import numpy as np
import matplotlib.pyplot as plt

"""
    plots for ATLAS 
        hadron: pi0
        y = 1.75
        root(s_NN) = 5.02 TeV
"""
pt1_th = [2.0, 2.5, 3.0, 3.5, 4.0, 4.25, 4.375, 4.53125]
cs1_th = [5.71750794605, 5.1879776782, 4.3578704379, 3.58916300251, 2.94870507051, 2.67696298019, 2.55214870641, 2.40485988424] * np.full(1, (1/10))

pt1_ex = [2.1069518716577535, 2.374331550802139, 2.716577540106952, 3.283422459893048, 3.7647058823529416, 4.256684491978609, 4.7486631016042775, 5.497326203208557, 6.502673796791443]
cs1_ex = [0.5060871986919672, 0.30121267254059286, 0.1792755760934971, 0.07468604063867895, 0.032139732183655324, 0.013389382961641544, 0.00699955487298442, 0.0036591505792521735, 0.0010329623848370942]

"""
    plots for BRAHMS
        hadron: pi0
        y = 4.0
        root(s_NN) = 200 GeV
"""
pt2_th = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
cs2_th = [0.0207549155555, 0.00257231968398, 0.00019720236829, 8.85773575025e-06, 1.50724266754e-07, 8.14671905573e-11]

pt2_ex = [1.0434782608695652, 1.2111801242236024, 1.3788819875776397, 1.5527950310559009, 1.7204968944099377, 1.894409937888199]
cs2_ex = [0.04365158322401665, 0.012022644346174156, 0.004130475019901619, 0.0013182567385564101, 0.0005058246620031147 , 0.00020137242498623895]


################################################ plot
pt_th = pt1_th
cs_th = cs1_th

pt_ex = pt1_ex
cs_ex = cs1_ex

plt.scatter(pt_th, cs_th, marker='o', label="theoretical")
plt.scatter(pt_ex, cs_ex, marker='^', label="experimental")

plt.xlabel("p_T (GeV)")
plt.ylabel("d^2N/dy/d^2p_T (GeV^-2)")

plt.title("pi0 yields in deuteron-gold collisions")

plt.xlim(1.0, 5.0)
plt.ylim(10.e-6, 10.e2)
plt.yscale("log")

plt.legend()
plt.show()






