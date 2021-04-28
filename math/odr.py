#!/usr/bin/env python3
"""Orthogonal (geometrical) linear regression example script.

    See also
    --------
    The example is in analogy to
    http://docs.scipy.org/doc/scipy/reference/odr.html

    The various fields of `output` are explained here:
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.Output.html

    Author
    ------
    Andreas Anhaeuser (2016)
"""

import numpy as np
import scipy.odr as odr
import matplotlib.pyplot as plt

from misc.math.linreg import lin_reg

#======== SETUP ========#
m_true = 1.3
t_true = -0.1
xnoise = 0.2
ynoise = 0.2
N = 10**3
#=======================#

# create ideal data set
x = np.random.random(N)
y = m_true * x + t_true

# add noise:
x = x + np.random.normal(0, xnoise, N)
y = y + np.random.normal(0, ynoise, N)

# vertical least square regression:
# alternativeley, use:     a, b = 2, 0
a, b, Da, Db, N = lin_reg(x, y)

# orthogonal least square regression:
def f(beta, x):
    return beta[0] * x + beta[1]           # fitting function

# odr needs a reasonable first guess. Use the vertical regression result. An
# alternative is to use beta0 = [1., 0.]
beta0 = [a, b]                           # first guess
model = odr.Model(f)
mydata = odr.Data(x, y)
myodr = odr.ODR(mydata, model, beta0)
output = myodr.run()
m, t = output.beta                       # best estimate
Dm, Dt = output.sd_beta                  # uncertainties

# plot:
text_title = 'Orthogonal vs. vertical linear regression'
text_truth = 'Truth: m = %1.2f, t= %1.2f' % (m_true, t_true)
text_vert = (
        'Vertical: m = %1.2f (+-%1.2f), t = %1.2f (+-%1.2f)'
        % (a, Da, b, Db)
        )
text_orth = (
        'Orthogonal: m = %1.2f (+-%1.2f), t = %1.2f (+-%1.2f)'
        % (m, Dm, t, Dt)
        )
plt.plot(x, y, 'k.', label=text_truth)
plt.plot(x, a * x + b, 'b-', label=text_vert)
plt.plot(x, m * x + t, 'r-', label=text_orth)
plt.legend(loc='upper left', fontsize=9)
plt.title(text_title)
plt.axis('square')
plt.show()
