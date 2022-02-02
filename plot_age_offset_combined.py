"""Summary plot showing both Neds Creek and Rock Creek data, and 
paleoearthquake record
"""

import os
import numpy as np
import matplotlib.pyplot as plt                                                                             
import matplotlib
from adjustText import adjust_text


xvals = np.arange(-10000, 60000, 1000)
xvals_short = np.arange(-10000, 40000, 1000)
# Rock creek short term
yvals_rock_short = -0.1983433732041342 + 0.00014382463319831778*xvals_short
plt.plot(xvals_short, yvals_rock_short, c='k', linestyle = 'dashed')

#Rock creek long term
yvals_rock_average = -1.4440933992406588 + 0.00019889438886251442*xvals
plt.plot(xvals, yvals_rock_average, c='k', label='Hyde')

#Hyde eq ages
h_eq_x = [10.3, 10.3, 22.8, 22.8, 38.8, 38.8, 47.2, 47.2]
h_eq_y = [0, 2., 2., 4., 4.0, 6.0, 6.0, 8.0]
h_eq_x = np.array(h_eq_x)*1000
plt.plot(h_eq_x, h_eq_y, c='k', linestyle='dotted')

# Neds Creek short term
xvals_short = np.arange(-10000, 27000, 1000) 
yvals_neds_since25 = 0.0 + 0.00021229974162166333*xvals_short
plt.plot(xvals_short, yvals_neds_since25, c='b', linestyle = 'dashed')

xvals_pre25 = np.arange(26000, 40000, 1000)
yvals_neds_pre25 = 4.648695516517866 - 1.04 + 7.256672601451256e-05*xvals_pre25
plt.plot(xvals_pre25, yvals_neds_pre25, c='b', linestyle = 'dashed')

yvals_neds_average = 2.5295892563814872 + 8.202296863648996e-05*xvals
plt.plot(xvals, yvals_neds_average, c='b', label='Dunstan')

# Neds paleoearthquake ages
d_eq_x = [14.1, 14.1, 15.3, 15.3, 16.6, 16.6, 19.0, 19.0, 26.0] 
d_eq_x = np.array(d_eq_x)*1000
d_eq_y = [0, 1.375, 1.375, 2.75, 2.75, 4.125, 4.125, 5.5, 5.5]
plt.plot(d_eq_x, d_eq_y, c='b', linestyle='dotted')


plt.xlim(000, 35000)
plt.ylim(-4, 8)
plt.plot([0,0],[-4, 12], c='0.5', linestyle = 'dotted')
plt.plot([-10000,50000],[0, 0], c='0.5', linestyle = 'dotted') 

ax = plt.gca()
ax.legend(loc=2) 
ax.set_xlabel('Age (years)')
ax.set_ylabel('Vertical displacement (m)') 
figname = './plots/combined_age_offset.png'
plt.savefig(figname)
