""""Plot fault offset against age of surface. Sample uncertainty space to produce uplift rate curves, following
method of Gold and Cowgill (2011) Deriving fault-slip histories to test for secular variation in slip,
 with examples from the Kunlun and Awatere faults, EPSL 301 (1-2).
"""

import os
import numpy as np
import matplotlib.pyplot as plt                                                                             
import matplotlib
from adjustText import adjust_text
from scipy.stats import truncnorm
from scipy import optimize
from pylr2 import regress2

data = np.genfromtxt('./data/age_offset_rock_creek_sampled.csv', delimiter=',', skip_header=1) [:,1:]
labels = np.genfromtxt('./data/age_offset_rock_creek_sampled.csv', delimiter=',', skip_header=1, usecols=0, dtype=str)
data_all = np.genfromtxt('./data/age_offset_rock_creek.csv', delimiter=',', skip_header=1) [:,1:]
labels_all = np.genfromtxt('./data/age_offset_rock_creek.csv', delimiter=',', skip_header=1, usecols=0, dtype=str) 
#print(data)
#print(labels)

def connect(ends):
    """From https://stackoverflow.com/questions/47704008/fastest-way-to-get-all-the-points-between-two-x-y-coordinates-in-python"""
    d0, d1 = np.abs(np.diff(ends, axis=0))[0]
    if d0 > d1: 
        return np.c_[np.linspace(ends[0, 0], ends[1, 0], d0+1),
                     np.linspace(ends[0, 1], ends[1, 1], d0+1)]
    else:
        return np.c_[np.linspace(ends[0, 0], ends[1, 0], d1+1),
                     np.linspace(ends[0, 1], ends[1, 1], d1+1)]

# Make list of colours depending on whether constructional or erosional surface
colours = []
i=0
for j,lab in enumerate(labels):
    if lab[0]=='B':
        colours.append('k')
    else:
        colours.append('k')
# Normal run
plt.scatter(data[:,0], data[:,2], c=colours, zorder=10)
plt.errorbar(data[:,0], data[:,2], xerr=data[:,1]*2, yerr=data[:,3], fmt='none', c=colours, zorder=9)

# Now plot raw data
plt.scatter(data_all[:,0], data_all[:,2], c='0.8', marker = 's', s=12)
plt.errorbar(data_all[:,0], data_all[:,2], xerr=data_all[:,1]*2, yerr=data_all[:,3],
             fmt='none', c='0.8')

data_age = data[:,0]
data_age_sigma = data[:,1]
slip_data_age = np.concatenate([np.array([0]),data_age])

data_slip = data[:,2]
data_slip_sigma = data[:,3] 
slip_data_offset = np.concatenate([np.array([0]),data_slip])
#print(data_slip)
#print(data_slip_sigma)

# Now we do the incremental slip rate calculations
# Make lists of truncated normal distributions for each observation
# of offset and ages. Distribution is truncated at 2 sigma
age_dist = []
offset_dist = []
for i, age in enumerate(data_age):
    # Get bounds to truncate distribution with
    a, b = -2*data_age_sigma[i] / data_age_sigma[i], \
        2*data_age_sigma[i] / data_age_sigma[i]
    age_tn = truncnorm(a, b, loc = age, scale=data_age_sigma[i])
    age_dist.append(age_tn)
    
# Now do equivalent for offset measurements
for i, offset in enumerate(data_slip):  
    a, b = -2*data_slip_sigma[i]/2 / (data_slip_sigma[i]/2), \
        2*data_slip_sigma[i]/2 / (data_slip_sigma[i]/2)    
    offset_tn = truncnorm(a, b, loc = data_slip[i], scale=(data_slip_sigma[i]/2)) 
    offset_dist.append(offset_tn)

# Now sample from distributions
n = 1000 # Number of samples
age_samples = [np.zeros(n)] # Start with origin point 
offset_samples = [np.zeros(n)] # Start with origin point 
for i, age_tn in enumerate(age_dist):
    sample_age = age_tn.rvs(size=n)
    age_samples.append(sample_age)
    sample_offset = offset_dist[i].rvs(size=n)
    offset_samples.append(sample_offset)

age_samples = np.array(age_samples)
offset_samples = np.array(offset_samples)      
# Need to check curve is monotonic
mono_age = np.diff(age_samples.T) > 0

# Get indices of monotonically increasing samples for age 
mono_ind1 = np.where(np.all(mono_age, axis=1))[0]
# Now for offset
mono_offset = np.diff(offset_samples.T) > 0 
mono_ind2 = np.where(np.all(mono_offset, axis=1))[0]
ind = np.intersect1d(mono_ind1, mono_ind2)
age_samples = age_samples.T[ind]
offset_samples = offset_samples.T[ind]

# Now we want to do a least squares fit to each slip rate curve sample and find median value
ls_slopes = []
ls_intercepts = []
ls_slopes_a = []
ls_intercepts_a = []
ls_slopes_b = []
ls_intercepts_b = []
for i, x in enumerate(age_samples):
    y = offset_samples[i]
    results = regress2(x, y, _method_type_2="reduced major axis")
    ls_slopes.append(results['slope'])
    ls_intercepts.append(results['intercept'])
    results = regress2(x[:4], y[:4], _method_type_2="reduced major axis")
    ls_slopes_a.append(results['slope'])
    ls_intercepts_a.append(results['intercept'])
    results = regress2(x[3:], y[3:], _method_type_2="reduced major axis")
    ls_slopes_b.append(results['slope'])
    ls_intercepts_b.append(results['intercept']) 
bval = np.median(ls_slopes)
bval_median = bval
aval = np.median(ls_intercepts)
print('bval', bval)
print('aval', aval)
xvals = np.arange(0, 125000, 500)
yvals = bval*xvals + aval
plt.plot(xvals, yvals, c='k')
plt.xlim([0, 1.1*max(xvals)])

# Print some of the fits
bval_a = np.median(ls_slopes_a)
bval_a_median = bval_a
aval_a = np.median(ls_intercepts_a)
break_point = np.mean(age_samples.T[3])
print('break_point', break_point)
xvals_a = np.arange(0, break_point, 500)  
yvals_a = bval_a*xvals_a + aval_a
print('bval_a', bval_a)
print('aval_a', aval_a)
xvals_b = np.arange(break_point, 125000, 500)
bval_b = np.median(ls_slopes_b)
bval_b_median = bval_b
aval_b = np.median(ls_intercepts_b)
yvals_b = bval_b*xvals_b + aval_b
print('bval_b', bval_b)
print('aval_b', aval_b) 
plt.plot(xvals_a, yvals_a, c='r')
plt.plot(xvals_b, yvals_b, c='r')


# Now do confidence intervals
perc = [2.5, 16, 84, 97.5]
linestyles = ['dashed', 'dashdot','dashdot','dashed']
for i, p in enumerate(perc):
    bval = np.percentile(ls_slopes, p)
    aval = np.percentile(ls_intercepts, p)
    yvals = bval*xvals + aval
    plt.plot(xvals, yvals, c='0.2', linestyle = linestyles[i], linewidth = '0.8')
    # Now for piecwise bits
    bval_a = np.percentile(ls_slopes_a, p)
    bval_b = np.percentile(ls_slopes_b, p)
    # Get limits on b value
    if i == 0:
        bval_ll = bval
        bval_ll_a = bval_a
        bval_ll_b = bval_b
    if i == len(perc)-1:
        bval_ul = bval
        bval_ul_a = bval_a
        bval_ul_b = bval_b
ax = plt.gca()
print(bval_ll, bval_ul)

# Add some labels
txt = r'${:6.2f}_{{:4.2f}}^{{:4.2f}}$ mm/yr'.format(bval_median*1000, bval_ll*1000, bval_ul*1000)
txt = r'${%.3f}_{{-%.3f}}^{+%.3f}$ mm/yr' % (bval_median*1000, (bval_median - bval_ll)*1000,
                                             (bval_ul-bval_median)*1000) 
ax.annotate(txt, (80000, 14.5), xytext = (40000, 31),
            arrowprops=dict(arrowstyle="->"), fontsize=10)
txt = r'${%.3f}_{{-%.3f}}^{+%.3f}$ mm/yr' % (bval_a_median*1000, (bval_a_median - bval_ll_a)*1000,
                                             (bval_ul_a-bval_a_median)*1000) 
ax.annotate(txt, (30000, 4.2), xytext = (5000, 16),
            arrowprops=dict(arrowstyle="->", color='r'), fontsize=10, color='r')
txt = r'${%.3f}_{{-%.3f}}^{+%.3f}$ mm/yr' % (bval_b_median*1000, (bval_b_median - bval_ll_b)*1000,
                                             (bval_ul_b-bval_b_median)*1000)
ax.annotate(txt, (108000, 20), xytext = (80000, 33.0),
            arrowprops=dict(arrowstyle="->", color='r'), fontsize=10, color='r')  

# Here we do piecwise fits along each sample.
line_xvals = []
line_yvals = []
for i, x in enumerate(age_samples):
    y = offset_samples[i]
    for j in range(len(x)):
        if j==0:
            pass
        elif j == 1:
            ends = np.array([[x[j-1], y[j-1]],
                             [x[j], y[j]]])
            line_segment = connect(ends)
            line = line_segment
        else:
            ends = np.array([[x[j-1], y[j-1]],
                              [x[j], y[j]]])
            line_segment = connect(ends)
            line = np.append(line, line_segment, axis=0)
    if i ==0:
        line_xvals = line.T[0]
        line_yvals = line.T[1]
    else:
        line_xvals = np.hstack([line_xvals, line.T[0]])
        line_yvals = np.hstack([line_yvals, line.T[1]])

# Convert xvals to integers
line_xvals = np.around(line_xvals).astype(int)

# Get every 100th point
x100 = np.arange(0, max(data_all[:,0]), 100)# max(line_xvals)-10000, 100)
perc = [2.5, 16, 50, 84, 97.5]
yperc = []
linestyles = ['dashed', 'dashdot','solid', 'dashdot','dashed']
colours = ['wheat', 'wheat', 'orange', 'orange', 'wheat']
colours = ['darkorange', 'darkorange', 'maroon', 'maroon',  'darkorange']
colours = ['salmon', 'salmon', 'cornflowerblue', 'cornflowerblue', 'salmon']
for x in x100:
    ind = np.where(line_xvals == x)[0]
    ylist = []
    for p in perc:
        yp = np.percentile(line_yvals[ind], p)
        ylist.append(yp)
    yperc.append(ylist)
yperc = np.array(yperc).T

# Now plot each of them
for i, p in enumerate(perc):
#    plt.plot(x100, yperc[i], c='0.9', linewidth=0.3)
    if i == 0:
        pass
    else:
        ax.fill_between(x100, yperc[i-1], yperc[i], color=colours[i], alpha=0.7, lw=0)

# Now plot median line
plt.plot(x100, yperc[2], c='0.5', linewidth=2)  
# Add dummy data for legend
p1 = plt.scatter([10000],[350], c='r')
p2 = plt.scatter([10000],[350], c='k')
ax = plt.gca()
ax.set_xlabel('Age (years)')
ax.set_ylabel('Vertical displacement (m)')
ax.set_ylim([0,37])

# Annotate points
texts = []
for i,lab in enumerate(labels):
    text = ax.annotate(lab, (data[i,0], data[i,2]))
    texts.append(text)
adjust_text(texts)
if not os.path.exists('./plots'):
    os.mkdir('plots')
plt.savefig('plots/age_offset_rockcreek_sampled.png')
