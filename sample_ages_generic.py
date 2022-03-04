"""Sample age uncertainty distributions in order to get combined estimate
of age of surface"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib

#Rock Creek

data = np.genfromtxt('./data/test_ages.csv', delimiter=',', skip_header=1) [:,2:]
surface_labels = np.genfromtxt('./data/test_ages.csv', delimiter=',', skip_header=1, usecols=0, dtype=str)

# Get samples we want from each surface
# Edit the code below to add additional surfaces.
# Name must match surface labels in first column of data file
age_dict = {'S1':[]}#, 'T2':[], 'T3':[]}
sigma_dict = {'S1':[]}#, 'S2':[], 'S3':[]}
sample_dict = {'S1':[]}#, 'S2':[], 'S3':[]}
for i, label in enumerate(surface_labels):
    if label == 'S1':
        age_dict['S1'].append(data[i][0])
        sigma_dict['S1'].append(data[i][1])
#    elif label == 'S2':
#        age_dict['S2'].append(data[i][0])
#        sigma_dict['S2'].append(data[i][1]) 
#    elif label == 'S3':
#        age_dict['S3'].append(data[i][0])
#        sigma_dict['S3'].append(data[i][1])
    else:
        print('Sample from unknown surface')

# Check we've read it in correctly
print(age_dict)
print(sigma_dict)
for key, value in age_dict.items():
    for i, mu in enumerate(value):
        dnorm = norm(loc=mu, scale=sigma_dict[key][i])
        samples = dnorm.rvs(size=10000)
        sample_dict[key].append(samples)
    # Make an example plot
    if key=='S1':
        # Estimate bounds from data, may need to adjust manually (see line below)
        xvals = np.arange(min(value)-6*max(sigma_dict[key]), max(value)+6*max(sigma_dict[key]), 0.1)
#        xvals = np.arange(0, 40, 0.1)
        for i, mu in enumerate(value):
            dnorm = norm(loc=mu, scale=sigma_dict[key][i])
            plt.plot(xvals, dnorm.pdf(xvals), c='0.6')
            p1, = plt.fill(xvals, dnorm.pdf(xvals), facecolor='0.8', label='Individual sample') 
        mu_samples = np.mean(sample_dict[key])
        sig_samples = np.std(sample_dict[key])
        dnorm_samp = norm(loc=mu_samples, scale=sig_samples)
        plt.plot(xvals, dnorm_samp.pdf(xvals), c='0.3', zorder=11)
        p2, = plt.fill(xvals, dnorm_samp.pdf(xvals), facecolor='0.4',
                       zorder=10, label = 'Combined age')#, alpha=0.5)
        plt.xlabel('Age (ka)')
        plt.ylabel('Density')
        plt.legend(handles=[p1,p2])
        plt.savefig('plots/sampling_method_test_example.png', dpi=300)
        
print('\nTest data')
for key, value in sample_dict.items():
    mean = np.mean(value)
    sigma = np.std(value)
    print(key)
    print('Mean age (ka)', mean)
    print('Sigma (kyr)', sigma)
