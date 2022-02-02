"""Sample age uncertainty distributions in order to get combined estimate
of age of surface"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib

#Rock Creek

data = np.genfromtxt('./data/age_offset_rock_creek.csv', delimiter=',', skip_header=1) [:,1:]
# Convert to ka
data = data/1000
labels = np.genfromtxt('./data/age_offset_rock_creek.csv', delimiter=',', skip_header=1, usecols=0, dtype=str)

# Get samples we want from each surface
age_dict = {'T2':[], 'T4':[], 'T6':[]}
sigma_dict = {'T2':[], 'T4':[],  'T6':[]}
sample_dict = {'T2':[], 'T4':[],  'T6':[]}    
for i, label in enumerate(labels):
    if label[2] == '2':
        if label[3] == 'A':
            pass
        else:
            age_dict['T2'].append(data[i][0])
            sigma_dict['T2'].append(data[i][1])  
    elif label == 'BT301' or label == 'BS401' or label == 'BS402':
        age_dict['T4'].append(data[i][0])
        sigma_dict['T4'].append(data[i][1]) 
    elif label[2] == '6':
        if label[4]=='4' or label[4]=='5':
            pass # Ignore exhumed bedrock samples
        else:
            age_dict['T6'].append(data[i][0])
            sigma_dict['T6'].append(data[i][1])

for key, value in age_dict.items():
    for i, mu in enumerate(value):
        dnorm = norm(loc=mu, scale=sigma_dict[key][i])
        samples = dnorm.rvs(size=10000)
        sample_dict[key].append(samples)
print('Rock Creek')
for key, value in sample_dict.items():
    mean = np.mean(value)
    sigma = np.std(value)
    print(key)
    print('mean', mean)
    print('sigma', sigma)

# Neds Creek
data = np.genfromtxt('./data/age_offset_neds_creek.csv', delimiter=',', skip_header=1) [:,1:]
# Convert to ka
data = data/1000
labels = np.genfromtxt('./data/age_offset_neds_creek.csv', delimiter=',', skip_header=1, usecols=0, dtype=str)

# Get samples we want from each surface
age_dict = {'T1':[], 'T3':[], 'T4':[], 'T5':[]}
sigma_dict = {'T1':[], 'T3':[], 'T4':[], 'T5':[]}
sample_dict = {'T1':[], 'T3':[], 'T4':[], 'T5':[]}

for i, label in enumerate(labels):
    if label[3] == '1':
        age_dict['T1'].append(data[i][0])
        sigma_dict['T1'].append(data[i][1])
    elif label[3] == '3':
        age_dict['T3'].append(data[i][0])
        sigma_dict['T3'].append(data[i][1]) 
    elif label == 'NCT403' or label == 'NCT406' or label == 'NCFW404':
        age_dict['T4'].append(data[i][0])
        sigma_dict['T4'].append(data[i][1]) 
    elif label == 'NCT502' or label == 'NCT503':
        age_dict['T5'].append(data[i][0])
        sigma_dict['T5'].append(data[i][1])
        
for key, value in age_dict.items():
    for i, mu in enumerate(value):
        dnorm = norm(loc=mu, scale=sigma_dict[key][i])
        samples = dnorm.rvs(size=10000)
        sample_dict[key].append(samples)
    # Make an example plot
    if key=='T3':
        xvals = np.arange(60, 120, 0.1)
        for i, mu in enumerate(value):
            dnorm = norm(loc=mu, scale=sigma_dict[key][i])
            plt.plot(xvals, dnorm.pdf(xvals), c='0.6')
            plt.fill(xvals, dnorm.pdf(xvals), facecolor='0.8') 
        mu_samples = np.mean(sample_dict[key])
        print('musamples', mu_samples)
        sig_samples = np.std(sample_dict[key])
        print(sig_samples)
        dnorm_samp = norm(loc=mu_samples, scale=sig_samples)
        plt.plot(xvals, dnorm_samp.pdf(xvals), c='0.3', zorder=11)
        plt.fill(xvals, dnorm_samp.pdf(xvals), facecolor='0.4', zorder=10)#, alpha=0.5)
        plt.xlabel('Age (ka)')
        plt.ylabel('Density')
        plt.savefig('plots/sampling_method.png', dpi=300)
        
print('\nNeds Creek')
for key, value in sample_dict.items():
    mean = np.mean(value)
    sigma = np.std(value)
    print(key)
    print('mean', mean)
    print('sigma', sigma)
