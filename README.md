# OtagoCosmo_publish
Clean code for for publishing calculations of slip rates from cosmo data in Otago. Supports:
Griffin, Stirling, Wilcken and Barrell: Late Quaternary slip rates for the Hyde and Dunstan faults,
southern New Zealand: implications for strain migration in a slowly deforming continental plate margin

Contains the following code:

`sample_ages.py`
This combines multiple CRN ages to estimate surface ages by taking Monte Carlo samples from
the uncertainty distribution of each individual age, combining these, and then taking the
mean and standard deviation assuming a normal distribution.

`age_offset_samples_rockcreek.py and age_offset_samples_nedscreek.py`
Take Monte Carlo samples for each age-offset combination, and then fit
long-term and piecwise linear fits to the data, to estimate uplift rates
and their uncertainties.