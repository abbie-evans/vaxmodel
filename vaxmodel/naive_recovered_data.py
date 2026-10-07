import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

goel_data = pd.read_csv('Goel long data.csv')
# remove rows with Y values in recovered
goel_data = goel_data[~goel_data['recovered'].isin(['N'])]
goel_data = goel_data[~goel_data['vaccine_type'].isin(['Moderna'])]
goel_data = goel_data[~goel_data['timepoint_factor'].isin(['t7', 't8', 't9', 't10'])]
goel_data = goel_data[goel_data['days_dose_1'] >= 0]
# remove rows where 'serum_spike_IgG' is NaN
goel_data = goel_data[~goel_data['serum_spike_IgG'].isna()]
grouped = goel_data.groupby('days_dose_1')['serum_spike_IgG']
# for each days_dose_1, compute the geometric mean by multiplying the values and taking the nth root
geometric_means = {}
for name, group in grouped:
    geometric_means[name] = group.prod()**(1/len(group))  # Geometric mean by multiplying values and taking the nth root
means = pd.Series(geometric_means)
goel_antibody_data = means.to_numpy()
goel_antibody_days = means.index.to_numpy()
# at day 21, split the means into two parts
goel_antibody_days_split = [means[means.index <= 21].index.to_numpy(),
                             means[(means.index > 21) & (means.index <= 200)].index.to_numpy()]
goel_antibody_data_split = [means[means.index <= 21].to_numpy(),
                             means[(means.index > 21) & (means.index <= 200)].to_numpy()]
goel_antibody_data_r = goel_antibody_data[goel_antibody_days <= 200]
goel_antibody_days_r = goel_antibody_days[goel_antibody_days <= 200]

goel_data = pd.read_csv('Goel long data.csv')
# remove rows with Y values in recovered
goel_data = goel_data[~goel_data['recovered'].isin(['Y'])]
goel_data = goel_data[~goel_data['vaccine_type'].isin(['Moderna'])]
goel_data = goel_data[~goel_data['timepoint_factor'].isin(['t7', 't8', 't9', 't10'])]
goel_data = goel_data[goel_data['days_dose_1'] >= 0]
# remove rows where 'serum_spike_IgG' is NaN
goel_data = goel_data[~goel_data['serum_spike_IgG'].isna()]
grouped = goel_data.groupby('days_dose_1')['serum_spike_IgG']
# for each days_dose_1, compute the geometric mean by multiplying the values and taking the nth root
geometric_means = {}
for name, group in grouped:
    geometric_means[name] = group.prod()**(1/len(group))  # Geometric mean by multiplying values and taking the nth root
means = pd.Series(geometric_means)
goel_antibody_data = means.to_numpy()
goel_antibody_days = means.index.to_numpy()
# at day 21, split the means into two parts
goel_antibody_days_split = [means[means.index <= 21].index.to_numpy(),
                             means[(means.index > 21) & (means.index <= 200)].index.to_numpy()]
goel_antibody_data_split = [means[means.index <= 21].to_numpy(),
                             means[(means.index > 21) & (means.index <= 200)].to_numpy()]
goel_antibody_data = goel_antibody_data[goel_antibody_days <= 200]
goel_antibody_days = goel_antibody_days[goel_antibody_days <= 200]

goel_mod_data = pd.read_csv('Goel moderna data.csv')
# remove rows with Y values in recovered
goel_mod_data = goel_mod_data[~goel_mod_data['recovered'].isin(['N'])]
goel_mod_data = goel_mod_data[~goel_mod_data['vaccine_type'].isin(['Pfizer'])]
goel_mod_data = goel_mod_data[~goel_mod_data['timepoint_factor'].isin(['t7', 't8', 't9', 't10'])]
goel_mod_data = goel_mod_data[goel_mod_data['days_dose_1'] >= 0]
# remove rows where 'serum_spike_IgG' is NaN
goel_mod_data = goel_mod_data[~goel_mod_data['serum_spike_IgG'].isna()]
grouped = goel_mod_data.groupby('days_dose_1')['serum_spike_IgG']
# for each days_dose_1, compute the geometric mean by multiplying the values and taking the nth root
geometric_means = {}
for name, group in grouped:
    geometric_means[name] = group.prod()**(1/len(group))  # Geometric mean by multiplying values and taking the nth root
mod_means = pd.Series(geometric_means)
goel_mod_antibody_data = mod_means.to_numpy()
goel_mod_antibody_days = mod_means.index.to_numpy()
# at day 21, split the means into two parts
goel_mod_antibody_days_split = [mod_means[mod_means.index <= 28].index.to_numpy(),
                             mod_means[(mod_means.index > 28) & (mod_means.index <= 200)].index.to_numpy()]
goel_mod_antibody_data_split = [mod_means[mod_means.index <= 28].to_numpy(),
                             mod_means[(mod_means.index > 28) & (mod_means.index <= 200)].to_numpy()]
goel_mod_antibody_data_r = goel_mod_antibody_data[goel_mod_antibody_days <= 200]
goel_mod_antibody_days_r = goel_mod_antibody_days[goel_mod_antibody_days <= 200]

goel_mod_data = pd.read_csv('Goel moderna data.csv')
# remove rows with Y values in recovered
goel_mod_data = goel_mod_data[~goel_mod_data['recovered'].isin(['Y'])]
goel_mod_data = goel_mod_data[~goel_mod_data['vaccine_type'].isin(['Pfizer'])]
goel_mod_data = goel_mod_data[~goel_mod_data['timepoint_factor'].isin(['t7', 't8', 't9', 't10'])]
goel_mod_data = goel_mod_data[goel_mod_data['days_dose_1'] >= 0]
# remove rows where 'serum_spike_IgG' is NaN
goel_mod_data = goel_mod_data[~goel_mod_data['serum_spike_IgG'].isna()]
grouped = goel_mod_data.groupby('days_dose_1')['serum_spike_IgG']
# for each days_dose_1, compute the geometric mean by multiplying the values and taking the nth root
geometric_means = {}
for name, group in grouped:
    geometric_means[name] = group.prod()**(1/len(group))  # Geometric mean by multiplying values and taking the nth root
mod_means = pd.Series(geometric_means)
goel_mod_antibody_data = mod_means.to_numpy()
goel_mod_antibody_days = mod_means.index.to_numpy()
# at day 21, split the means into two parts
goel_mod_antibody_days_split = [mod_means[mod_means.index <= 28].index.to_numpy(),
                             mod_means[(mod_means.index > 28) & (mod_means.index <= 200)].index.to_numpy()]
goel_mod_antibody_data_split = [mod_means[mod_means.index <= 28].to_numpy(),
                             mod_means[(mod_means.index > 28) & (mod_means.index <= 200)].to_numpy()]
goel_mod_antibody_data = goel_mod_antibody_data[goel_mod_antibody_days <= 200]
goel_mod_antibody_days = goel_mod_antibody_days[goel_mod_antibody_days <= 200]

goel_antibody_days_weeks = goel_antibody_days[np.isin(goel_antibody_days, [0, 7, 14, 21, 28, 42, 60, 91, 120, 150, 180, 200])]
index = [np.where(goel_antibody_days == day)[0][0] for day in goel_antibody_days_weeks]
goel_antibody_data_weeks = goel_antibody_data[index]

goel_antibody_days_r_weeks = goel_antibody_days_r[np.isin(goel_antibody_days_r, [0, 7, 14, 21, 28, 42, 60, 91, 112, 150, 180, 200])]
index_r = [np.where(goel_antibody_days_r == day)[0][0] for day in goel_antibody_days_r_weeks]
goel_antibody_data_r_weeks = goel_antibody_data_r[index_r]

goel_mod_antibody_days_weeks = goel_mod_antibody_days[np.isin(goel_mod_antibody_days, [0, 7, 14, 21, 28, 35, 42, 60, 90, 96, 120, 150, 181, 200])]
index_mod = [np.where(goel_mod_antibody_days == day)[0][0] for day in goel_mod_antibody_days_weeks]
goel_mod_antibody_data_weeks = goel_mod_antibody_data[index_mod]

goel_mod_antibody_days_r_weeks = goel_mod_antibody_days_r[np.isin(goel_mod_antibody_days_r, [0, 7, 14, 21, 28, 42, 60, 90, 120, 150, 181, 200])]
index_mod_r = [np.where(goel_mod_antibody_days_r == day)[0][0] for day in goel_mod_antibody_days_r_weeks]
goel_mod_antibody_data_r_weeks = goel_mod_antibody_data_r[index_mod_r]

plt.figure(figsize=(5, 6))
plt.plot(goel_antibody_days_weeks, goel_antibody_data_weeks, '-', label='Naïve', color='blue')
plt.plot(goel_antibody_days_r_weeks, goel_antibody_data_r_weeks, '-', label='Recovered', color='red')
for time in [0, 21]:
    plt.axvline(x=time, ymin=0, ymax=1, color='gray', linestyle='--', alpha=0.5)
plt.axhline(0.48, color='black', linestyle='--', label='LOD', lw=2)
plt.yscale('log')
plt.xlim(0, 100)
plt.xlabel('Days post initial vaccine', fontsize=20, labelpad=10)
plt.ylabel(r'IgG Antibody ($\mu$g/mL)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=16)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

plt.figure(figsize=(5, 6))
plt.plot(goel_mod_antibody_days_weeks, goel_mod_antibody_data_weeks, '--', label='Naïve', color='blue')
plt.plot(goel_mod_antibody_days_r_weeks, goel_mod_antibody_data_r_weeks, '--', label='Recovered', color='red')
for time in [0, 28]:
    plt.axvline(x=time, ymin=0, ymax=1, color='gray', linestyle='--', alpha=0.5)
plt.axhline(0.48, color='black', linestyle='--', label='LOD', lw=2)
plt.yscale('log')
plt.xlim(0, 100)
plt.xlabel('Days post initial vaccine', fontsize=20, labelpad=10)
plt.ylabel(r'IgG Antibody ($\mu$g/mL)', fontsize=20, labelpad=10)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(fontsize=16)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()
