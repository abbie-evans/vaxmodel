import numpy as np
import matplotlib.pyplot as plt
import json
from model import Model

params = json.load(open('vaxmodel/parameters.json', 'r'))

model = Model()
times = [0, 21, 100] # Timings of doses
params['M'] = 2
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.figure(figsize=(8, 6))
plt.semilogy(t, antibodies, label='Naïve', color='blue', lw=2)
t, antibodies, antigen = model.simulate(times, params, doses=params['M'], naive=False)
plt.semilogy(t, antibodies, label='Recovered', color='red', lw=2, linestyle='--')
plt.xlim(0, times[params['M']])
for time in times[:params['M']]:
    plt.axvline(x=time, ymin=0, ymax=1, color='gray', linestyle='--', alpha=0.5, lw=2)
plt.xlabel('Days post initial vaccine', labelpad=10, fontsize=20)
plt.ylabel(r'IgG Antibody ($\mu$g/ml)', labelpad=10, fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.axhline(y=0.48, color='black', linestyle='--', label='LOD', lw=2)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.legend(fontsize=16)
plt.tight_layout()
plt.show()

model = Model()
times = [0, 21, 100] # Timings of doses
params['M'] = 2
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
# split into the two doses
naive_doses_response = []
start_idx = 0
for i in range(1, len(times)):
    end_time = times[i]
    end_idx = np.searchsorted(t, end_time, side='right')
    naive_doses_response.append((t[start_idx:end_idx], antibodies[start_idx:end_idx]))
    start_idx = end_idx
t, antibodies, antigen = model.simulate(times, params, doses=params['M'], naive=False)
recovered_doses_response = []
start_idx = 0
for i in range(1, len(times)):
    end_time = times[i]
    end_idx = np.searchsorted(t, end_time, side='right')
    recovered_doses_response.append((t[start_idx:end_idx], antibodies[start_idx:end_idx]))
    start_idx = end_idx
first_fold_increase = max(naive_doses_response[0][1])/naive_doses_response[0][1][0] # ratio of max antibody after first dose to initial antibody before first dose for naive
second_fold_increase = max(naive_doses_response[1][1])/naive_doses_response[1][1][0] # ratio of max antibody after second dose to max antibody after first dose for naive
first_fold_increase_rec = max(recovered_doses_response[0][1])/recovered_doses_response[0][1][0] # ratio of max antibody after first dose to max antibody before first dose for recovered
second_fold_increase_rec = max(recovered_doses_response[1][1])/recovered_doses_response[1][1][0] # ratio of max antibody after second dose to max antibody after first dose for recovered

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
fig.subplots_adjust(hspace=0.05)
ax1.bar([1, 3], [first_fold_increase, second_fold_increase], width=1.5, color='blue', label='Naïve')
ax1.bar([6, 8], [first_fold_increase_rec, second_fold_increase_rec], width=1.5, color='red', label='Recovered')
ax2.bar([1, 3], [first_fold_increase, second_fold_increase], width=1.5, color='blue')
ax2.bar([6, 8], [first_fold_increase_rec, second_fold_increase_rec], width=1.5, color='red')
ax1.set_ylim(35, 50)
ax2.set_ylim(0, 10)
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)
ax2.xaxis.tick_bottom()
ax1.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

d = .5
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0], [0], transform=ax1.transAxes, **kwargs)
ax2.plot([0], [1], transform=ax2.transAxes, **kwargs)

fig.supylabel('Fold change in antibody titre', fontsize=18, x=0.02)
plt.xticks([2, 7], ['Prime', 'Boost'], fontsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax2.tick_params(axis='y', labelsize=18)
ax2.yaxis.get_major_locator().set_params(integer=True)
fig.legend(fontsize=16)
plt.tight_layout()
plt.show()