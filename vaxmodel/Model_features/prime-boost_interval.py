import numpy as np
import matplotlib.pyplot as plt
import json
from model import Model

params = json.load(open('vaxmodel/parameters.json', 'r'))

model = Model()
params['M'] = 2
times = [0, 21, 200] # Timings of doses
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.figure(figsize=(8, 6))
plt.semilogy(t, antibodies, lw=2, label='3 weeks')
times = [0, 28, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='4 weeks')
times = [0, 35, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='5 weeks')
times = [0, 42, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='6 weeks')
times = [0, 56, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='8 weeks')
times = [0, 70, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='10 weeks')
times = [0, 84, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='12 weeks')
times = [0, 126, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='18 weeks')
times = [0, 168, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
plt.semilogy(t, antibodies, lw=2, label='24 weeks')
plt.axhline(y=0.48, color='black', linestyle='--', label='LOD')
plt.ylabel(r'IgG Antibody ($\mu$g/ml)', labelpad=10, fontsize=20)
plt.xlabel('Days post initial vaccine', labelpad=10, fontsize=20)
plt.legend(fontsize=14, ncol=2)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim(0, times[params['M']])
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# Naive

plt.figure(figsize=(8, 6))
model = Model()
params['M'] = 2
times = [0, 28, 200] # Timings of doses
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
four_week = max(antibodies)
times = [0, 42, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
six_week = max(antibodies)
times = [0, 56, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
eight_week = max(antibodies)
times = [0, 70, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
ten_week = max(antibodies)
times = [0, 84, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'])
twelve_week = max(antibodies)
plt.bar([4, 6, 8, 10, 12], [four_week, six_week, eight_week, ten_week, twelve_week], width=1.5)
plt.ylabel(r'IgG Antibody ($\mu$g/ml)', labelpad=10, fontsize=20)
plt.xlabel('Dosing interval (weeks)', labelpad=10, fontsize=20)
plt.yscale('log')
plt.ylim(10**1, 10**3)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# Data

plt.figure(figsize=(8, 6))
y = [10**5.239061486504216, 10**5.5372379737877155, 10**5.414109850782877, 10**5.414109850782877, 10**5.585337110074525]
plt.bar([4, 6, 8, 10, 12], y, width=1.5)
# add error bars
error = [[10**4.92276201274602, 10**5.367485268872713, 10**5.26631965854948, 10**5.113710403370401, 10**5.376777952801864],
        [10**5.451700950579383, 10**5.633854058801847, 10**5.643607904867345, 10**5.742079977713819, 10**5.8726668238040105]]
plt.errorbar([4, 6, 8, 10, 12], y, yerr=[[y[i] - error[0][i] for i in range(len(y))], [error[1][i] - y[i] for i in range(len(y))]],
            fmt='none', ecolor='black', capsize=5)
plt.ylabel(r'AU/ml ($log_{10}$)', labelpad=10, fontsize=20)
plt.xlabel('Dosing interval (weeks +/-7 days)', labelpad=10, fontsize=20)
plt.yscale('log')
plt.ylim(10**4, 10**6)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# Recovered

plt.figure(figsize=(8, 6))
model = Model()
params['M'] = 2
times = [0, 28, 200] # Timings of doses
t, antibodies, antigen = model.simulate(times, params, doses=params['M'], naive=False)
four_week = max(antibodies)
times = [0, 56, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'], naive=False)
eight_week = max(antibodies)
times = [0, 70, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'], naive=False)
ten_week = max(antibodies)
times = [0, 84, 200]
t, antibodies, antigen = model.simulate(times, params, doses=params['M'], naive=False)
twelve_week = max(antibodies)
plt.bar([4, 6, 8, 10], [four_week, eight_week, ten_week, twelve_week], width=1.5, color='red')
plt.ylabel(r'IgG Antibody ($\mu$g/ml)', labelpad=10, fontsize=20)
plt.xlabel('Dosing interval (weeks)', labelpad=10, fontsize=20)
plt.yscale('log')
plt.ylim(10**1, 10**3)
plt.xticks([4, 6, 8, 10], ['4', '8', '10', '12'], fontsize=18)
plt.yticks(fontsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# Data

plt.figure(figsize=(8, 6))
y = [10**5.6415094339622645, 10**5.528301886792453, 10**5.528301886792453, 10**5.6415094339622645]
plt.bar([4, 6, 8, 10], y, width=1.5, color='red')
# add error bars
error = [[10**5.283018867924528, 10**5.452830188679245, 10**5.39622641509434, 10**5.433962264150944],
        [10**6, 10**5.69811320754717, 10**5.735849056603773, 10**5.867924528301886]]
plt.errorbar([4, 6, 8, 10], y, yerr=[[y[i] - error[0][i] for i in range(len(y))], [error[1][i] - y[i] for i in range(len(y))]],
            fmt='none', ecolor='black', capsize=5)
plt.ylabel(r'AU/ml ($log_{10}$)', labelpad=10, fontsize=20)
plt.xlabel('Dosing interval (weeks +/-7 days)', labelpad=10, fontsize=20)
plt.yscale('log')
plt.ylim(10**4, 10**7)
plt.xticks([4, 6, 8, 10], ['4', '8', '10', '12'], fontsize=18)
plt.yticks(fontsize=18)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.show()