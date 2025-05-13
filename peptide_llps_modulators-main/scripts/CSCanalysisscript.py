'''
nephelostar data analysis
'''

import pandas as pd
import peptides as pep
import matplotlib.pyplot as plt
import re
import numpy as np
import matplotlib as mpl
import seaborn as sns
mpl.style.use('default')

def timereadout(time_value):
    minutes = '0'
    seconds = '0'
    # parse hours
    match = re.search(r'\d+\s*min', time_value)
    if match:
        minutes = re.search(r'\d+', match.group()).group()
    # parse minutes
    match = re.search(r'\d+\s*s', time_value)
    if match:
        seconds = re.search(r'\d+', match.group()).group()
    # returns hours * 60 + minutes 
    return int(minutes) * 60 + int(seconds)

# curdir = r'\\huckdfs-srv.science.ru.nl\huckdfs\RobotLab\AsCatPep\AsCatPep_PSM\AsCatPep_PSM0009\salttitrations'
curdir = './data'
outdir = './out'
df = pd.read_csv(rf'{curdir}/AsCatPep_PSM0009_1.txt', header=0, sep='\t')

#df = df.rename(columns={0: 'sequence'})
time = df.loc[0][2:].values
df = df.loc[1:]
df = df.set_index('Well')

dflayout = pd.read_csv(rf'{curdir}/AsCatPep_PSM0009_1.csv', header=0, sep=',')
dflayout = dflayout.dropna()
dflayout = dflayout.rename(columns={'titration plate position':'Well'})
dflayout = dflayout.set_index('Well')

# todo
# load sequecnes file from below here, merge/add column of 'sequences' to this dflayout
# def map_OT2wellID_to_reactorPos(OT2wellid):
#     ''' maps OT2-formulation well plate #ID to 
#         biotage synth well ID'''
#     # return wellID
# def map_reactorPos_to_sequence(reactorBlockId):
#     ''' maps biotage synth well ID to sequence '''
#     # return wellID
# def map_id_to_seq(reactorPos):
#     wellid = map_OT2wellID_to_reactorPos(reactorPos)
#     seq = map_reactorPos_to_sequence(wellid)
#     return seq
# df['sequence'] = df['reactor_block_position'].apply(map_id_to_seq)


df = dflayout.join(df)
df = df.drop(columns = 'Content')

print(df)
print(df.shape)
print('\n'*5)

samplelist = df.index
nosamples = len(df.index)

#going through the samples and calculating integral (timestep*signal, normalised by area of rectangle 2e6 in height and time(end)-time(start) broad

times = [timereadout(x) for x in time]
timestep = times[1] - times[0]
timespan = times[-1] - times[0]
A_max = 2000000*timespan

# this one keeps track of outputs
integrals = []

print(df)
print(dflayout)

for i in samplelist:
    turb = df.loc[i][1:].values
    print(turb)
    turb = [int(x) for x in turb]
    normturb = [x/max(turb) for x in turb]
    
    integral = sum(turb*timestep)
    integrals.append(integral)
    
    plt.plot(times, turb)
    
print('num integrals', len(integrals), integrals)
plt.savefig(f'{outdir}/turbidity_traces.png')
plt.clf()

#normalisation by the reference (no peptide) and with background (buffer) subtracted
df['integral'] = integrals

A_ref = np.average(df.loc[df['reactor block position']=='reference']['integral'])
A_blank = np.average(df.loc[df['reactor block position']=='blank']['integral'])

normintegrals = [(x-A_blank)/(A_ref-A_blank) for x in integrals]

df['normintegral'] = normintegrals

sns.barplot(data=df, x='integral', y='reactor block position')
plt.show()
sns.barplot(data=df, x='normintegral', y='reactor block position')
plt.show()

# ? missing: join reactor_block_pos / well_id --> sequence
# read_csv(sequences_file)
# join this df with df from above

df.to_csv(f'./out/csc_integrals_sequences_df.csv')
