import pandas as pd

column_names = ['x', 'y', 'z', 'Px', 'Py', 'Pz', 't', 'PDGid', 'EventID', 'TrackID', 'ParentID', 'Weight']
df = pd.read_csv('initial.dat', sep='\s+', skiprows=3, names=column_names)

t_shift = 0.5 * 1/(325e6)

print(df)

for i in range(len(df['t'].values)):
    df['t'] = df['t']-t_shift
    df['PDGid'] *= -1

df = df[df.PDGid == 13]

df.to_csv('initial_muminus.txt', sep='\t', index=False, header=True) 