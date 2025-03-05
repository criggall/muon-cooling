import subprocess
from matplotlib import pyplot as plt

file = '/Users/criggall/Documents/muon-cooling/HFOFO-G4bl/sol_place7_31.txt'

fig_dir = '/Users/criggall/Documents/muon-cooling/Figures/'

period_length = 4200 #mm

with open(file, 'r') as f:
        lines = f.readlines()
        
        current_vals = []
        z_vals = []
        for i in range(len(lines)):
            if i > 10 and i < 210:
                start_ind = lines[i].find("current=") + len("current=")
                stop_ind = lines[i].find(" rotation=") - 1
                current = lines[i][start_ind:stop_ind]
                if current != '':
                    # current_vals.append(float(current))
                    current_vals.append(abs(float(current)))

                start_ind_z = lines[i].find("z=") + len("z=")
                stop_ind_z = lines[i].find("+$period*") - 1
                start_ind_pnum = lines[i].find("+$period*") + len("+$period*") 
                stop_ind_pnum = lines[i].find(" current") 
                z = lines[i][start_ind_z:stop_ind_z]
                pnum = lines[i][start_ind_pnum:stop_ind_pnum]
                if z != '' and pnum != '':
                    z_vals.append((float(z)+float(pnum)*period_length)/1000)

plt.scatter(z_vals,current_vals)
plt.ylabel('Magnitude of solenoid current')
plt.xlabel('z (m)')
plt.savefig(fig_dir+'current_vals_vs_z.png',dpi=300)