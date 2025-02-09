import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
np.random.seed(44)


# Simulation parameters
N_beads = 10  # Number of beads in one system
N_systems = 5 # Number of bead systems
m = 1.0  # Mass of each bead, maybe 1 microgram
r = 0.1  # Radius of each bead
R = 2  # intializing R
t = 0.0 #initializing time
counter = 0 #initializing counter
# Initialize positions, velocities and accelerations for all systems
positions = np.zeros((N_systems, N_beads, 2))
velocities = np.zeros((N_systems, N_beads, 2))
ac = np.zeros((N_systems, N_beads, 2))   # accelerations


def read_from_backup(f_name):
    global t,R,positions,velocities,ac,freq,phi,counter
    with open(f_name,"r") as f1:
        s = f1.readline().strip().split(",") ; t = float(s[0]) ; counter = int(s[1]) ;R = float(s[2])
        f1.readline();f1.readline()
        for i in range(N_systems):
            pos = f1.readline().strip().split(",")
            vel = f1.readline().strip().split(",")
            acc = f1.readline().strip().split(",")
            for j in range(N_beads):
                positions[i][j][0] = float(pos[2*j]) ; positions[i][j][1] = float(pos[2*j+1])
                velocities[i][j][0] = float(vel[2*j]) ; velocities[i][j][1] = float(vel[2*j+1])
                ac[i][j][0] = float(acc[2*j]) ; ac[i][j][1] = float(acc[2*j+1])


# reading data
read_from_backup("out.txt")

# calculating occupied volume fraction
vf = (N_systems*N_beads*np.pi*r**2)*100 / (np.pi*(R+r)**2)

# plot the Final Configuration
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-R - 5, R + 5)
ax.set_ylim(-R - 5, R + 5)
circle = Circle((0, 0), R, color='black', fill=False)
ax.add_patch(circle)
# plot each bead-spring system
for j in range(N_systems):
    x = positions[j, :, 0]
    y = positions[j, :, 1]
    ax.plot(x, y, 'b-')  # plot springs
    for i in range(N_beads):
        circle = plt.Circle((positions[j][i][0],positions[j][i][1]), r, color='r', fill=True)
        ax.add_artist(circle)

plt.title(f"Bead-Spring Systems {vf}% volume fraction at time {t}s")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()



