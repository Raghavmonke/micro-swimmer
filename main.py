# '''
# BTP project: simulating biological glass using 2D taylor line model for microswimmers.
#
# There are 3 code files: main.py, propagate.py, plot.py
# -> main.py file has the system initialization code, and after initialization, data is stored into data.txt file in same directory. Also, eventually a plot of the intialized system is shown.
# -> propagate.c is the C code that is required to be compiled first with command 'gcc propagate.c -o <output_filename> -lm' which produces a binary executable file in same directory. The binary code when executed with arguements:'<iterations_count> <input_datafile> <output_datafilename>' reads the data from specified input data file(which stores time elapsed,total iterations done,current radius,freq,phi,positions,velocities,accelerations) and runs the simulation iterations over data for specified number of iterations and then output the data(in same format as input) to specified output_filename. Note: for compiling, run the code in same directory which stores the propagate.c file, else instead of just propagate.c, you need to specify complete location with filename of 'propagate.c'. Also, remember to use different output data filename than input filename, else input data would be lost.
# -> plot.py reads the input data file and then plots the data.
#
# --- format of datafile
# line1:'<elapsed time>comma<counter till now>comma<current radius>newline'
# line2:'<swimmer1_freq>comma<swimmer2_freq>comma......<lastswimmer_freq>'
# line3:'<swimmer1_phase>comma<swimmer2_phase>comma......<lastswimmer_phase>'
# line4:'<swimmer1_bead1_posx>comma<swimmer1_bead1_posy>comma<swimmer1_bead2_posx>comma<swimmer1_bead2_posy>......<swimmer1_lastbead_posx>comma<posy>\newline'
# line5:'<swimmer1_bead1_velx>comma<swimmer1_bead1_vely>comma<swimmer1_bead2_velx>comma<swimmer1_bead2_vely>......<swimmer1_lastbead_velx>comma<vely>\newline'
# line6:'<swimmer1_bead1_accx>comma<swimmer1_bead1_accy>comma<swimmer1_bead2_accx>comma<swimmer1_bead2_accy>......<swimmer1_lastbead_accx>comma<accy>\newline'
# line7:'<swimmer2_bead1_posx>comma<swimmer2_bead1_posy>comma<swimmer2_bead2_posx>comma<swimmer2_bead2_posy>......<swimmer2_lastbead_posx>comma<posy>\newline'
# line8:'<swimmer2_bead1_velx>comma<swimmer2_bead1_vely>comma<swimmer2_bead2_velx>comma<swimmer2_bead2_vely>......<swimmer2_lastbead_velx>comma<vely>\newline'
# line9:'<swimmer2_bead1_accx>comma<swimmer2_bead1_accy>comma<swimmer2_bead2_accx>comma<swimmer2_bead2_accy>......<swimmer2_lastbead_accx>comma<accy>\newline'
# .... so on
# '''


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
np.random.seed(44)


# Simulation parameters
N_beads = 10  # Number of beads in one system
N_systems = 5 # Number of bead systems
K = 1000000.0  # Spring constant (Hooke's Law)
m = 1.0  # Mass of each bead, maybe 1 microgram
r = 0.1  # Radius of each bead
collisionR = 2.0*r   # center to center spacing for collision
L0 = 2.0 * r  # Natural length of the spring
K_bending = 100000.0     # Bending rigidity constant
dt = 0.0005  # Time step
R = 2  #confinement initial radius
dRdt = -0.5    # confinement shrink rate
d2Rdt = 0.01   # rate of change of shrink rate
freq = np.random.uniform(20, 30 , N_systems)  # Random frequencies for each system, intial velocities of COM of each microwimmer depends on this
phi = np.random.uniform(0, 2 * np.pi, N_systems)  # Random angles for each system
b = 0.7   # amplitude parameter
L0b = L0*b
min_distance = 3.0 * r # minimum distance between two system intially
KbT = 2
var = KbT/m     # variance for normal sampling dist
t = 0  # Initial global time
counter = 0  # global counter
nu = 100    # coupling strength
nuxdt = nu*dt   # andersen collision probability
vf  = (N_systems*N_beads*np.pi*r**2)*100 / (np.pi*(R+r)**2) # intial occupied volume percent

# Initialize positions, velocities and accelerations for all systems
positions = np.zeros((N_systems, N_beads, 2))
velocities = np.zeros((N_systems, N_beads, 2))
ac = np.zeros((N_systems, N_beads, 2))   # accelerations

for i in range(N_systems):
    a = np.random.random(1)
    if (a<0.5):
        velocities[i,:,0] = -1*freq[i]*2*L0
    else:
        velocities[i,:,0] = freq[i]*2*L0




def generate_random_position_within_circle(R, system_length):
    """Generate a random center for the COM microswimmer within the circle of radius R,
    ensuring that no bead extends beyond the boundary."""
    while True:
        # Random angle and radius for the center
        theta = np.random.uniform(0, 2 * np.pi)
        radius = np.random.uniform(0, R - system_length / 2)  # Ensure no bead is outside the boundary
        
        # Convert polar coordinates to Cartesian for the center
        center_x = radius * np.cos(theta)
        center_y = radius * np.sin(theta)
        
        # Ensure the system length fits within the circular region
        if np.sqrt(center_x ** 2 + center_y ** 2) + system_length / 2 <= R:
            return np.array([center_x, center_y])

def check_minimum_distance(new_bead_positions, existing_positions, min_distance):
    """Check if the new bead positions are at least `min_distance` away from existing beads."""
    for existing_system in existing_positions:
        for existing_bead in existing_system:
            for new_bead in new_bead_positions:
                distance = np.linalg.norm(existing_bead - new_bead)
                if distance < min_distance:
                    return False
    return True

# Randomize the positions for each bead-spring system inside the circle
for i in range(N_systems):
    system_length = (N_beads - 1) * L0  # Total length of the bead-spring system
    
    while True:
        center = generate_random_position_within_circle(R, system_length)
        
        # Generate x-positions for the beads within the system, centered around the system's center
        x_positions = np.linspace(center[0] - system_length / 2, center[0] + system_length / 2, N_beads)
        new_bead_positions = np.column_stack((x_positions, np.full(N_beads, center[1])))
        
        # Check if the new bead positions are at least `min_distance` away from existing beads
        if check_minimum_distance(new_bead_positions, positions[:i], min_distance):
            break  # The new positions are valid, break the loop

    # Assign x and y positions for all beads in the system
    positions[i, :, 0] = x_positions
    positions[i, :, 1] = center[1]  # Keep the y-position the same for all beads in the system


# If you want to see animation, uncomment all commented parts. Note: only small system N_systems<10 could be animated easily.
# Initialize plot
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-R - 5, R + 5)
ax.set_ylim(-R - 5, R + 5)
circle = Circle((0, 0), R, color='black', fill=False)
ax.add_patch(circle)
bead_systems = [ax.plot([], [], 'b-')[0] for _ in range(N_systems)];bead_systems.append(circle)

# def read_from_backup(f_name):
#     '''this method reads the simulation values from filename given as parameter'''
#     global t,R,positions,velocities,ac,freq,phi,counter
#     with open(f_name,"r") as f1:
#         s = f1.readline().strip().split(",") ; t = float(s[0]) ; counter = int(s[1]) ;R = float(s[2])
#         v = f1.readline().strip().split(",")
#         p = f1.readline().strip().split(",")
#         for i in range(N_systems):
#             freq[i] = float(v[i])
#             phi[i] = float(p[i])
#         for i in range(N_systems):
#             pos = f1.readline().strip().split(",")
#             vel = f1.readline().strip().split(",")
#             acc = f1.readline().strip().split(",")
#             for j in range(N_beads):
#                 positions[i][j][0] = float(pos[2*j]) ; positions[i][j][1] = float(pos[2*j+1])
#                 velocities[i][j][0] = float(vel[2*j]) ; velocities[i][j][1] = float(vel[2*j+1])
#                 ac[i][j][0] = float(acc[2*j]) ; ac[i][j][1] = float(acc[2*j+1])


def write_to_backup(f_name,mode="w"):
    '''this method writes the simulation values(elapsed time,iterations,current radius, positions,velocities,accelerations) into filename given as first parameter
     In second parameter, if "a" is given then it will append the new date at the end of specified file. else by default it will erase earlier data and then write new data.'''
    with open(f_name,mode) as f1:
        f1.write(f"{t},{counter},{R}\n")
        for i in range(N_systems-1):
            f1.write(f"{freq[i]},")
        f1.write(f"{freq[N_systems-1]}\n")
        for i in range(N_systems-1):
            f1.write(f"{phi[i]},")
        f1.write(f"{phi[N_systems-1]}\n")
        for i in range(N_systems):
            for j in range(N_beads-1):
                f1.write(f"{positions[i][j][0]},{positions[i][j][1]},")
            f1.write(f"{positions[i][N_beads-1][0]},{positions[i][N_beads-1][1]}\n")
            for j in range(N_beads-1):
                f1.write(f"{velocities[i][j][0]},{velocities[i][j][1]},")
            f1.write(f"{velocities[i][N_beads-1][0]},{velocities[i][N_beads-1][1]}\n")
            for j in range(N_beads-1):
                f1.write(f"{ac[i][j][0]},{ac[i][j][1]},")
            f1.write(f"{ac[i][N_beads-1][0]},{ac[i][N_beads-1][1]}\n")

#write_to_backup("data.txt","w") # writing intial system values to data.txt file

def propel():
    '''this method calculates the hooke's and bending forces, then it detects and handles collision dynamics.'''
    global ac,velocities
    col_vel = np.zeros((N_systems, N_beads, 2)) #aggregates velocity changes during collisions
    for j in range(N_systems):
        # handling collison dynamics
        for k in range(j+1,N_systems):
            for a in range(N_beads):
                for b in range(N_beads):
                    rab = positions[j][a]-positions[k][b]
                    r_mag = np.linalg.norm(rab)
                    if (r_mag<collisionR):
                        bab = np.dot(rab,(velocities[j][a]-velocities[k][b]))
                        if (bab<0):
                            # moving half time step back
                            positions[j][a][0] -= (velocities[j][a][0] + (ac[j][a][0]*dt/4)) * dt/2
                            positions[j][a][1] -= (velocities[j][a][1] + (ac[j][a][1]*dt/4)) * dt/2
                            positions[k][b][0] -= (velocities[k][b][0] + (ac[k][b][0]*dt/4)) * dt/2
                            positions[k][b][1] -= (velocities[k][b][1] + (ac[k][b][1]*dt/4)) * dt/2
                            #projecting center-to-center spacing to collisonR
                            rab_ = positions[j][a]-positions[k][b]
                            r_mag_ = np.linalg.norm(rab_)
                            dab = ( (collisionR/r_mag_ - 1) /2)*rab_ 
                            positions[k][b] -= dab
                            positions[j][a] += dab
                            # calculating change in velocities due to collision
                            vab = ( bab / (r_mag*r_mag) ) * rab
                            col_vel[j][a] -= vab
                            col_vel[k][b] += vab
    # adding changes to velocities due to collision and moving half time step further
    for j in range(N_systems):
        for i in range(N_beads):
            if (col_vel[j][i][0]!=0.0 and col_vel[j][i][1]!=0.0):
                velocities[j][i] += col_vel[j][i]
                positions[j][i][0] += (velocities[j][i][0] + (ac[j][i][0]*dt/4)) * dt/2
                positions[j][i][1] += (velocities[j][i][1] + (ac[j][i][1]*dt/4)) * dt/2
    
    # calculating forces
    for j in range(N_systems):
        
        # Calculation of hookean and bending force starts here
        F1 = np.zeros((N_beads, 2))  # Hookean force
        F2 = np.zeros((N_beads, 2))  # Bending force
        ta = np.zeros((N_beads, 2))  # Vector distance between adjacent beads
        alpha = np.zeros(N_beads)  # Rotation angle
        cos_alpha = np.zeros(N_beads)
        sin_alpha = np.zeros(N_beads)


        # Calculate alpha using the provided formula
        for i in range(N_beads):
            alpha[i] = L0b * np.sin(phi[j] + 2 * np.pi * (freq[j] * t + i * (2.0 / N_beads)))
            cos_alpha[i] = np.cos(alpha[i])
            sin_alpha[i] = np.sin(alpha[i])

        # Calculating tangents and Hookean forces
        for i in range(N_beads-1):
            ta[i][0] = positions[j][i+1][0] - positions[j][i][0]
            ta[i][1] = positions[j][i+1][1] - positions[j][i][1]
            length = (ta[i][0]**2+ta[i][1]**2)**0.5
            Kdl_n = K*(length - L0) / length
            f1 = Kdl_n*ta[i][0] ; f2 = Kdl_n * ta[i][1]
            F1[i][0] += f1 ; F1[i][1] += f2
            F1[i+1][0] -= f1 ; F1[i+1][1] -= f2


        #bending force calculations
        F2[0][0] -= ta[0][0] ; F2[0][1] -= ta[0][1]
        F2[1][0] += ta[0][0] ; F2[1][1] += ta[0][1]
        for i in range(0,N_beads-2):

            rtix = cos_alpha[i]*ta[i][0] - sin_alpha[i]*ta[i][1]
            rtiy = cos_alpha[i]*ta[i][1] + sin_alpha[i]*ta[i][0]
            F2[i+2][0] += rtix ; F2[i+2][1] += rtiy
            F2[i+1][0] -= rtix ; F2[i+1][1] -= rtiy

            rttiux = 2*ta[i][0] - cos_alpha[i]*ta[i+1][0] - sin_alpha[i]*ta[i+1][1]
            rttiuy = 2*ta[i][1] - cos_alpha[i]*ta[i+1][1] + sin_alpha[i]*ta[i+1][0]
            F2[i][0] += rttiux ; F2[i][1] += rttiuy
            F2[i+1][0] -= rttiux ; F2[i+1][1] -= rttiuy

        F2[N_beads-2][0] += ta[N_beads-2][0] ; F2[N_beads-2][1] += ta[N_beads-2][1]
        F2[N_beads-1][0] -= ta[N_beads-2][0] ; F2[N_beads-1][1] -= ta[N_beads-2][1]
        F2 *= K_bending

        # aggregating Hooke's and bending force'
        for i in range(N_beads):
            ac[j][i][0] = (F1[i][0]+F2[i][0]) / m
            ac[j][i][1] = (F1[i][1]+F2[i][1]) / m
            




def update_positions():
    '''this function uses propel function to calculate forces and handle collisions and then using new forces calculates velocity and then positions.
It also updates the value of radius based on given shrink rate. at the end it incremements the time and counter.'''

    global positions,velocities,ac,t,R,dRdt,counter

    #calculate new radius of system
    R += dRdt*dt
    dRdt += d2Rdt*dt
    circle.radius = R
    # calculate new accelerations
    propel()
    # Update velocities and then positions
    for j in range(N_systems):
        for i in range(N_beads):
            velocities[j][i][0] += ac[j][i][0] * dt
            velocities[j][i][1] += ac[j][i][1] * dt
            positions[j][i][0] += (velocities[j][i][0] + (ac[j][i][0]*dt/2)) * dt
            positions[j][i][1] += (velocities[j][i][1] + (ac[j][i][1]*dt/2)) * dt
            bounce_forward_boundary(positions[j][i],velocities[j][i])

        if(counter % 100 == 0):   # andersen thermostat step
            if (np.random.rand() < nuxdt):
                print(counter,j)
                velocities[j] = np.random.normal(loc=0, scale=var, size=(N_beads,2))  # Mean of 0, std dev = KbT/m
    t += dt
    counter+=1




def bounce_forward_boundary(bead_position, bead_velocity):
    '''This function handles the dynamics of collisions with the bounding walls, using bounce forward rule,adapted for moving bounding wall'''
    # Calculate the distance of the bead from the center
    distance_from_center = np.linalg.norm(bead_position)

    # Check if the bead is outside the boundary
    if distance_from_center > R:
        bead_position -= bead_velocity*dt/2
        r = np.linalg.norm(bead_position)
        # project the bead back to the boundary
        bead_position *= R / r  # Scale position to stay on the boundary

        # Reflect the velocity along the normal direction to "bounce back" but keep the tangential component
        # Using the formula for velocity reflection: v = v - 2 * (v ⋅ n) * n
        bead_velocity -=  ( (2 * (np.dot(bead_velocity, bead_position)/R) - dRdt )/R) * bead_position
        bead_position += bead_velocity*dt/2

def animate(i):
    update_positions()
    for j in range(N_systems):
        bead_systems[j].set_data(positions[j, :, 0], positions[j, :, 1])
    return bead_systems

# Set up animation
ani = animation.FuncAnimation(fig, animate, frames=200, blit=True,cache_frame_data=True, repeat=True, interval=dt*1000)
plt.show()
plt.close()

# dRdt = -0.0001
# d2Rdt = 0
# ani = animation.FuncAnimation(fig, animate, frames=200, blit=True,cache_frame_data=True, repeat=True, interval=dt*1000)

vf  = (N_systems*N_beads*np.pi*r**2)*100 / (np.pi*(R+r)**2)
print(vf,R,type(vf))
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
        circle = plt.Circle((positions[j][i][0],positions[j][i][1]), 0.1, color='r', fill=True)
        ax.add_artist(circle)

plt.title(f"Bead-Spring Systems {vf}% volume fraction at time {t}")
plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.show()
plt.close()


