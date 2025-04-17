MAIN.PY
 '''
 BTP project: simulating biological glass using 2D taylor line model for microswimmers.

 There are 3 code files: main.py, propagate.py, plot.py
 -> main.py file has the system initialization code, and after initialization, data is stored into data.txt file in same directory. Also, eventually a plot of the intialized system is shown.
 -> propagate.c is the C code that is required to be compiled first with command 'gcc propagate.c -o <output_filename> -lm' which produces a binary executable file in same directory. The binary code when executed with arguements:'<iterations_count> <input_datafile> <output_datafilename>' reads the data from specified input data file(which stores time elapsed,total iterations done,current radius,freq,phi,positions,velocities,accelerations) and runs the simulation iterations over data for specified number of iterations and then output the data(in same format as input) to specified output_filename. Note: for compiling, run the code in same directory which stores the propagate.c file, else instead of just propagate.c, you need to specify complete location with filename of 'propagate.c'. Also, remember to use different output data filename than input filename, else input data would be lost.
 -> plot.py reads the input data file and then plots the data.

 --- format of datafile
 line1:'<elapsed time>comma<counter till now>comma<current radius>newline'
 line2:'<swimmer1_freq>comma<swimmer2_freq>comma......<lastswimmer_freq>'
 line3:'<swimmer1_phase>comma<swimmer2_phase>comma......<lastswimmer_phase>'
 line4:'<swimmer1_bead1_posx>comma<swimmer1_bead1_posy>comma<swimmer1_bead2_posx>comma<swimmer1_bead2_posy>......<swimmer1_lastbead_posx>comma<posy>\newline'
 line5:'<swimmer1_bead1_velx>comma<swimmer1_bead1_vely>comma<swimmer1_bead2_velx>comma<swimmer1_bead2_vely>......<swimmer1_lastbead_velx>comma<vely>\newline'
 line6:'<swimmer1_bead1_accx>comma<swimmer1_bead1_accy>comma<swimmer1_bead2_accx>comma<swimmer1_bead2_accy>......<swimmer1_lastbead_accx>comma<accy>\newline'
 line7:'<swimmer2_bead1_posx>comma<swimmer2_bead1_posy>comma<swimmer2_bead2_posx>comma<swimmer2_bead2_posy>......<swimmer2_lastbead_posx>comma<posy>\newline'
 line8:'<swimmer2_bead1_velx>comma<swimmer2_bead1_vely>comma<swimmer2_bead2_velx>comma<swimmer2_bead2_vely>......<swimmer2_lastbead_velx>comma<vely>\newline'
 line9:'<swimmer2_bead1_accx>comma<swimmer2_bead1_accy>comma<swimmer2_bead2_accx>comma<swimmer2_bead2_accy>......<swimmer2_lastbead_accx>comma<accy>\newline'
 .... so on
 '''

PROPAGATE.C
This file can do one of two things, depending on the arguments provided
1. Use the file from main.py and run the simulation with a reducing radius, and output a the final configuration.
2. Sample the system with a fixed radius by letting the simulation for some given iterations, and output both the sampled positions and the final state of the system.

CALC.C
This file is used to perform various calculations on the sampled data, behaviour is dependent on arguments provided
1. We can calculate MSD
2. We can calculate the intermediate scattering function for a given wavevector value
3. Given the final state of the system (not the sampled data), make a csv with vx,vy pairs (can plot velocity distributions with it)

PLOT.PY
Used to plot the state of the system

MOVIEMAKER.PY
This code takes the sampled data and creates and saves snapshots of the system.
use appropriate ffmpeg commands to make a movie out of it then
