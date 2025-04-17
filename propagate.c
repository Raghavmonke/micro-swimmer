#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// constants of simulation
#define N_beads 20  //Number of beads used in taylor line
#define N_systems 120 //Number of systems
#define m 1.0    //mass of each bead
#define K 1000000.0 //hookes' constant
#define r 0.1   // radius of beads
#define L0 0.2      //natural length of bonds
#define dt 0.0005 // time step         //
#define dt_half 0.00025 //half timestep
#define K_bending 100000.0 // bending rigidity constant
#define collisionR 0.2
#define collisionR_sq 0.04
#define L0b 0.14 // L0 times b
#define TWO_PI 6.28318530718
#define Var 2 // KbT/m
#define d2Rdt 0.002  // rate of change of shrink rate
#define nu 100.0    // andersen coupling strength

int iter;   //number of iterations
double R;   // radius of confinement
double dRdt = -0.8; //intial shrink rate
double positions[N_systems][N_beads][2];
double velocities[N_systems][N_beads][2];
double ac[N_systems][N_beads][2];
double freq[N_systems];   //frequencies of swimmers
double phi[N_systems];  // initial phase of swimmers
double t;   // global time
int counter;    // global counter
double nudt = nu * dt;
FILE *input;
FILE *output;




// Generates two random numbers from normal distribution wiht variance KbT/M and stores in provided array as arguement.
void rand_normal(double* ar) {
    double u1 = (double)rand() / RAND_MAX;
    if (u1 < 1e-10) u1 = 1e-10;
    double u2 = ((double)rand() / RAND_MAX) * TWO_PI;
    double z0 = sqrt(-2.0 * Var * log(u1));
    ar[0] = z0 * cos(u2);
    ar[1] = z0 * sin(u2);
}
// Calculates the bending and hooke's forces on each beads, then detects collisions and handles collision dynamics
void propel() {
    double col_vel[N_systems][N_beads][2] = {0.0};      //aggregates change in velocity due to collisions
    for (int j = 0; j < N_systems; ++j) {
        double F1[N_beads][2] = {0.0}; // Hookean force
        double F2[N_beads][2] = {0.0}; // Bending force
        double ta[N_beads][2] = {0.0}; // Vector distance between adjacent beads
        double alpha[N_beads] = {0.0}; // Rotation angle
        double cos_alpha[N_beads] = {0.0};
        double sin_alpha[N_beads] = {0.0};

        // Handling Collision dynamics
        for (int k = j + 1; k < N_systems; ++k) {
            for (int x = 0; x < N_beads; ++x) {
                for (int y = 0; y < N_beads; ++y) {
                  double diff[2] = {positions[j][x][0] - positions[k][y][0],positions[j][x][1] - positions[k][y][1]};
                  double mag_sq = diff[0]*diff[0] + diff[1]*diff[1];
                  if (mag_sq < collisionR_sq) {
                    double bxy = ( diff[0]*(velocities[j][x][0]-velocities[k][y][0]) + diff[1]*(velocities[j][x][1]-velocities[k][y][1]) ) / mag_sq ;
                    if (bxy<0){
                        //moving half time step back
                        // positions[j][x][0] -= (velocities[j][x][0] + (ac[j][x][0]*dt/4)) * dt/2;
                        // positions[j][x][1] -= (velocities[j][x][1] + (ac[j][x][1]*dt/4)) * dt/2;
                        // positions[k][y][0] -= (velocities[k][y][0] + (ac[k][y][0]*dt/4)) * dt/2;
                        // positions[k][y][1] -= (velocities[k][y][1] + (ac[k][y][1]*dt/4)) * dt/2;

                        //projecting center-to-center spacing to collisionR, so overlapping beads become just-touching beads
                        double diff_new[2] = {positions[j][x][0] - positions[k][y][0],positions[j][x][1] - positions[k][y][1]};
                        double mag_sq_new  = diff_new[0]*diff_new[0] + diff_new[1]*diff_new[1];
                        double dxy[2] = {((collisionR/sqrt(mag_sq_new) -1)/2)*diff_new[0],((collisionR/sqrt(mag_sq_new) -1)/2)*diff_new[1]};
                        positions[k][y][0]-=dxy[0];
                        positions[k][y][1]-=dxy[1];
                        positions[j][x][0]+=dxy[0];
                        positions[j][x][1]+=dxy[1];
                        // printf("bounce by %f  , %f",dxy[0],dxy[1]);
                        //calculating change in velocities, following elastic collision rules
                        double vxy[2] = {bxy * diff[0], bxy * diff[1]};
                        col_vel[j][x][0] -= vxy[0];
                        col_vel[j][x][1] -= vxy[1];
                        col_vel[k][y][0] += vxy[0];
                        col_vel[k][y][1] += vxy[1];
                        
                    }
                  }

                }
            }
        }
        
        // Calculating alpha based on prescribed curvature
        for (int i = 0; i < N_beads; ++i) {
            alpha[i] = L0b * sin(phi[j] + TWO_PI * (freq[j] * t + i * (2.0 / N_beads)));
            cos_alpha[i] = cos(alpha[i]);
            sin_alpha[i] = sin(alpha[i]);
        }

        // Calculating tangents and Hookean forces
        for (int i = 0; i < N_beads - 1; ++i) {
            ta[i][0] = positions[j][i + 1][0] - positions[j][i][0];
            ta[i][1] = positions[j][i + 1][1] - positions[j][i][1];
            double length = sqrt(ta[i][0] * ta[i][0] + ta[i][1] * ta[i][1]);
            double Kdl_n = K * (length - L0) / length;
            double f1 = Kdl_n * ta[i][0];
            double f2 = Kdl_n * ta[i][1];
            F1[i][0] += f1;
            F1[i][1] += f2;
            F1[i + 1][0] -= f1;
            F1[i + 1][1] -= f2;
        }

        // Bending force calculations
        F2[0][0] -= ta[0][0];
        F2[0][1] -= ta[0][1];
        F2[1][0] += ta[0][0];
        F2[1][1] += ta[0][1];

        for (int i = 0; i < N_beads - 2; ++i) {
            double rtix = cos_alpha[i] * ta[i][0] - sin_alpha[i] * ta[i][1];
            double rtiy = cos_alpha[i] * ta[i][1] + sin_alpha[i] * ta[i][0];
            F2[i + 2][0] += rtix;
            F2[i + 2][1] += rtiy;
            F2[i + 1][0] -= rtix;
            F2[i + 1][1] -= rtiy;

            double rttiux = 2 * ta[i][0] - cos_alpha[i] * ta[i + 1][0] - sin_alpha[i] * ta[i + 1][1];
            double rttiuy = 2 * ta[i][1] - cos_alpha[i] * ta[i + 1][1] + sin_alpha[i] * ta[i + 1][0];
            F2[i][0] += rttiux;
            F2[i][1] += rttiuy;
            F2[i + 1][0] -= rttiux;
            F2[i + 1][1] -= rttiuy;
        }

        F2[N_beads - 2][0] += ta[N_beads - 2][0];
        F2[N_beads - 2][1] += ta[N_beads - 2][1];
        F2[N_beads - 1][0] -= ta[N_beads - 2][0];
        F2[N_beads - 1][1] -= ta[N_beads - 2][1];

        // Scaling F2 by bending constant
        for (int i = 0; i < N_beads; ++i) {
            F2[i][0] *= K_bending;
            F2[i][1] *= K_bending;
        }

        // Aggregating Hookean and bending forces
        for (int i = 0; i < N_beads; ++i) {
            ac[j][i][0] = (F1[i][0] + F2[i][0]) / m;
            ac[j][i][1] = (F1[i][1] + F2[i][1]) / m;
        }

        
    }
    for (int j=0; j<N_systems ; ++j){   // adding velocity changes during collisions to the velocities
        for (int i = 0; i < N_beads; ++i){
            velocities[j][i][0] += col_vel[j][i][0];
            velocities[j][i][1] += col_vel[j][i][1];
            
        }
    }
    
}

// Handles collision with the moving bounding wall
void bounce_forward_boundary(double* pos,double* vel){
    // Calculate the distance of the bead from the center
    double distance_from_center = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

    // Check if the bead is outside the boundary
    if (distance_from_center > R) {
        pos[0] -= vel[0]*dt_half;
        pos[1] -= vel[1]*dt_half;
        double rad = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

        // Project the bead back to the boundary
        double fac = R / rad ;
        pos[0] *= fac ; pos[1] *= fac ; // Scale position to stay on the boundary

        // Reflect the velocity along the normal direction to "bounce back" but keep the tangential component
        // Using the formula for velocity reflection: v = v - 2 * (v ⋅ n) * n
        // note: pos now is equal to normal vector
        double dot_prod = vel[0]*pos[0] + vel[1]*pos[1] ;
        fac = (( 2 * dot_prod/ R ) - dRdt) / R ;
        vel[0] -=   fac  * pos[0] ; vel[1] -= fac * pos[1] ;
        pos[0] += vel[0]*dt_half ; pos[1] += vel[1]*dt_half ;
    }
}


// uses propel function to calculate forces and handle collisions and then
// updates positions,velocites and increments time and counter.
// After updation,implements andersen thermostat
void update_positions(){
    // calculate new accelerations
    propel() ;
    // calculate new radius of system
    R += dRdt * dt;
    dRdt += d2Rdt*dt ;
    // Update velocities and then positions
    for (int j = 0; j<N_systems ; ++j){
        for (int i = 0; i<N_beads ; ++i){
            velocities[j][i][0] += ac[j][i][0] * dt ;
            velocities[j][i][1] += ac[j][i][1] * dt ;

            positions[j][i][0] += (velocities[j][i][0] + ac[j][i][0]*dt_half) * dt ;
            positions[j][i][1] += (velocities[j][i][1] + ac[j][i][1]*dt_half) * dt ;
            bounce_forward_boundary(positions[j][i],velocities[j][i]) ;
        }
        if ( (counter%100 == 0) && (((double)rand()/RAND_MAX)<nudt) ){
        //andersen thermostat step
        for (int i = 0; i<N_beads ; ++i){
                double normal[2] = {0};
                rand_normal(normal);
                velocities[j][i][0] = normal[0];
                velocities[j][i][1] = normal[1];
            }
        }
    }

    t += dt ;
    counter+=1;
}
void update_positions_fix(){
    // calculate new accelerations
    propel() ;
    //same function as update_positions(), but without radius changing
    // Update velocities and then positions
    for (int j = 0; j<N_systems ; ++j){
        for (int i = 0; i<N_beads ; ++i){
            velocities[j][i][0] += ac[j][i][0] * dt ;
            velocities[j][i][1] += ac[j][i][1] * dt ;

            positions[j][i][0] += (velocities[j][i][0] + ac[j][i][0]*dt_half) * dt ;
            positions[j][i][1] += (velocities[j][i][1] + ac[j][i][1]*dt_half) * dt ;
            bounce_forward_boundary(positions[j][i],velocities[j][i]) ;
        }
        if ( (counter%100 == 0) && (((double)rand()/RAND_MAX)<nudt) ){
        //andersen thermostat step
        //scope for further experimentation here
        for (int i = 0; i<N_beads ; ++i){
                double normal[2] = {0};
                rand_normal(normal);
                velocities[j][i][0] = normal[0];
                velocities[j][i][1] = normal[1];
            }
        }
    }

    t += dt ;
    counter+=1;
}

// reads values intially from the data file whose name is provided as arguement
void read_from_backup(char* f_name) {
    input = fopen(f_name, "r");
    // reading values from file
    char temp_str[1024];int s = 0;
    int x = 0 ;
    char c;
    while (1) {   // this loop stores time and R value
        char c = fgetc(input);
        if (c==',' || c=='\n'){
            temp_str[s] = '\0';
            s = 0;
            if (x==0){
                t = strtof(temp_str,NULL);  //time
                x = 1 ;
            }
            else if(x==1){
                counter = atoi(temp_str);   // counter
                x = 2 ;
            }
            else{
                R = strtod(temp_str,NULL );     //radius
                x = 0 ;
                break;
            }
        }
        else temp_str[s++] = c;
    }
    for (int i = 0; i < N_systems;++i) {   // this loop reads frequencies
        while (1) {   // this loop reads one complete number
            char c = fgetc(input);
            if (c==',' || c=='\n'){
                temp_str[s] = '\0';
                s = 0;
                freq[i] = strtod(temp_str,NULL );
                break;
            }
            else temp_str[s++] = c;
        }
    }
    for (int i = 0; i < N_systems;++i) {   // this loop reads phases
        while (1) {   // this loop reads one complete number
            char c = fgetc(input);
            if (c==',' || c=='\n'){
                temp_str[s] = '\0';
                s = 0;
                phi[i] = strtod(temp_str,NULL );
                break;
            }
            else temp_str[s++] = c;
        }
    }
    //reading positions, velocities and accelerations
    for (int i = 0; i < N_systems; ++i) {
      for (int j = 0; j < N_beads ; ++j) {
          while (1) {   // this loop reads one complete number
              char c = fgetc(input);
              if (c==',' || c=='\n'){
                  temp_str[s] = '\0';
                  positions[i][j][x] = strtod(temp_str,NULL);
                  x = (x+1) % 2 ;   //x is flipped
                  s = 0;
                  if (x==0)  break;
              }
              else temp_str[s++] = c;
          }
      }
      for (int j = 0; j < N_beads ; ++j) {
          while (1) {   // this loop reads one complete number
              char c = fgetc(input);
              if (c==',' || c=='\n'){
                  temp_str[s] = '\0';
                  velocities[i][j][x] = strtod(temp_str,NULL);
                  x = (x+1) % 2 ;
                  s = 0;
                  if (x==0)  break;
              }
              else temp_str[s++] = c;
          }
      }
      for (int j = 0; j < N_beads ; ++j) {
          while (1) {   // this loop reads one complete number
              char c = fgetc(input);
              if (c==',' || c=='\n'){
                  temp_str[s] = '\0';
                  ac[i][j][x] = strtod(temp_str,NULL);
                  x = (x+1) % 2 ;
                  s = 0;
                  if (x==0)  break;
              }
              else temp_str[s++] = c;
          }
      }
    }
    fclose(input); // input file closed
}


// writes the final values,after all iterations, to the file whose filename is provided as arguement
void write_to_backup(char *f_name) {
    // saving output to output file
    output = fopen(f_name, "w");
    fprintf(output, "%.16f,%d,%.16f\n", t,counter, R); // storing 't,R\n'
    //
    for (int i = 0 ; i< (N_systems-1) ; ++i){   //storing comma separated frequencies
        fprintf(output,"%.16f,",freq[i]);
    }fprintf(output, "%.16f\n", freq[N_systems - 1]);

    for (int i = 0 ; i< (N_systems-1) ; ++i){   //storing comma separated phases
        fprintf(output,"%.16f,",phi[i]);
    }fprintf(output, "%.16f\n", phi[N_systems - 1]);

    for (int i = 0; i < N_systems ; ++i){   // storing positions, velocities and accelerations
        for (int j = 0; j < N_beads-1 ; ++j) {    //storing positions
            fprintf(output,"%.16f,%.16f,",positions[i][j][0],positions[i][j][1]);
        }
        fprintf(output, "%.16f,%.16f\n", positions[i][N_beads - 1][0],positions[i][N_beads - 1][1]);
        for (int j = 0; j < N_beads-1 ; ++j) {    //storing velocities
            fprintf(output,"%.16f,%.16f,",velocities[i][j][0],velocities[i][j][1]);
        }
        fprintf(output,"%.16f,%.16f\n",velocities[i][N_beads-1][0],velocities[i][N_beads-1][1]);
        for (int j = 0; j < N_beads-1 ; ++j) {    //storing accelerations
            fprintf(output,"%.16f,%.16f,",ac[i][j][0],ac[i][j][1]);
        }
        fprintf(output,"%.16f,%.16f\n",ac[i][N_beads-1][0],ac[i][N_beads-1][1]);
    }

    //closing output file
    fclose(output);
}

// This function takes a sample(that includes time,counter,radius,positions,velocites) and appends it to the end of filename provided as argument
void sample(char *f_name){  // samples t,R,counter,positions,velocities and appends to specified filename
    // saving output to output file
    output = fopen(f_name, "a");
    fprintf(output, "%.16f,%d,%.16f\n", t,counter, R); // storing 't,R\n'


    for (int i = 0; i < N_systems ; ++i){   // storing positions, velocities and accelerations
        for (int j = 0; j < N_beads-1 ; ++j) {    //storing positions
            fprintf(output,"%.16f,%.16f,",positions[i][j][0],positions[i][j][1]);
        }
        fprintf(output, "%.16f,%.16f\n", positions[i][N_beads - 1][0],positions[i][N_beads - 1][1]);
        // for (int j = 0; j < N_beads-1 ; ++j) {    //storing velocities
        //     fprintf(output,"%.16f,%.16f,",velocities[i][j][0],velocities[i][j][1]);
        // }
        // fprintf(output,"%.16f,%.16f\n",velocities[i][N_beads-1][0],velocities[i][N_beads-1][1]);
    }

    //closing output file
    fclose(output);
}


// main program involving reading input values, running simulation, writing output values.
int main(int argc, char* argv[]) {//program behavior depends on number of arguments
    iter = atoi(argv[1]); // number of iterations of simulation
    read_from_backup(argv[2]);//reading values
    if(argc==4){//reducing radius to reach initial state of high density
        
        // running simulation
        for (int i = 0 ; i < iter ; ++i){
            update_positions();
        }
        // saving output to backup
        write_to_backup(argv[3]);
        //an example for a command line is
        //propagate.exe 10000 data.txt out.txt
    }
    
    // else{//sampling with constant radius, with linear interval 
    //     int sfreq = atoi(argv[4]);
    //     for(int i=0;i< iter;++i){
    //         update_positions_fix();
    //         if(i % sfreq ==0){ //dumping sample
    //             sample(argv[3]);
    //         }
    //     }
    //     //saving final configuration once again(use it for velocity distr)
    //     write_to_backup(argv[5]);
    // }
    else { // sampling with log-spaced intervals using n_per_decade, at fixed radius
        int n_per_decade = atoi(argv[4]);  // number of samples per decade
        double log_base = log(10.0);       // natural log of 10
        double step_size = log_base / n_per_decade;
        int next_sample = 0;
        int sample_index = atoi(argv[6]); //use for continuing with a previous sample, for new sample just use 1
        // Compute next sample time: floor(exp(step_size * sample_index))
        next_sample = (int)floor(exp(step_size * sample_index));
        for (int i = 0; i < iter; ++i) {
            update_positions_fix();
            if (i == next_sample) {// sample at this log-spaced instant
                sample(argv[3]);  
                sample_index++;
                while ((int)floor(exp(step_size * sample_index))==next_sample){
                    sample_index++;
                }
                next_sample = (int)floor(exp(step_size * sample_index));
                if (next_sample >= iter) break;  // ensure we stay within bounds
            }
        }
        // save final configuration (for velocity distribution, etc.)
        write_to_backup(argv[5]);
        //example for command line is 
        //propagate.exe 100000 out.txt sample.txt 4 final.txt 10
    }
    
    
}

