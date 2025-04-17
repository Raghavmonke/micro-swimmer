//This file is used to perform various calculations on the sampled data, behaviour is dependent on arguments provided
//1. We can calculate MSD
//2. We can calculate the intermediate scattering function for a given wavevector value
//3. Given the final state of the system (not the sampled data), make a csv with vx,vy pairs (can plot velocity distributions with it)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#define N_beads 20  //Number of beads used in taylor line
#define N_systems 120 //Number of systems
#define PI 3.141592653589793
#define VALUES_PER_LINE (N_beads * 2)  // x and y per bead
#define LINE_BUFFER 2048

double xref[N_systems][N_beads], yref[N_systems][N_beads];
double x_ref[N_systems][N_beads];
double y_ref[N_systems][N_beads];
void compute_scattering_function(const char *filename, double q) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open input file");
        return;
    }

    char output_filename[128];

    // Round q to 2 decimal places and convert to int (e.g. 6.28 → 628)
    int q_int = (int)(q * 100 + 0.5);
    int q_whole = q_int / 100;
    int q_frac  = q_int % 100;

    // Create filename like scattering_q_6p28.csv
    snprintf(output_filename, sizeof(output_filename),
            "scattering_q_%dp%02d.csv", q_whole, q_frac);

    FILE *out = fopen(output_filename, "w");
    if (!out) {
        perror("Failed to open output file");
        fclose(fp);
        return;
    }

    char line[LINE_BUFFER];
    double time;
    int timestep = 0;

    fprintf(out, "Time,f(q,t)\n");

    while (fgets(line, sizeof(line), fp)) {
        sscanf(line, "%lf", &time);  // First line of the block

        double complex fq_sum = 0.0 +I*0.0;

        for (int i = 0; i < N_systems; ++i) {
            if (!fgets(line, sizeof(line), fp)) break;

            char *token = strtok(line, ",");
            for (int j = 0; j < N_beads; ++j) {
                if (token == NULL) break;
                double x = atof(token);
                token = strtok(NULL, ",");
                if (token == NULL) break;
                double y = atof(token);

                if (timestep == 0) {
                    x_ref[i][j] = x;
                    y_ref[i][j] = y;
                } else {
                    double dx = x - x_ref[i][j];
                    double dy = y - y_ref[i][j];

                    // Assuming q vector is along (1, 1)
                    double dot = (q / sqrt(2)) * (dx + dy);
                    fq_sum += cexp(I * dot);
                }

                if (j < N_beads - 1) token = strtok(NULL, ",");
            }
        }

        if (timestep > 0) {
            double fq = creal(fq_sum) / (N_systems * N_beads);
            fprintf(out, "%.8f, %.8f\n", time, fq);
        }

        timestep++;
    }

    fclose(fp);
    fclose(out);
    printf("Scattering function written to scattering_output.csv\n");
}

void compute_msd(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Failed to open file");
        return;
    }

    FILE *out = fopen("msd_output.csv", "w");
    if (!out) {
        perror("Failed to open output file");
        fclose(fp);
        return;
    }

    char line[LINE_BUFFER];
    double time = 0.0;
    int timestep = 0;
    int total = N_systems * N_beads;
    fprintf(out, "Time,MSD\n");

    while (fgets(line, sizeof(line), fp)) {
        // First line of block: time info
        sscanf(line, "%lf", &time);

        double msd = 0.0;

        for (int i = 0; i < N_systems; ++i) {
            if (!fgets(line, sizeof(line), fp)) break;

            char *token = strtok(line, ",");
            for (int j = 0; j < N_beads; ++j) {
                if (token == NULL) break;
                double x = atof(token);
                token = strtok(NULL, ",");
                if (token == NULL) break;
                double y = atof(token);

                if (timestep == 0) {
                    xref[i][j] = x;
                    yref[i][j] = y;
                } else {
                    double dx = x - xref[i][j];
                    double dy = y - yref[i][j];
                    msd += dx * dx + dy * dy;
                }

                if (j < N_beads - 1) token = strtok(NULL, ",");
            }
        }

        if (timestep > 0) {
            
            fprintf(out, "%.8f, %.8f\n", time, msd / total);
        }

        timestep++;
    }

    fclose(fp);
    fclose(out);
    printf("MSD written to msd_output.csv\n");
}
void extract_velocities(const char *input_file, const char *output_csv) {
    
    //extracts velocities from final state and dumps into csv as vx,vy pairs
    FILE *in = fopen(input_file, "r");
    if (!in) {
        perror("Error opening input file");
        return;
    }

    FILE *out = fopen(output_csv, "w");
    if (!out) {
        perror("Error opening output file");
        fclose(in);
        return;
    }

    char line[8192];

    // Skip first 3 header lines
    for (int i = 0; i < 3; ++i) {
        if (!fgets(line, sizeof(line), in)) {
            fprintf(stderr, "Unexpected end of file while skipping headers.\n");
            fclose(in);
            fclose(out);
            return;
        }
    }

    int line_counter = 0;
    while (fgets(line, sizeof(line), in)) {
        line_counter++;

        // Only process every 3rd line (velocity lines)
        if (line_counter % 3 == 2) {
            // Parse the velocity line
            char *token = strtok(line, ",");
            int bead_counter = 0;

            while (token && bead_counter < N_beads) {
                double vx = atof(token);
                token = strtok(NULL, ",");
                if (!token) break;

                double vy = atof(token);
                fprintf(out, "%.6f,%.6f\n", vx, vy);
                bead_counter++;
                token = strtok(NULL, ",");
            }
        }
    }

    fclose(in);
    fclose(out);
    printf("Velocities extracted and saved to %s\n", output_csv);
}
int main(int argc, char* argv[]){
    if(argc == 3){
        extract_velocities(argv[1],argv[2]);
        //example
        //<executable> <finalstatefilename> <vel csv>
    }
    else{compute_msd("sample.txt");
    compute_scattering_function("sample.txt",PI*10);//change value of q or filename here. 
    //example
    //<executable>
}
    //note: time dumped in the msd and scattering files does not start from 0.
    //Time taken to reach the required volume fraction needs to be manualy subtracted (significant for log plots)
}