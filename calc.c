#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define VALUES_PER_LINE (N_beads * 2)  // x and y per bead
#define LINE_BUFFER 2048
#define N_beads 20  //Number of beads used in taylor line
#define N_systems 120 //Number of systems
double xref[N_systems][N_beads], yref[N_systems][N_beads];

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
int main(int argc, char* argv[]){
    compute_msd("sample.txt");
}