import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Circle
#This code takes the sampled data and creates and saves snapshots of the system.
#use appropriate ffmpeg commands to make a movie out of it then
# Constants
N_systems = 120
N_beads = 20
r = 0.1  # bead radius

def plot_frame(t, R, positions, frame_number, out_dir="frames"):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(-R - 5, R + 5)
    ax.set_ylim(-R - 5, R + 5)
    ax.set_aspect('equal')
    ax.grid(True)

    # draw system boundary
    ax.add_patch(Circle((0, 0), R, color='black', fill=False))

    # plot swimmers
    for j in range(N_systems):
        x = [positions[j][i][0] for i in range(N_beads)]
        y = [positions[j][i][1] for i in range(N_beads)]
        ax.plot(x, y, 'b-')
        for i in range(N_beads):
            ax.add_patch(Circle((x[i], y[i]), r, color='r', fill=True))

    ax.set_title(f"t = {t:.3f}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    os.makedirs(out_dir, exist_ok=True)
    filename = os.path.join(out_dir, f"frame_{frame_number:05d}.png")
    plt.savefig(filename, dpi=200)
    plt.close(fig)

def main():
    input_file = "sample.txt"
    output_dir = "frames"
    frame_number = 0

    with open(input_file, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break  # EOF

            # read time and R
            tokens = header.strip().split(",")
            if len(tokens) < 3:
                break
            t = float(tokens[0])
            R = float(tokens[2])

            positions = []

            for _ in range(N_systems):
                line = f.readline()
                if not line:
                    break
                vals = list(map(float, line.strip().split(",")))
                body = []
                for j in range(N_beads):
                    x = vals[2*j]
                    y = vals[2*j + 1]
                    body.append((x, y))
                positions.append(body)

            plot_frame(t, R, positions, frame_number, output_dir)
            frame_number += 1

if __name__ == "__main__":
    main()
