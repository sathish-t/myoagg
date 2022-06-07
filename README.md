# Simulation of cytokinetic rings in isolated fission yeast protoplasts

Simulation of a fully three-dimensional coarse-grained mathematical model of the
fission yeast cytokinetic ring is presented here. The goal is to observe a stable,
tense ring under normal conditions and to investigate ring phenotypes when component
turnover is absent. A brief overview of the simulation follows. For complete details,
refer to the Results and the Methods section of
https://www.biorxiv.org/content/10.1101/2021.03.18.436017v1 or later version(s) if
available. Please cite this preprint or a peer-reviewed, published version whichever
is later if you use the software in your research.

The simulation has four components - the myosin-II molecule Myo2, actin filaments,
the formin Cdc12, and the actin crosslinker alpha-actinin. The first three components
are represented as beads whereas the last is represented as a simple spring. Myo2 and Cdc12
molecules are anchored to the plasma membrane, represented as the inner surface of a cylinder
of fixed radius. They can slide along the membrane, maintaining fixed distances from it,
and experience a drag force from the membrane. One bead of myosin represents an anchored
cluster of 16 heads. One actin bead represents a 100 nm stretch of actin.

Force producing interactions in the simulation are: a myosin cluster binding actin filament(s)
with a capture force, myosin clusters exerting pulling forces on actin filaments along their length
and dependent on their relative velocity, actin crosslinking forces, actin bending forces, drag
forces from the cytosol, lateral and perpendicular forces from the membrane on anchored components,
excluded volume forces between myosins, and restoring forces that push straying components back into the
interior of the cylinder. Components evolve according to first-order dynamics.

On and off rates of Myo2, Cdc12, alpha-actinin, the actin synthesis rate, and the actin severing
rate produced by cofilin are all tuned to faithfully reproduce experimentally measured
amounts and turnover times. After 10 minutes of ring simulation, dissociation rates can be
greatly reduced and association rates abolished to investigate ring phenotypes in the
absence of turnover.

A brief overview of the MATLAB files follows. Simulation can be run using the Slurm workload
manager on a computing cluster using runv9Terre.sh, or standalone simulations can be initiated
using the command task_yeti(num) where num is supplied as the seed to the random number generator.

Filename | Description 
--- | ---
bm_nbr.m | Identify which myosin cluster interacts with which actin bead
get_a.m | Construct one part of the drag matrix for first-order dynamics
get_force.m | Calculate forces on all simulation components
get_m.m | Construct the other part of the drag matrix for first-order dynamics
get_pcpq_xy.m | Calculate derivatives of constraint equations
get_rdiff.m | Obtain vectors along actin filament lengths
initial_circle_rand.m | Initialize the ring
parameters_v7.m | Set parameters for the simulation
ptcap.m | Adjust which actin bead interacts with which myosin cluster
rk1.m | Organizing script that executes each part of the simulation in order
run_pos_z_cyl.m | Set some more simulation parameters and initialize the ring, discarding defective initial conditions
runv9Terre.sh | Shell script to execute the simulation using a Slurm workload manager
stall_force_pf.m | Calculate forces exerted by myosin clusters
task_yeti_v9.m | Called by runv9Terre. Main coordinating script that executes many cycles of rk1 and saves output data from each cycle.
update.m | Adjust ring component numbers
velocity.m | Obtain constraint forces and component velocities

Henceforth are the analysis scripts in the analysis/ folder.

Filename | Description 
--- | ---
tens_vs_time.m | Calculate tension vs time for a set of runs
tension_circ.m | Helper script of tens_vs_time.m
get_segtens.m | Helper script of tens_vs_time.m
get_tens_ring.m | Helper script of tens_vs_time.m
runv9TerreTiffCreate.sh | Shell script to generate TIFF files from simulation outputs using a slurm manager
ringplot_prepare_tiff.m | Main coordinating script to generate TIFF files
fluo_emu_myo.m | Generate one frame of a TIFF file
corr_len_vs_time_sim.m | Calculate correlation function for a set of runs
corrcalcCirc.m | Helper script of corr_len_vs_time_sim.m
corr_len_vs_time_expt.m | Calculate correlation function for experimental kymographs
corrcalc.m | Helper script of corr_len_vs_time_expt.m
track_aggregates.m | Script to plot mean distance b/w aggregates vs time, both simulation and experiment
st_plot_ring_2d_3feb21.m | Plot the simulated ring
make_kymo.ijm | Fiji macro to make kymographs.
prom2.roi | Region of interest (ROI) used by Fiji macro to make kymographs.

## Notes
1. Use make_kymo.ijm with prom2.roi and fake confocal images created from simulation data to generate kymographs.
2. Use make_kymo.ijm with the respective ROI file and the respective real confocal image from experiment to generate
real experimental kymographs. The width parameter in make_kymo.ijm may have to be adjusted, but it is usually in the range
10-18.