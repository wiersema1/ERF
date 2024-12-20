# ------------------  INPUTS TO MAIN PROGRAM  -------------------
max_step = 10

amrex.fpe_trap_invalid = 1

fabarray.mfiter_tile_size = 2048 1024 2048

# PROBLEM SIZE & GEOMETRY
geometry.prob_lo     = -15000.   0.    0.
geometry.prob_hi     =  15000. 100. 24000.
amr.n_cell           =    300    1    240    # dx=dy=dz=100 m

geometry.is_periodic = 0 1 0

xlo.type = "Open"
xhi.type = "Open"
zlo.type = "SlipWall"
zhi.type = "HO_Outflow"

# TIME STEP CONTROL
erf.fixed_dt       = 0.25      # fixed time step [s] -- Straka et al 1993
erf.fixed_fast_dt  = 0.125     # fixed time step [s] -- Straka et al 1993

# DIAGNOSTICS & VERBOSITY
erf.sum_interval   = 1       # timesteps between computing mass
erf.v              = 1       # verbosity in ERF.cpp
amr.v              = 1       # verbosity in Amr.cpp

# REFINEMENT / REGRIDDING
amr.max_level       = 0       # maximum level number allowed

# CHECKPOINT FILES
amr.check_file      = chk        # root name of checkpoint file
amr.check_int       = 1000       # number of timesteps between checkpoints
#amr.restart         = chk09000

# PLOTFILES
erf.plot_file_1         = plt        # root name of plotfile
erf.plot_int_1          = 100         # number of timesteps between plotfiles
erf.plot_vars_1         = density rhotheta rhoQ1 rhoQ2 rhoQ3 x_velocity y_velocity z_velocity pressure theta temp qv qc qrain rain_accum pert_dens

# SOLVER CHOICE
erf.use_gravity = true
erf.buoyancy_type = 1
erf.use_coriolis = false

#erf.les_type = "Smagorinsky"
#erf.Cs       = 0.25
erf.les_type = "None"

#
# diffusion coefficient from Straka, K = 75 m^2/s
#
erf.molec_diff_type   = "ConstantAlpha"
erf.dynamic_viscosity = 200.0 # [kg/(m-s)]
erf.alpha_T           = 200.0 # [m^2/s]
erf.alpha_C           = 200.0

erf.moisture_model = "Kessler"
erf.use_moist_background = true

# PROBLEM PARAMETERS (optional)
prob.z_tr = 12000.0
prob.height = 1200.0
prob.theta_0 = 300.0
prob.theta_tr = 343.0
prob.T_tr = 213.0
prob.x_c = 0.0
prob.z_c = 2000.0
prob.x_r = 10000.0
prob.z_r = 1500.0
prob.theta_c = 3.0
