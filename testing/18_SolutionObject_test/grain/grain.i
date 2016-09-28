# This simulation predicts GB migration of a 2D copper polycrystal with 15 grains
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 12 # Number of elements in the x-direction
  ny = 12 # Number of elements in the y-direction
  nz = 0 # Number of elements in the z-direction
  xmax = 1000 # maximum x-coordinate of the mesh
  ymax = 1000 # maximum y-coordinate of the mesh
  zmax = 0
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 2 # Initial uniform refinement of the mesh
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 4 # Number of grains
  var_name_base = gr # Base name of grains
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  active = 'PolycrystalVariables c'
  [./PolycrystalVariables]
    # Custom action that created all of the grain variables
    order = FIRST # element type used by each grain variable
    family = LAGRANGE
  [../]
  [./c]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalVoronoiIC]
      grain_num = 4
      advanced_op_assignment = true
    [../]
  [../]
[]

[AuxVariables]
  # active = ''
  # Dependent variables
  active = 'bnds unique_grains ghost_regions halos var_indices'
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
    order = FIRST
    family = LAGRANGE
  [../]
  [./sigma11_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma22_aux]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ghost_regions]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./halos]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  active = 'PolycrystalKernel bndsDiff time_c'
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
  [../]
  [./bndsDiff]
    type = MatDiffusion
    variable = c
    args = bnds
    conc = c
  [../]
  [./time_c]
    type = TimeDerivative
    variable = c
  [../]
  [./time_x]
    type = TimeDerivative
    variable = disp_x
  [../]
  [./time_y]
    type = TimeDerivative
    variable = disp_y
  [../]
[]

[AuxKernels]
  # active = ''
  # AuxKernel block, defining the equations used to calculate the auxvars
  active = 'unique_grains bnds_aux ghosted_entities halos var_indices'
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./matl_sigma11]
    type = RankTwoAux
    variable = sigma11_aux
    rank_two_tensor = stress
    index_j = 0
    index_i = 0
  [../]
  [./matl_sigma22]
    type = RankTwoAux
    variable = sigma22_aux
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
  [../]
  [./halos]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
  [../]
[]

[BCs]
  # Boundary Condition block
  active = 'left'
  [./left]
    type = DirichletBC
    variable = c
    boundary = left
    value = 1
  [../]
  [./right]
    type = DirichletBC
    variable = c
    boundary = right
    value = 0
  [../]
  [./fixed_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0
  [../]
  [./fixed_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./top_disp]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 100
  [../]
[]

[Materials]
  active = 'CuGrGr diffusionMat'
  [./CuGrGr]
    # Material properties
    type = GBEvolution # Quantitative material properties for copper grain growth.  Dimensions are nm and ns
    GBmob0 = 2.5e-6 # Mobility prefactor for Cu from Schonfelder1997
    GBenergy = 0.708 # GB energy for Cu from Schonfelder1997
    Q = 0.23 # Activation energy for grain growth from Schonfelder 1997
    T = 450 # K   #Constant temperature of the simulation (for mobility calculation)
    wGB = 13 # nm      #Width of the diffuse GB
  [../]
  [./diffusionMat]
    type = ParsedMaterial
    function = if(bnds>0.99,0.00001,500)
    f_name = D
    args = bnds
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '15 15'
    fill_method = symmetric_isotropic
    elasticity_tensor_prefactor = 2
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./dependence]
    type = DerivativeParsedMaterial
    function = if(bnds>0.99,1,0.9)
    f_name = var_dep
    args = bnds
    derivative_order = 2
  [../]
  [./eigenstrain]
    type = ComputeVariableEigenstrain
    eigen_base = '1 1 1 0 0 0'
    args = bnds
    prefactor = var_dep
  [../]
  [./conc_elas]
    type = ComputeConcentrationDependentElasticityTensor
    c = c
    C1_ijkl = '7 7'
    C0_ijkl = '15 15'
    fill_method1 = symmetric_isotropic
    fill_method0 = symmetric_isotropic
  [../]
[]

[Postprocessors]
  # Scalar postprocessors
  active = 'dt grain_tracker'
  [./ngrains]
    # Counts the number of grains in the polycrystal
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./grain_tracker]
    type = GrainTracker
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre    boomeramg      101                ds'
  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_abs_tol = 1e-11 # Relative tolerance for nonlienar solves
  nl_rel_tol = 1e-8 # Absolute tolerance for nonlienar solves
  start_time = 0.0
  end_time = 10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]
  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    max_h_level = 3 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  [./console]
    type = Console
    max_rows = 20
  [../]
[]

