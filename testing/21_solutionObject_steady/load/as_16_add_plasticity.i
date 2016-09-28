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
  elem_type = TRI3 # Type of elements used in the mesh
  uniform_refine = 2 # Initial uniform refinement of the mesh
  skip_partitioning = true
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 4 # Number of grains
  var_name_base = gr # Base name of grains
[]

[MeshModifiers]
  [./seed]
    type = BoundingBoxNodeSet
    top_right = '10 650 0'
    new_boundary = seed
    bottom_left = '-1 550 0'
  [../]
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./c]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[AuxVariables]
  # active = ''
  # Dependent variables
  active = 'failed sigma11_aux sigma22_aux plastic_strain von_mises cAux bndsSol'
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
  [./failed]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./von_mises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./cAux]
  [../]
  [./bndsSol]
  [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  active = 'TensorMechanics'
  [./bndsDiff]
    type = MatDiffusion
    variable = c
    args = bnds
    conc = c
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    use_displaced_mesh = true
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
  active = 'matl_sigma11 eps_kernel failed matl_sigma22 von_mises_kernel cSolAux bndsSolAux'
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
  [./failed]
    type = ParsedAux
    variable = failed
    function = if(von_mises>0.5,if(bnds>0.8,1,0),0)
    args = 'von_mises bnds'
    execute_on = timestep_end
  [../]
  [./von_mises_kernel]
    type = RankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./eps_kernel]
    type = RankTwoScalarAux
    variable = plastic_strain
    scalar_type = EquivalentPlasticStrain
    rank_two_tensor = stress
    execute_on = timestep_end
  [../]
  [./cSolAux]
    type = SolutionAux
    variable = cAux
    from_variable = c
    solution = cSoln
    execute_on = initial
  [../]
  [./bndsSolAux]
    type = SolutionAux
    variable = bndsSol
    from_variable = bnds
    solution = bndsSoln
    execute_on = initial
  [../]
[]

[BCs]
  # Boundary Condition block
  active = 'Pressure fixed_x fixed_y'
  [./left]
    type = FunctionDirichletBC
    variable = c
    boundary = seed
    function = if(t<10,t/10,1)
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
    boundary = right
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
  [./Pressure]
    [./top_pressure]
      boundary = top
      factor = -10
      disp_y = disp_y
      disp_x = disp_x
    [../]
  [../]
[]

[Materials]
  active = 'conc_elas finiteStrain diffusionMat isotropic_plasticity_recompute radial_return_stress'
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
    function = if(bnds>0.8,grainDiff,boundDiff)
    f_name = D
    args = bnds
    constant_expressions = '0.00001 500'
    constant_names = 'grainDiff boundDiff'
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
  [./conc_elas]
    type = ComputeConcentrationDependentElasticityTensor
    c = c
    C1_ijkl = '7 7'
    C0_ijkl = '15 15'
    fill_method1 = symmetric_isotropic
    fill_method0 = symmetric_isotropic
  [../]
  [./finiteStrain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y'
  [../]
  [./isotropic_plasticity_recompute]
    type = RecomputeRadialReturnIsotropicPlasticity
    relative_tolerance = 1e-10
    yield_stress = 10
    hardening_constant = 5
    absolute_tolerance = 1e-12
    output_iteration_info_on_error = true
    max_iterations = 50
  [../]
  [./radial_return_stress]
    type = ComputeReturnMappingStress
    return_mapping_models = isotropic_plasticity_recompute
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
  [./sigma11_post]
    type = ElementIntegralVariablePostprocessor
    variable = sigma11_aux
  [../]
[]

[UserObjects]
  [./yield]
    type = TensorMechanicsHardeningConstant
    value = 0.5
  [../]
  [./isoPlasticity]
    type = TensorMechanicsPlasticIsotropicSD
    c = -0.779422863
    b = -0.2
    use_custom_cto = false
    yield_function_tolerance = 1e-5
    yield_strength = yield
    internal_constraint_tolerance = 1e-9
    use_custom_returnMap = false
  [../]
  [./bndsSoln]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = bnds
    system = aux0
    mesh = ../grain/peacock_run_tmp_out.e-s125
    execute_on = initial
  [../]
  [./cSoln]
    type = SolutionUserObject
    timestep = LATEST
    system_variables = c
    mesh = ../grain/peacock_run_tmp_out.e-s125
    execute_on = initial
  [../]
[]

[Executioner]
  # Preconditioned JFNK (default)
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  dt = 10
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre    boomeramg      31                ds'
  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 10 # Max number of nonlinear iterations
  nl_abs_tol = 1e-11 # Relative tolerance for nonlienar solves
  nl_rel_tol = 1e-8 # Absolute tolerance for nonlienar solves
  start_time = 0.0
  end_time = 1000
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1 # Initial time step.  In this simulation it changes.
    growth_factor = 1.1
  [../]
  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 3 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    max_h_level = 5
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

