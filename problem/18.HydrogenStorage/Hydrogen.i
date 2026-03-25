# H2 injection in a brine aquifer
# basic case is from cpgr/hystre github

# Injection schedule
# Inject 10,000 kg H2 over 5 days, then store for 30 days, then produce for 10 days
inject_mass = 1E4
inject_time = 5
store_time = 30
produce_mass = 1E4
produce_time = 10
end_time = '${fparse inject_time + store_time + produce_time}'

start_injection = 0
end_injection = '${fparse start_injection + inject_time}'

start_production = '${fparse end_injection + store_time}'
end_production = '${fparse start_production + produce_time}'

# Time in seconds (used in inputs below)
inject_time_s = '${fparse inject_time * 86400}'
produce_time_s = '${fparse produce_time * 86400}'
start_injection_s = '${fparse start_injection * 86400}'
end_injection_s = '${fparse end_injection * 86400}'
start_production_s = '${fparse start_production * 86400}'
end_production_s = '${fparse end_production * 86400}'
end_time_s = '${fparse end_time * 86400}'

# Mesh details
# borehole radus (m)
bh_r = 0.1
# model radius (m)
max_r = 1000
# aquifer base
min_y = 0
# aquifer top (m)
max_y = 30
# number of elements in radial direction
num_r = 100
# mesh bias in vertical direction in aquifer top
bias_r = 1.04
# number of elements in vertical direction
num_y = 30
# mesh bias in vertical direction in aquifer top
bias_y = 0.97
# depth (m) - for hydrostatic pressure initial condition
depth = 1000

# Aquifer properties
# Insitu temperature at centre of aquifer (degC)
temp = 50
# Vertical geothermal gradient (K/m).
geothermal_gradient = 0
# Gravity
gravity = -9.81
# Porosity
porosity = 0.25
# Permeability (m^2)
hor_perm = 1E-13
ver_perm = 1E-14
# Salinity (mass fraction)
xnacl = 0.1

# Input file
[Mesh]
    [aquifer]
        type = GeneratedMeshGenerator
        dim = 2
        nx = ${num_r}
        xmin = ${bh_r}
        xmax = ${max_r}
        bias_x = ${bias_r}
        bias_y = ${bias_y}
        ny = ${num_y}
        ymin = ${min_y}
        ymax = ${max_y}
    []
    coord_type = RZ
[]

[GlobalParams]
    PorousFlowDictator = dictator
    gravity '0 ${gravity} 0'
[]

[UserObjects]
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'pgas zi'
        number_fluid_phases = 2
        number_fluid_components = 2
    []
    [pc]
        type = PorousFlowCapillaryPressureVG
        alpha = 1e-3
        m = 0.4
        sat_lr = 0.1
    []
    [fs]
        type = PorousFlowBrineH2
        brine_fp = brine
        h2_fp = h2tab
        capillary_pressure = pc
    []
[]

[Variables]
    [pgas]
    []
    [zi]
        initial_condition = 0
        scaling = 1e4
    []
[]

[ICs]
    [pgas]
        type = FunctionIC
        variable = pgas
        function = insitu_pressure
    []
    
