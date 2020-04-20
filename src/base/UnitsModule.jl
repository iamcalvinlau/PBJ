module UnitsModule

#> "Unit" type holds these things
mutable struct Units
#> Set basic units 
    unit_density::Float64
    unit_velocity::Float64
    unit_charge::Float64
    unit_mass::Float64
#> Derive unit measurements from basic units
    unit_length::Float64
    unit_temperature::Float64
    unit_time::Float64
    unit_potential::Float64  
#> Normalized units in simulation
    norm_density::Float64
    norm_velocity::Float64
    norm_charge::Float64
    norm_mass::Float64
    norm_length::Float64
    norm_temperature::Float64
    norm_time::Float64
    norm_potential::Float64
    norm_permittivity::Float64
end

# function SetUnits(density)
#     #> Default set units
#     #>> electron mass
#     mass = 9.1093837015e-31
#     #>> electron charge
#     charge = 1.602176634e-19
#     #>> speed of light
#     speed_of_light_c = 299792458.0
#     #>> 1 eV = ? joules
#     eV_to_joules = 1.602e-19
#     #>> permittivity of free space
#     permittivity = 8.85418782e-12
    
#     #> Derived units
#     velocity = speed_of_light_c
#     time = sqrt(permittivity*mass/(charge^2)/(density))
#     length = velocity*time
#     potential = mass*(velocity^2)/charge
#     temperature = mass*(velocity^2)/eV_to_joules
#     #>> unit of permittivity
#     unit_permittivity = (charge^2)*(time^2)/(mass)/(length^3)
    
#     units = Units(
#         density,
#         velocity,charge,mass,
#         length,temperature,time,potential,
#         permittivity/unit_permittivity
#     )
#     return units
# end
function SetUnits(density,temperature)
    #> Input density assumed to be [m^-3]
    #> Input temperture assumed to be [eV]
    #> Default set units
    #>> electron mass (kg)
    mass = 9.1093837015e-31
    #>> electron charge (coulomb)
    charge = 1.602176634e-19
    #>> 1 eV to X joules
    eV_to_joules = 1.602e-19
    #>> electron thermal velocity
    thermal_velocity = sqrt((temperature*eV_to_joules)/mass)
    #>> permittivity of free space
    permittivity = 8.85418782e-12
    #>> plasma frequency
    omega_pe = sqrt((density*charge*charge)/(mass*permittivity))
    #>> debye length
    lambda_de = thermal_velocity/omega_pe
    
    #> Derived units
    velocity = thermal_velocity
    time = 1.0/omega_pe
    length = lambda_de
    potential = mass*(velocity^2)/charge
    #>> unit of permittivity
#    unit_permittivity = (charge^2)*(time^2)/(mass)/(length^3)
    unit_permittivity = (charge)/(potential*length)
    
    units = Units(
        density,
        velocity,charge,mass,
        length,temperature,time,potential,
        1.0,
        1.0,1.0,1.0,
        1.0,1.0,1.0,1.0,
        permittivity/unit_permittivity
    )
    return units
end

function GetThermalVelocity(temperature,mass)
    #>> 1 eV to X joules
    eV_to_joules = 1.602e-19
    #>> electron thermal velocity
    thermal_velocity = sqrt((temperature*eV_to_joules)/mass)
    return thermal_velocity
end

end