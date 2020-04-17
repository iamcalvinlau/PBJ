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
    norm_permittivity::Float64
end

function SetUnits(density)
    #> Default set units
    #>> electron mass
    mass = 9.1093837015e-31
    #>> electron charge
    charge = 1.602176634e-19
    #>> speed of light
    speed_of_light_c = 299792458.0
    #>> 1 eV = ? joules
    eV_to_joules = 1.602e-19
    #>> permittivity of free space
    permittivity = 8.85418782e-12
    
    #> Derived units
    velocity = speed_of_light_c
    time = sqrt(permittivity*mass/(charge^2)/(density))
    length = velocity*time
    potential = mass*(velocity^2)/charge
    temperature = mass*(velocity^2)/eV_to_joules
    #>> unit of permittivity
    unit_permittivity = (charge^2)*(time^2)/(mass)/(length^3)
    
    units = Units(
        density,
        velocity,charge,mass,
        length,temperature,time,potential,
        permittivity/unit_permittivity
    )
    return units
end

end