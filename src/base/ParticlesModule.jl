module ParticlesModule

using LinearAlgebra
using Random
using ElasticArrays

#> "particle" type holds these things
mutable struct Particle
#> Phase-space position
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    vx_old::Float64
    vy_old::Float64
    vz_old::Float64
#> Diagnostics    
    energy::Float64
    momentum::Float64
    f_over_g::Float64
#> Electric field
    Ex::Float64
    Ey::Float64
    Ez::Float64
#> Magnetic field
    Bx::Float64
    By::Float64
    Bz::Float64
#> Time step for actions
    dt::Float64
end

#> "particle array" type to hold "particles" and the species parameters
mutable struct Particle_Array
    number_particles::Int
    mass::Float64
    charge::Float64
    total_kinetic_energy::Float64
    particle_shape_x::Float64
    particle_shape_y::Float64
    particle_shape_z::Float64
    particle::Array{Particle,1}
end

#> function to increase the size of the particle array size
function IncreaseParticleArraySize!(
        particle_array_struct,integer_increase
    )
    new_particle_size = length(particle_array_struct.particle)+integer_increase
    new_particle_array = Array{Particle,1}(undef,new_particle_size)
    for ip in 1:length(particle_array_struct.particle)
       new_particle_array[ip]=particle_array_struct.particle[ip]
    end 
    for ip in length(particle_array_struct.particle):new_particle_size
       new_particle_array[ip]=Particle()
    end 
    particle_array_struct.particle=new_particle_array
end

# outer constructor that zeroes fields
function Particle() 
    particle = Particle(
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0
    )
    return particle
end

# outer constructor that initializes RANDOM positions with RANDOM thermal velocities
function ParticleRandomInit(x_length,y_length,z_length,vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,time_step)
    particle = Particle(
        x_length*(rand()),
        y_length*(rand()),
        z_length*(rand()),
        vx_thermal_speed*randn(),
        vy_thermal_speed*randn(),
        vz_thermal_speed*randn(),
        vx_thermal_speed*randn(),
        vy_thermal_speed*randn(),
        vz_thermal_speed*randn(),
        0.0,0.0,1.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0,
        time_step
    )
    return particle
end

# outer constructor that initializes SPECIFIED positions with RANDOM thermal velocities
function ParticleOrderedInit(x,y,z,vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,time_step)
    particle = Particle(
        x,y,z,
        vx_thermal_speed*randn(),
        vy_thermal_speed*randn(),
        vz_thermal_speed*randn(),
        vx_thermal_speed*randn(),
        vy_thermal_speed*randn(),
        vz_thermal_speed*randn(),
        0.0,0.0,1.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0,
        time_step
    )
    return particle
end

#
function ParticleRandomInit(number_particles,x_length,y_length,z_length,vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,time_step)
    particles=Array{Particle,1}(undef,number_particles)
    for ip in 1:number_particles
       particles[ip]=ParticleRandomInit(
            x_length,y_length,z_length,
            vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
            time_step
        )
    end
    return particles
end

function ParticleOrderedInit(number_particles,x_length,y_length,z_length,vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,time_step)
    particles=Array{Particle,1}(undef,number_particles)
    for ip in 1:number_particles
       particles[ip]=ParticleOrderedInit(
            float(ip-1)*x_length/float(number_particles),0.0,0.0,
            vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
            time_step
        )
    end
    return particles
end

function ParticleVelocityOrderedInit_1D(number_particles,vx_thermal_speed)
    v_lim = 3.0*vx_thermal_speed
    dv_this_grid = 2.0*(v_lim)/(number_particles)
    #> Get a v-array
    v_this_grid = zeros(number_particles+1)
    v_this_grid[1] = -v_lim
    for i in 2:number_particles+1
        v_this_grid[i] = v_this_grid[i-1]+dv_this_grid
    end
    #> Weight based on Maxwellian
    w_this_grid = zeros(number_particles)
    vx_this_grid = zeros(number_particles)
    for i in 1:number_particles
        w_this_grid[i] = 0.5*dv_this_grid*(exp(-((v_this_grid[i]/vx_thermal_speed)^2))+exp(-((v_this_grid[i+1]/vx_thermal_speed)^2)))
        vx_this_grid[i] = 0.5*(v_this_grid[i]+v_this_grid[i+1])
    end
    w_this_grid[:] .*= number_particles/sum(w_this_grid[:])
    ind_this_grid = collect(Int,range(1,stop=length(w_this_grid),length=length(w_this_grid)))
    Random.shuffle!(ind_this_grid)
    tmp = zero(w_this_grid)
    for i in 1:length(w_this_grid)
        tmp[i] = w_this_grid[ind_this_grid[i]]
    end
    w_this_grid[:] = tmp[:]
    for i in 1:length(w_this_grid)
        tmp[i] = vx_this_grid[ind_this_grid[i]]
    end
    vx_this_grid[:] = tmp[:]
    return w_this_grid,vx_this_grid
end

function ParticleVelocityOrderedInit_1D_equal_weight(number_particles,vx_thermal_speed,nvx_cell=10)
    if(iseven(nvx_cell))
        nvx_cell += 1
    end
    v_lim = 3.0*vx_thermal_speed
    dv_this_grid = 2.0*(v_lim)/(nvx_cell)
    #> Get a v-array
    v_this_grid = zeros(nvx_cell+1)
    v_this_grid[1] = -v_lim
    for i in 2:nvx_cell+1
        v_this_grid[i] = v_this_grid[i-1]+dv_this_grid
    end
    
    #> Weight based on Maxwellian
    w_this_grid = zeros(nvx_cell)
    vx_this_grid = zeros(nvx_cell)
    for i in 1:nvx_cell
        w_this_grid[i] = 0.5*dv_this_grid*(exp(-((v_this_grid[i]/vx_thermal_speed)^2))+exp(-((v_this_grid[i+1]/vx_thermal_speed)^2)))
        vx_this_grid[i] = 0.5*(v_this_grid[i]+v_this_grid[i+1])
    end
    
    n_this_grid = copy(w_this_grid)
    n_this_grid[:] ./= sum(n_this_grid[:])
    n_this_grid[:] .*= number_particles
    n_this_grid[:] = round.(Int,n_this_grid)
    for i in 1:nvx_cell
        n_this_grid[i] = max(n_this_grid[i],1)
    end
    i_mid = Int(((nvx_cell-1)/2)+1)
    n_mismatch = number_particles-sum(n_this_grid)
    #if(isodd(n_mismatch))
    #    n_mismatch_start = (n_mismatch-1)
    #end
    n_this_grid[i_mid] += number_particles-sum(n_this_grid)
    
    w_particles = zeros(number_particles)
    v_particles = zeros(number_particles)
    ind_cell = 1
    for ip in 1:number_particles
        if(ip>sum(n_this_grid[1:ind_cell]))
            ind_cell +=1
        end
        v_particles[ip]=vx_this_grid[ind_cell]
        w_particles[ip]=w_this_grid[ind_cell]/n_this_grid[ind_cell]
    end
    w_particles[:] ./= sum(w_particles)
    
    ind_particles = collect(Int,range(1,stop=length(w_particles),length=length(w_particles)))
    Random.shuffle!(ind_particles)
    tmp = zero(w_particles)
    for i in 1:length(w_particles)
        tmp[i] = w_particles[ind_particles[i]]
    end
    w_particles[:] = tmp[:]
    for i in 1:length(w_particles)
        tmp[i] = v_particles[ind_particles[i]]
    end
    v_particles[:] = tmp[:]
    
    return w_particles,v_particles
end

function ParticlePhasespaceOrderedInit(number_particles,x_length,y_length,z_length,vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,time_step,gather_length)
    particles=Array{Particle,1}(undef,number_particles)
    #> First, load the particles uniformly in space
    for ip in 1:number_particles
       particles[ip]=ParticleOrderedInit(
            float(ip-1)*x_length/float(number_particles),0.0,0.0,
            0.0,0.0,0.0,
            time_step
        )
    end
    #> Make a temporary cell array
    cell_x = collect(range(0,stop=x_length,length=gather_length))
    pop!(cell_x)
    cell_x .+= 0.5*(cell_x[2]-cell_x[1])
    #> indices for the particles
    particle_indices = zeros(Int,length(cell_x)+1)
    particle_indices[1] = 1
    cell_ind = 1
    particle_ind = 1
    for ip in 1:number_particles
        particle_ind = findmin(abs.(cell_x.-particles[ip].x))[2]
        if(particle_ind == cell_ind)
            continue
        else
            cell_ind += 1
            particle_indices[cell_ind]=ip-1
        end
    end
    particle_indices[end]=number_particles
    #println("pic/",length(particle_indices),"/",particle_indices)
    #> Return evenly velocity-spaced unequally weighted particles
    for ic in 1:length(particle_indices)-1
        p_weights,p_velocities = ParticleVelocityOrderedInit_1D(particle_indices[ic+1]-particle_indices[ic]+1,vx_thermal_speed)
        for ip in 1:length(p_weights)
            particles[ip+particle_indices[ic]-1].vx = p_velocities[ip]
            particles[ip+particle_indices[ic]-1].f_over_g = p_weights[ip] 
        end
    end
    return particles
end
    
function ParticlePhasespaceOrderedInit_v2(number_particles,x_length,y_length,z_length,vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,time_step,gather_length,nvx_grid)
    particles=Array{Particle,1}(undef,number_particles)
    #> First, load the particles uniformly in space
    for ip in 1:number_particles
       particles[ip]=ParticleOrderedInit(
            float(ip-1)*x_length/float(number_particles),0.0,0.0,
            0.0,0.0,0.0,
            time_step
        )
    end
    #> Make a temporary cell array
    cell_x = collect(range(0,stop=x_length,length=gather_length))
    pop!(cell_x)
    cell_x .+= 0.5*(cell_x[2]-cell_x[1])
    #> indices for the particles
    particle_indices = zeros(Int,length(cell_x)+1)
    particle_indices[1] = 1
    cell_ind = 1
    particle_ind = 1
    for ip in 1:number_particles
        particle_ind = findmin(abs.(cell_x.-particles[ip].x))[2]
        if(particle_ind == cell_ind)
            continue
        else
            cell_ind += 1
            particle_indices[cell_ind]=ip-1
        end
    end
    particle_indices[end]=number_particles
    #println("pic/",length(particle_indices),"/",particle_indices)
    #> Return evenly velocity-spaced unequally weighted particles
    for ic in 1:length(particle_indices)-1
        p_weights,p_velocities = ParticleVelocityOrderedInit_1D_equal_weight(particle_indices[ic+1]-particle_indices[ic]+1,vx_thermal_speed,nvx_grid)
        for ip in 1:length(p_weights)
            particles[ip+particle_indices[ic]-1].vx = p_velocities[ip]
            particles[ip+particle_indices[ic]-1].f_over_g = p_weights[ip] 
        end
    end
    return particles
end

#
function ParticleArrayRandomInit(;
        number_particles=1,mass=1.0,charge=1.0,
        x_length=1.0,y_length=1.0,z_length=1.0,
        vx_thermal_speed=0.0,vy_thermal_speed=0.0,vz_thermal_speed=0.0,
        particle_shape_x=0.0,particle_shape_y=0.0,particle_shape_z=0.0,
        time_step=0.0
    )
    particle=ParticleRandomInit(
        number_particles,
        x_length,y_length,z_length,
        vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
        time_step
    )
    particle_array=Particle_Array(
        number_particles,mass,charge,0.0,
        particle_shape_x,particle_shape_y,particle_shape_z,
        particle
    )
    return particle_array
end

function ParticleArrayOrderedInit(;
        number_particles=1,mass=1.0,charge=1.0,
        x_length=1.0,y_length=1.0,z_length=1.0,
        vx_thermal_speed=0.0,vy_thermal_speed=0.0,vz_thermal_speed=0.0,
        particle_shape_x=0.0,particle_shape_y=0.0,particle_shape_z=0.0,
        time_step=0.0
    )
    particle=ParticleOrderedInit(
        number_particles,
        x_length,y_length,z_length,
        vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
        time_step
    )
    particle_array=Particle_Array(
        number_particles,mass,charge,0.0,
        particle_shape_x,particle_shape_y,particle_shape_z,
        particle
    )
    return particle_array
end
    
function ParticleArrayPhasespaceOrderedInit(;
        number_particles=1,mass=1.0,charge=1.0,
        x_length=1.0,y_length=1.0,z_length=1.0,
        vx_thermal_speed=0.0,vy_thermal_speed=0.0,vz_thermal_speed=0.0,
        particle_shape_x=0.0,particle_shape_y=0.0,particle_shape_z=0.0,
        time_step=0.0,
        gather_length=10,nvx_grid=11
    )
    particle=ParticlePhasespaceOrderedInit_v2(
        number_particles,
        x_length,y_length,z_length,
        vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
        time_step,gather_length,nvx_grid
    )
    particle_array=Particle_Array(
        number_particles,mass,charge,0.0,
        particle_shape_x,particle_shape_y,particle_shape_z,
        particle
    )
    return particle_array
end

function AddParticleDrift_x!(particle_array::Particle_Array,v_x)
    for ip in 1:particle_array.number_particles
       particle_array.particle[ip].vx += v_x
    end
end

function AddParticleRandomPositionPerturbation_x!(particle_array::Particle_Array,perturbation_size)
    for ip in 1:particle_array.number_particles
       particle_array.particle[ip].x += perturbation_size*randn()
    end
end
    
function AddParticleRandomVelocityPerturbation_x!(particle_array::Particle_Array,perturbation_size)
    for ip in 1:particle_array.number_particles
       particle_array.particle[ip].vx += perturbation_size*randn()
    end
end

function AddParticleCosinePositionPerturbation_x!(particle_array::Particle_Array,perturbation_size,perturbation_k)
    for ip in 1:particle_array.number_particles
       particle_array.particle[ip].x += perturbation_size*cos(particle_array.particle[ip].x*perturbation_k)
    end
end

function AddParticleCosineVelocityPerturbation_x!(particle_array::Particle_Array,velocity,perturbation_size,perturbation_k)
    for ip in 1:particle_array.number_particles
       particle_array.particle[ip].vx += velocity*(1.0+(perturbation_size*(rand()-0.5)*2.0))*cos(particle_array.particle[ip].x*perturbation_k)
    end
end

#
function UpdateParticlePosition_single!(particle::Particle)
    particle.x += particle.vx*particle.dt
    particle.y += particle.vy*particle.dt
    particle.z += particle.vz*particle.dt
end

#
function UpdateParticleVelocity_Boris_single!(particle::Particle,particle_array::Particle_Array)
    particle.vx_old = particle.vx
    particle.vy_old = particle.vy
    particle.vz_old = particle.vz    
#>
    boris_coef = 0.5*particle.dt*particle_array.charge/particle_array.mass
    tVec = boris_coef*[particle.Bx,particle.By,particle.Bz]
    sVec = tVec*2.0/(1.0+(tVec[1]^2)+(tVec[2]^2)+(tVec[3]^2))
#>
    particle.vx += boris_coef*particle.Ex
    particle.vy += boris_coef*particle.Ey
    particle.vz += boris_coef*particle.Ez
#>
    vdum = [particle.vx,particle.vy,particle.vz]
    vPrimeVec = vdum + LinearAlgebra.cross(vdum,tVec)
    vdum += LinearAlgebra.cross(vPrimeVec,sVec)
#>
    particle.vx = vdum[1]+boris_coef*particle.Ex
    particle.vy = vdum[2]+boris_coef*particle.Ey
    particle.vz = vdum[3]+boris_coef*particle.Ez
end

function UpdateParticlePosition!(particle_array::Particle_Array)
    for ip in 1:particle_array.number_particles
        UpdateParticlePosition_single!(particle_array.particle[ip])
    end
end

function UpdateParticleVelocity!(particle_array::Particle_Array)
    for ip in 1:particle_array.number_particles
        UpdateParticleVelocity_Boris_single!(particle_array.particle[ip],particle_array)
    end
end

function UpdateParticleKineticEnergy_single!(particle::Particle,particle_array::Particle_Array)
    particle.energy = 0.5*particle_array.mass*(
                        (particle.vx*particle.vx_old)
                       +(particle.vy*particle.vy_old)
                       +(particle.vz*particle.vz_old)
    )
end

function UpdateParticleKineticEnergy!(particle_array::Particle_Array)
    for ip in 1:particle_array.number_particles
        UpdateParticleKineticEnergy_single!(particle_array.particle[ip],particle_array)
    end
    particle_array.total_kinetic_energy = GetParticleKineticEnergy_total(particle_array::Particle_Array)
end

function GetParticleKineticEnergy_average(particle_array::Particle_Array)
    output_quantity=0.0
    for ip in 1:particle_array.number_particles
        output_quantity += particle_array.particle[ip].energy*particle_array.particle[ip].f_over_g
    end
    output_quantity /= float(particle_array.number_particles)
    return output_quantity
end

function GetParticleKineticEnergy_total(particle_array::Particle_Array)
    output_quantity=0.0
    for ip in 1:particle_array.number_particles
        output_quantity += (particle_array.particle[ip].energy*particle_array.particle[ip].f_over_g)
    end
    return output_quantity
end
    
function InjectionParticles_fromTheLeft(
        grid,density;
        particles_per_cell=10,
        vx_thermal_speed=1.0,vy_thermal_speed=0.0,vz_thermal_speed=0.0,
        buffer_fraction = 0.01,
        time_step=1.0
    )
    x_length = (vx_thermal_speed*3.0*time_step)
    dx = grid.x[2]-grid.x[1]
    x_cells = Int(round(ceil(x_length/dx)))
    x_length = x_cells*dx
    
    #> The weight of each particle is determined by the
    #> number of physical particles divided by the 
    #> number of simulation particles
    N_physical = (density*x_length)
    N_marker = Int(round(particles_per_cell*x_cells))
    f_over_g = float(N_physical)/float(N_marker)
    
    #> These are the POSSIBLE particles which may be injected
    #> into the simulation domain.
    possible_particles=ParticleRandomInit(
        N_marker,
        x_length,maximum(grid.y),maximum(grid.z),
        vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
        time_step
    )
    ip_injection = ElasticArray{Int}(undef, 1, 0)
    for ip in 1:N_marker
        x_step=possible_particles[ip].vx*possible_particles[ip].dt
        possible_particles[ip].x+=x_step-x_length+(dx*buffer_fraction)
        possible_particles[ip].f_over_g=f_over_g
        if(possible_particles[ip].x>= 0.0)
            append!(ip_injection,[ip])
        end
    end
    injected_particles = Array{Particle,1}(undef,length(ip_injection))
    for ip in 1:length(ip_injection)
        injected_particles[ip]=possible_particles[ip_injection[ip]]
    end
    possible_particles=injected_particles
    return possible_particles
end
        
# # #> Use a half-maxwellian
# function InjectionParticles_fromTheLeft(
#         grid,density;
#         particles_per_cell=10,
#         vx_thermal_speed=1.0,vy_thermal_speed=0.0,vz_thermal_speed=0.0,
#         buffer_fraction = 0.01,
#         time_step=1.0
#     )
#     x_length = (vx_thermal_speed*time_step)
#     dx = grid.x[2]-grid.x[1]
#     #x_cells = (x_length/dx)
    
#     #> The weight of each particle is determined by the
#     #> number of physical particles divided by the 
#     #> number of simulation particles
#     #> NOTE: the factor of 1/2 is because this is initializing only
#     #> the positive half of the velocities
#     N_physical = (density*x_length*0.5)
#     #N_physical = (density*x_length)
#     N_marker = Int(round(particles_per_cell*dx))
#     f_over_g = float(N_physical)/float(N_marker)
    
#     #> These are the POSSIBLE particles which may be injected
#     #> into the simulation domain.
#     possible_particles=ParticleRandomInit(
#         N_marker,
#         x_length,maximum(grid.y),maximum(grid.z),
#         vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
#         time_step
#     )
#     for ip in 1:N_marker
#         possible_particles[ip].vx=abs(possible_particles[ip].vx)
#         x_step=(possible_particles[ip].vx*possible_particles[ip].dt)
#         possible_particles[ip].x=x_step*rand()
#         possible_particles[ip].f_over_g=f_over_g
#     end
#     return possible_particles
# end
        
function ParticleInjectionFromLeft_Init(
        grid,density;
        particles_per_cell=1,mass=1.0,charge=1.0,
        vx_thermal_speed=0.0,vy_thermal_speed=0.0,vz_thermal_speed=0.0,
        particle_shape_x=0.0,particle_shape_y=0.0,particle_shape_z=0.0,
        buffer_fraction = 0.01,
        time_step=0.0
    )
    
    particle=InjectionParticles_fromTheLeft(
        grid,density,
        particles_per_cell=particles_per_cell,
        vx_thermal_speed=vx_thermal_speed,
        vy_thermal_speed=vy_thermal_speed,
        vz_thermal_speed=vz_thermal_speed,
        buffer_fraction = buffer_fraction,
        time_step=time_step
    );
    if(length(particle)==0)
        particle=Array{Particle,1}(undef,1)
        particle[1]=ParticleRandomInit(
            buffer_fraction*(grid.x[2]-grid.x[1]),
            buffer_fraction*(grid.y[2]-grid.y[1]),
            buffer_fraction*(grid.z[2]-grid.z[1]),
            vx_thermal_speed,vy_thermal_speed,vz_thermal_speed,
            time_step
        )
    end
    particle_array=Particle_Array(
        length(particle),mass,charge,0.0,
        particle_shape_x,particle_shape_y,particle_shape_z,
        particle
    );
    return particle_array
end
    
function ParticleInjectionFromLeft_Continue!(
        particle_array_struct,grid,density;
        particles_per_cell=1,
        vx_thermal_speed=0.0,
        vy_thermal_speed=0.0,
        vz_thermal_speed=0.0,
        time_step=0.0,
        buffer_fraction = 0.01,
        buffer_multiplier=10
    )
    injected_particles=InjectionParticles_fromTheLeft(
        grid,density,
        particles_per_cell=particles_per_cell,
        vx_thermal_speed=vx_thermal_speed,
        vy_thermal_speed=vy_thermal_speed,
        vz_thermal_speed=vz_thermal_speed,
        buffer_fraction = buffer_fraction,
        time_step=time_step
    );      
    N_injected = length(injected_particles)
    N_current = particle_array_struct.number_particles
    N_current_array = length(particle_array_struct.particle)
    if((N_injected+N_current)>N_current_array)
        IncreaseParticleArraySize!(
            particle_array_struct,Int(round(N_injected*buffer_multiplier))
        )
    end
    for ip in N_current+1:N_current+N_injected
        particle_array_struct.particle[ip]=injected_particles[ip-N_current]
    end
    particle_array_struct.number_particles=N_current+N_injected
end
    
end