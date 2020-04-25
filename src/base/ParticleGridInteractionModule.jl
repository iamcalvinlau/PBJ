module ParticleGridInteractionModule

using FFTW
using ParticlesModule
using GridModule
using ElasticArrays

function ScatterFieldToParticle_1D_single(
        field_x::Array{Float64,1},
        particle_x::Float64,
        grid_x::Array{Float64,1}
    )
    ind = searchsortedlast(grid_x, particle_x)
    ind = min(ind,length(grid_x)-1)
    ind = max(ind,1)
    weight_upper = (particle_x-grid_x[ind])/(grid_x[ind+1]-grid_x[ind])
    weight_lower = 1.0-weight_upper
    output_quantity_tmp = (weight_lower*field_x[ind])+(weight_upper*field_x[ind+1])
    return output_quantity_tmp
end

function ScatterFieldToParticle_Ex!(
        electric_field_x::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    for ip in 1:particle_array.number_particles
        particle_array.particle[ip].Ex = ScatterFieldToParticle_1D_single(electric_field_x,particle_array.particle[ip].x,grid_array.x)
    end
end

function ScatterParticleToGrid_1D_single(
        particle_x::Float64,
        grid_x::Array{Float64,1}
    )
    output_quantity_tmp = zero(grid_x)
    ind = searchsortedlast(grid_x, particle_x)
    ind = min(ind,length(grid_x)-1)
    ind = max(ind,1)
    weight_upper = (particle_x-grid_x[ind])/(grid_x[ind+1]-grid_x[ind])
    weight_lower = 1.0-weight_upper
    output_quantity_tmp[ind] = weight_lower
    output_quantity_tmp[ind+1] = weight_upper
    return output_quantity_tmp
end

function ScatterParticleToGrid_2D_single(
        particle_x::Float64,
        grid_x::Array{Float64,1},
        particle_y::Float64,
        grid_y::Array{Float64,1},
    )
    output_quantity_tmp = zeros(Float64,length(grid_x),length(grid_y))
    
    ind_x = searchsortedlast(grid_x, particle_x)
    ind_x = min(ind_x,length(grid_x)-1)
    ind_x = max(ind_x,1)
    weight_upper_x = (particle_x-grid_x[ind_x])/(grid_x[ind_x+1]-grid_x[ind_x])
    weight_lower_x = 1.0-weight_upper_x
    
    ind_y = searchsortedlast(grid_y, particle_y)
    ind_y = min(ind_y,length(grid_y)-1)
    ind_y = max(ind_y,1)
    weight_upper_y = (particle_y-grid_y[ind_y])/(grid_y[ind_y+1]-grid_y[ind_y])
    weight_lower_y = 1.0-weight_upper_y
    
    output_quantity_tmp[ind_x,ind_y] = weight_lower_x*weight_lower_y
    output_quantity_tmp[ind_x+1,ind_y] = weight_upper_x*weight_lower_y
    output_quantity_tmp[ind_x+1,ind_y+1] = weight_upper_x*weight_upper_y
    output_quantity_tmp[ind_x,ind_y+1] = weight_lower_x*weight_upper_y
    return output_quantity_tmp
end

function ScatterParticleToGrid_x(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    output_quantity=zero(grid_array.x)
    for ip in 1:particle_array.number_particles
        output_quantity[:] += ScatterParticleToGrid_1D_single(
            particle_array.particle[ip].x,grid_array.x
        )*particle_array.particle[ip].f_over_g
    end
    return output_quantity
end

function return_one(x,y,z)
    return 1.0
end

function return_one(x)
    return 1.0
end

function ScatterParticleToGrid_x_speedlimited(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array;
        beta=return_one
    )
    output_quantity=zero(grid_array.x)
    for ip in 1:particle_array.number_particles
        output_quantity[:] += ScatterParticleToGrid_1D_single(
            particle_array.particle[ip].x,grid_array.x
        )*(
            particle_array.particle[ip].f_over_g
            *beta(
                particle_array.particle[ip].vx,
                particle_array.particle[ip].vy,
                particle_array.particle[ip].vz
            )
        )
    end
    return output_quantity
end

function ScatterParticleToGrid_x_vx(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array,
        vgrid_array::GridModule.Grid_Array
    )
    output_quantity = zeros(Float64,length(grid_array.x),length(vgrid_array.x))
    
    for ip in 1:particle_array.number_particles
        output_quantity[:,:] += ScatterParticleToGrid_2D_single(
            particle_array.particle[ip].x,grid_array.x,
            particle_array.particle[ip].vx,vgrid_array.x
        )*particle_array.particle[ip].f_over_g
    end
    output_quantity=transpose(output_quantity)
    return output_quantity
end

function ApplyFiniteParticleShape_x_periodic(
        field_x::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    field_fft=FFTW.rfft(field_x[1:end-1])
    
    k_fft=range(0,stop=length(field_fft)-1,length=length(field_fft))
    k_fft*=2.0*pi/(grid_array.x[end]-grid_array.x[1])
    
    field_fft_out = field_fft.*exp.(-(k_fft*particle_array.particle_shape_x).^2)
    
    field_out=FFTW.irfft(field_fft_out,length(field_x[1:end-1]))
    push!(field_out,field_out[1])

    return field_out
end

function ApplyFiniteParticleShape_x(
        field_x::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    field_fft=FFTW.rfft(field_x[1:end])
    
    k_fft=range(0,stop=length(field_fft),length=length(field_fft))
    k_fft*=2.0*pi/(grid_array.x[end]-grid_array.x[1])
    
    field_fft_out = field_fft.*exp.(-(k_fft*particle_array.particle_shape_x).^2)
    
    field_out=FFTW.irfft(field_fft_out,length(field_x[1:end]))
    return field_out
end

function ApplyFiniteParticleShape_x(
        field_x::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array,
        left_boundary_value::Float64,
        right_boundary_value::Float64
    )
    
    n_x = length(field_x)
    n_zeros = n_x
    n_ones = n_x
    field_tmp = fill(0.0,n_x+n_zeros+n_ones)
    for i in 1:n_ones
        field_tmp[i]=left_boundary_value
    end
    for i in n_ones+1:n_ones+n_x
        field_tmp[i]=field_x[i-(n_ones)]
    end
    for i in n_ones+n_x+1:n_ones+n_x+n_zeros
        field_tmp[i]=right_boundary_value
    end
    
    field_fft=FFTW.rfft(field_tmp)

    k_fft=range(0,stop=length(field_fft),length=length(field_fft))
    k_fft*=2.0*pi/((grid_array.x[end]-grid_array.x[1])*3.0)
    
    field_fft_out = field_fft.*exp.(-(k_fft*particle_array.particle_shape_x).^2)
    
    field_out=FFTW.irfft(field_fft_out,length(field_tmp))
    return field_out[n_ones+1:n_ones+n_x]
end

function ApplyFiniteParticleShape_x_121(
        n_smooth::Int,
        field_in::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array,
        left_boundary_value::Float64,
        right_boundary_value::Float64
    )
    
    n_x = length(field_in)
    field_x = copy(field_in)
    field_tmp = fill(0.0,n_x)
    for i_smooth in 1:n_smooth
        field_tmp = fill(0.0,n_x)
        for i in 2:n_x-1
            field_tmp[i]=(
                (field_x[i-1]*0.25)
                +(field_x[i]*0.5)
                +(field_x[i+1]*0.25)
            )
        end
        field_tmp[1]=(
            (left_boundary_value*0.25)
            +(field_x[1]*0.5)
            +(field_x[2]*0.25)
        )
        field_tmp[end]=(
            (field_x[end-1]*0.25)
            +(field_x[end]*0.5)
            +(right_boundary_value*0.25)
        )
        field_x=copy(field_tmp)
    end
    field_out=field_tmp
    return field_out
end

function ApplyFiniteParticleShape_x_121(
        n_smooth::Int,
        field_in::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    
    n_x = length(field_in)
    field_x = copy(field_in)
    field_tmp = fill(0.0,n_x)
    for i_smooth in 1:n_smooth
        for i in 2:n_x-1
            field_tmp[i]=(
                (field_x[i-1]*0.25)
                +(field_x[i]*0.5)
                +(field_x[i+1]*0.25)
            )
        end
        field_tmp[1]=field_x[1]
        field_tmp[end]=field_x[end]
        field_x=copy(field_tmp)
    end
    field_out=field_tmp
    return field_out
end

function ApplyFiniteParticleShape_x_121_periodic(
        n_smooth::Int,
        field_in::Array{Float64,1},
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    #> Assumes a duplicate point
    n_x = length(field_in)
    field_x = copy(field_in)
    field_tmp = fill(0.0,n_x)
    for i_smooth in 1:n_smooth
        field_tmp = fill(0.0,n_x)
        for i in 2:n_x-1
            field_tmp[i]=(
                (field_x[i-1]*0.25)
                +(field_x[i]*0.5)
                +(field_x[i+1]*0.25)
            )
        end
        field_tmp[1]=(
            (field_x[end-1]*0.25)
            +(field_x[1]*0.5)
            +(field_x[2]*0.25)
        )
        field_tmp[end]=(
            (field_x[end-1]*0.25)
            +(field_x[end]*0.5)
            +(field_x[2]*0.25)
        )
        field_x=copy(field_tmp)
    end
    field_out=field_tmp
    return field_out
end

function ApplyPeriodicParticleBoundary!(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    for ip in 1:particle_array.number_particles
        particle_array.particle[ip].x = mod(particle_array.particle[ip].x-grid_array.x[1],grid_array.x[end]-grid_array.x[1])+grid_array.x[1]
        particle_array.particle[ip].y = mod(particle_array.particle[ip].y-grid_array.y[1],grid_array.y[end]-grid_array.y[1])+grid_array.y[1]
        particle_array.particle[ip].z = mod(particle_array.particle[ip].z-grid_array.z[1],grid_array.z[end]-grid_array.z[1])+grid_array.z[1]
    end
end

function ApplyAbsorptionParticleBoundary_Right!(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    ip_absorbed = ElasticArray{Int}(undef, 1, 0)
    for ip in 1:particle_array.number_particles
        if(particle_array.particle[ip].x>=maximum(grid_array.x))
            append!(ip_absorbed,[ip])
        end
    end
    N_pre_absorption = particle_array.number_particles
    N_post_absorption = N_pre_absorption-length(ip_absorbed)
    particle_array.number_particles=N_post_absorption
    for ip_a in 1:length(ip_absorbed)
        ip = ip_absorbed[ip_a]
        #> Shift to remove the absorbed particles
        particle_array.particle[ip:N_pre_absorption-1]=particle_array.particle[ip+1:N_pre_absorption]
        #> Shift down 1 index because of above shift
        ip_absorbed[ip_a:end]+=fill(-1,length(ip_absorbed)-ip_a+1)
    end
    #> Zero out the removed particles
    for ip in N_post_absorption+1:N_pre_absorption
        particle_array.particle[ip]=ParticlesModule.Particle()
    end
end

function ApplyAbsorptionParticleBoundary_LeftRight!(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array
    )
    ip_absorbed = ElasticArray{Int}(undef, 1, 0)
    for ip in 1:particle_array.number_particles
        if(particle_array.particle[ip].x>=maximum(grid_array.x))
            append!(ip_absorbed,[ip])
        elseif(particle_array.particle[ip].x<=minimum(grid_array.x))
            append!(ip_absorbed,[ip])
        end
    end
    N_pre_absorption = particle_array.number_particles
    N_post_absorption = N_pre_absorption-length(ip_absorbed)
    particle_array.number_particles=N_post_absorption
    for ip_a in 1:length(ip_absorbed)
        ip = ip_absorbed[ip_a]
        #> Shift to remove the absorbed particles
        particle_array.particle[ip:N_pre_absorption-1]=particle_array.particle[ip+1:N_pre_absorption]
        #> Shift down 1 index because of above shift
        ip_absorbed[ip_a:end]+=fill(-1,length(ip_absorbed)-ip_a+1)
    end
    #> Zero out the removed particles
    for ip in N_post_absorption+1:N_pre_absorption
        particle_array.particle[ip]=ParticlesModule.Particle()
    end
end

function ApplyRefluxParticleBoundary_Left!(
        particle_array::ParticlesModule.Particle_Array,
        grid_array::GridModule.Grid_Array,
        vx_thermal_speed::Float64,
        vy_thermal_speed::Float64,
        vz_thermal_speed::Float64
    )
    ip_reflux = ElasticArray{Int}(undef, 1, 0)
    for ip in 1:particle_array.number_particles
        if(particle_array.particle[ip].x<=minimum(grid_array.x))
            append!(ip_reflux,[ip])
        end
    end
    N_reflux=length(ip_reflux)
    for ip_a in 1:N_reflux
        ip = ip_reflux[ip_a]
        #> Re-thermalize the particles' velocities
        particle_array.particle[ip].vx=vx_thermal_speed*randn()
        particle_array.particle[ip].vy=vy_thermal_speed*randn()
        particle_array.particle[ip].vz=vz_thermal_speed*randn()
        #> Re-flect the particle back into the system
        particle_array.particle[ip].x=(2.0*minimum(grid_array.x))-particle_array.particle[ip].x
    end
end

end