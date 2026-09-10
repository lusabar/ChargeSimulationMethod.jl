##########
# Timed Stuff
function time_domain_rms_field(cont::SVector, charges::Vector{LineCharge};
    ω::Float64=2π*60,
    nperiods::Int=20,
    samples_per_period::Int=200)

    println("Yeah, we also got here")
    T = 2π/ω
    ts = range(0.0, nperiods*T, length=nperiods*samples_per_period + 1)[1:(end-1)]

    Ex = zeros(Float64, length(ts))
    Ey = zeros(Float64, length(ts))

    geo = [fcoeff(cont, ch) for ch in charges]

    for (k, ch) in enumerate(charges)
        amp = √2 * abs(ch.q)
        phase = angle(ch.q)
        @. Ex += amp * cos(ω*ts + phase) * geo[k][1]
        @. Ey += amp * cos(ω*ts + phase) * geo[k][2]
    end

    Ex_rms = √(sum(Ex .^ 2) / length(Ex))
    Ey_rms = √(sum(Ey .^ 2) / length(Ey))
    return √(Ex_rms^2 + Ey_rms^2)
end

function ground_field_t(charges::Vector{<:Charge};
    HEIGHT=1)
    println("Yep, we got here")
    xs = LinRange(-100, 100, 41)
    y = []

    for x in xs
        push!(y, time_domain_rms_field(SVector(x, HEIGHT), charges))
    end
    return xs, y
end

#


function fcoeff(cont::SVector, charge::LineCharge)
    image_pos = SVector(charge.pos.x, -charge.pos.y)
    F = (1 / (2π*ϵ*norm(charge.pos - cont))) * normalize(charge.pos - cont)
    Fi = (1 / (2π*ϵ*norm(image_pos - cont))) * normalize(image_pos - cont)
    return F - Fi
end
fcoeff(charge::LineCharge, cont::SVector) = fcoeff(cont::SVector, charge::LineCharge)

function efield(pt::SVector, charges::Vector{<:Charge})
    E = SVector(0, 0)
    for ch in charges
        E += ch.q * fcoeff(pt, ch)
    end
    return norm(E)
end

function ground_field(charges::Vector{<:Charge};
    HEIGHT=1)
    xs = LinRange(-100, 100, 41)
    y = []

    for x in xs
        push!(y, efield(SVector(x, HEIGHT), charges))
    end
    return xs, y
end


##########################################################################################
function field_contribution(λ::Complex, zj::Complex, z::Complex)
    zj_image = conj(zj)          # mirror across y = 0
    return (λ / (2π * ϵ)) * (1/(z - zj) - 1/(z - zj_image))
end

"""
Total field phasor at point z from all charges Q at locations zloc.
"""
function total_field(Q::Vector{<:Complex}, zloc::Vector{<:Complex}, z::Complex)
    E = zero(ComplexF64)
    for i in eachindex(Q)
        E += field_contribution(Q[i], zloc[i], z)
    end
    return E  # Ex - i*Ey
end

function field_along_ground(Q::Vector{<:Complex}, zloc::Vector{<:Complex};
    y0::Float64=1.0,
    xs=range(-100.0, 100.0, length=41))
    Ex = Vector{Float64}(undef, length(xs))
    Ey = Vector{Float64}(undef, length(xs))

    for (k, x) in enumerate(xs)
        z = complex(x, y0)
        E = total_field(Q, zloc, z)
        Ex[k] = real(E)
        Ey[k] = -imag(E)   # note the sign flip from Ex - i*Ey
    end

    y = []
    for i in 1:length(Ex)
        push!(y, √(Ex[i]^2 + Ey[i]^2))
    end

    return xs, y
end
##########################################################################################
