module ChargeSimulationMethod
export pol
using StaticArrays
using LinearAlgebra

include("csm.jl")

pol(m, a) = m * cis(deg2rad(a))
polsvec(m, a) = SVector(m*cosd(a), m*sind(a))

abstract type Simulation end

struct SimulationGround <: Simulation
    rmul::AbstractFloat # Multiplier of the radius of charges
    ncharges::Int
    range
end

struct SimulationConductor <: Simulation
    rmul::AbstractFloat # Multiplier of the radius of charges
    ncharges::Int
    angs_range # 
    rads_mul # Radius of the E points
end

struct SimulationGeneral <: Simulation
    rmul::AbstractFloat # Multiplier of the radius of charges
    ncharges::Int
    first_ang::AbstractFloat # Angle of the first charge
    epoints::Vector{SVector} # Points in which to determine E 
end

const ϵ = 8.8541878188e-12
ln = log

abstract type Conductor end

struct SolidConductor <: Conductor
    pos::SVector{2,AbstractFloat}
    r::AbstractFloat
    v::Complex
end

include("efield.jl")


function creatematrices(conds::Vector{SolidConductor}, sim::Simulation)
    # j = Array{AbstractFloat}(undef, sim.ncharges * length(vec))
    # ϕ = Array{AbstractFloat}(undef, sim.ncharges * length(vec))
    # i = Array{AbstractFloat}(undef, sim.ncharges * length(vec))
    j::Array{SVector} = []
    ϕ::Array{Complex} = []
    i::Array{SVector} = []

    angs = Array{AbstractFloat}(undef, sim.ncharges)
    for i in 1:sim.ncharges
        angs[i] = (360/sim.ncharges) * (i - 1)
    end



    for c in conds
        for x in 1:sim.ncharges
            push!(i, c.pos + SVector(c.r * cosd(angs[x]), c.r * sind(angs[x])))
            push!(j, c.pos + sim.rmul*SVector(c.r * cosd(angs[x]), c.r * sind(angs[x])))
            push!(ϕ, c.v)
        end
    end

    j = identity(j)

    return (j, ϕ, i)

end

function calculate_coeff(p1::SVector, p2::SVector)
    x = p1.x
    y = p1.y
    xj = p2.x
    yj = p2.y
    Pj = (1 / (2π*ϵ)) * ln((√((y + yj)^2 + (x - xj)^2)) / (√((y - yj)^2 + (x - xj)^2)))

    return Pj
end
function calculate_coeff_matrix(i::Vector{SVector}, j::Vector{SVector})
    M = zeros(length(i), length(j))
    for x in 1:length(i)
        for y in 1:length(j)
            M[x, y] = calculate_coeff(i[x], j[y])
        end
    end
    return M
end


################################################################################
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
################################################################################

# Distance considering actual charge and image charge
function distance_with_image(j::SVector, p::SVector)
    image = SVector(j.x, -j.y)
    #image.y = -image.y
    return (1/norm(j - p)) - (1/norm(image - p))
end


function check_accuracy(conds::Vector{SolidConductor}, sim::Simulation, j::Vector{SVector}, Q::Vector{Complex})
    ϕ_acc_target::Vector{Complex} = []
    checkpoints::Vector{SVector} = []

    angs = 0:10:350

    for c in conds
        for x in eachindex(angs)
            push!(checkpoints, c.pos + SVector(c.r * cosd(angs[x]), c.r * sind(angs[x])))
            push!(ϕ_acc_target, c.v)
        end
    end

    P = calculate_coeff_matrix(checkpoints, j)
    ϕ_actual = P * Q
    errs_percent = abs.(ϕ_acc_target - ϕ_actual) ./ abs.(ϕ_acc_target)

    return reshape(errs_percent, length(errs_percent), 1)
end

function calc_ground(conds::Vector{<:Conductor}, sim::Simulation)
    (j, ϕ, i) = creatematrices(conds, sim)
    P = calculate_coeff_matrix(i, j)
    #Q = ϕ \ P
    Q = (P^-1) * ϕ


    errs_percent = check_accuracy(conds, sim, j, Q)

    jz = [complex(a.x, a.y) for a in j]
    (x, y) = field_along_ground(Q, jz)



    return (x, y, errs_percent)
end

function create_charges_and_contour_points(conds::Vector{<:SolidConductor}, sim::Simulation)
    charges::Vector{LineCharge} = []
    contour_points::Vector{ContourPoint} = []

    if :first_ang in fieldnames(typeof(sim))
        first_ang = sim.first_ang
    else 
        first_ang = 0
    end

    angs = Array{AbstractFloat}(undef, sim.ncharges)
    for i in 1:sim.ncharges
        angs[i] = first_ang + (360/sim.ncharges) * (i - 1)
    end

    for c in conds, ang in angs
        push!(charges, LineCharge(c.pos + sim.rmul*SVector(c.r*cosd(ang), c.r*sind(ang)), 0))
        push!(contour_points, ContourPoint(c.pos + SVector(c.r*cosd(ang), c.r*sind(ang)), c.v))
        #push!(charges, c.pos + sim.rmul*SVector(c.r * cosd(angs[x]), c.r * sind(angs[x])))
    end

    return charges, contour_points
end

function create_checkpoints(conds::Vector{<:SolidConductor})
    check_points::Vector{ContourPoint} = []

    angs = 0:10:350

    for c in conds, ang in angs
        push!(check_points, ContourPoint(c.pos + SVector(c.r*cosd(ang), c.r*sind(ang)), c.v))
    end

    return check_points

end

function calc_ground2(conds::Vector{<:Conductor}, sim::Simulation)
    charges, contour_points = create_charges_and_contour_points(conds, sim)
    check_points = create_checkpoints(conds)

    charges, errs_percent = solve_csm(contour_points, charges, check_pts=check_points)

    ## Remover dps
    # println(charges)

    # jz = [complex(a.pos.x, a.pos.y) for a in charges]
    # Q = [ch.q for ch in charges]
    # (x, y) = field_along_ground(Q, jz)
    (x, y) = ground_field(charges, sim)

    return (x, y, errs_percent)
end

function calc_ground_t(conds::Vector{<:Conductor}, sim::Simulation)
    charges, contour_points = create_charges_and_contour_points(conds, sim)
    check_points = create_checkpoints(conds)

    charges = solve_csm(contour_points, charges, check_pts=check_points)

    # jz = [complex(a.pos.x, a.pos.y) for a in charges]
    # Q = [ch.q for ch in charges]
    # (x, y) = field_along_ground(Q, jz)
    (x, y) = ground_field_t(charges)

    return (x, y)
end

function calc_cond(conds::Vector{<:Conductor}, sim::SimulationConductor)
    charges, contour_points = create_charges_and_contour_points(conds, sim)
    check_points = create_checkpoints(conds)

    charges, errs_percent = solve_csm(contour_points, charges, check_pts=check_points)
    res = cond_field(conds, charges, sim)
    return res, errs_percent
end

function calc_general(conds::Vector{<:Conductor}, sim::SimulationGeneral)
    charges, contour_points = create_charges_and_contour_points(conds, sim)
    check_points = create_checkpoints(conds)

    charges, errs_percent = solve_csm(contour_points, charges, check_pts=check_points)
    res = general_field(conds, charges, sim)
    return res, charges, errs_percent
end





end
