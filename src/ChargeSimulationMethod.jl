module ChargeSimulationMethod
using StaticArrays

const ϵ = 8.8541878188e-12
ln = log

abstract type Charge end

struct LineCharge <: Charge
    pos::SVector{2,AbstractFloat}
    q::Complex
end

abstract type Conductor end

struct SolidConductor <: Conductor
    pos::SVector{2,AbstractFloat}
    r::AbstractFloat
    v::Complex
end

struct Simulation
    rmul::AbstractFloat # Multiplier of the radius of charges
    ncharges::Int
end

function creatematrices(conds::Vector{Conductor}, sim::Simulation)
    # j = Array{AbstractFloat}(undef, sim.ncharges * length(vec))
    # ϕ = Array{AbstractFloat}(undef, sim.ncharges * length(vec))
    # i = Array{AbstractFloat}(undef, sim.ncharges * length(vec))
    j = []
    ϕ = []
    i = []

    angs = Array{AbstractFloat}(undef, sim.ncharges)
    for i in 1:sim.ncharges
        angs[i] = (360/sim.ncharges) * (i - 1)
    end

    for c in conds
        for i in 1:sim.ncharges
            push!(i, c.pos + SVector(c.r * cosd(angs[i]), c.r * sind(angs[i])))
            push!(j, c.pos + sim.rmul*SVector(c.r * cosd(angs[i]), c.r * sind(angs[i])))
            push!(ϕ, c.v)
        end
    end

    return (j, ϕ, i)

end



end
