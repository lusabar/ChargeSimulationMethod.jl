function fcoeff(cont::SVector, charge::LineCharge)
    image_pos = SVector(charge.pos.x, -charge.pos.y)
    F = (1 / (2π*ϵ*norm(charge.pos - cont))) * normalize(cont - charge.pos)
    Fi = (1 / (2π*ϵ*norm(image_pos - cont))) * normalize(cont - image_pos)
    return F - Fi
end
fcoeff(charge::LineCharge, cont::SVector) = fcoeff(cont::SVector, charge::LineCharge)

function efield(pt::SVector, charges::Vector{<:Charge})
    E = SVector(0, 0)
    for ch in charges
        E += ch.q * fcoeff(pt, ch)
    end

    return E # Returns the vector, the caller must use `norm` if desired
end

function ground_field(charges::Vector{<:Charge}, sim::SimulationGround;
    HEIGHT=1)
    xs = sim.range
    y = []

    for x in xs
        push!(y, efield(SVector(x, HEIGHT), charges) |> norm)
    end
    return xs, y
end

# Points where E will be determined
function create_epoints(conds::Vector{SolidConductor}, sim::SimulationConductor)
    epoints::Vector{SVector} = []
    for cond in conds, a in sim.angs_range, mul in sim.rads_mul 
        push!(epoints, cond.pos + polsvec(cond.r*mul, a))
    end
    return epoints
end

function cond_field(conds::Vector{<:Conductor}, charges::Vector{<:Charge}, sim::SimulationConductor)
    epoints = create_epoints(conds, sim)
    res = []

    for pt in epoints
        push!(res, [pt... norm(efield(pt, charges))])
    end

    return res
end

function general_field(conds::Vector{<:Conductor}, charges::Vector{<:Charge}, sim::SimulationGeneral)
    res = Tuple{SVector{2,Float64}, SVector{2,Float64}}[]

    for pt in sim.epoints
        push!(res, (pt, efield(pt, charges)))
        #push!(res, [pt efield(pt, charges)])
        #push!(res, [efield(pt, charges)])
    end

    return res
end

