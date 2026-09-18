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

    #println("Ex: $(E.x)")
    #println("Ey: $(E.y)")

    return norm(E)
    # The  original is above 
    #Ex, Ey = E[1], E[2]
    #Ax, Ay = abs(Ex), abs(Ey)
    #Δφ = angle(Ex) - angle(Ey)

    ## peak resultant magnitude of the elliptically-rotating field vector,
    ## not the RMS-combined sqrt(Ax^2+Ay^2)
    #disc = max(0.0, Ax^4 + Ay^4 + 2*Ax^2*Ay^2*cos(2Δφ))  # clamp tiny FP negatives
    #return sqrt((Ax^2 + Ay^2)/2 + 0.5*sqrt(disc))

end

function ground_field(charges::Vector{<:Charge}, sim::Simulation;
    HEIGHT=1)
    xs = sim.range
    y = []

    for x in xs
        push!(y, efield(SVector(x, HEIGHT), charges))
    end
    return xs, y
end

