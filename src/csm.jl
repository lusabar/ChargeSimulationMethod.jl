abstract type Charge end
mutable struct LineCharge <: Charge
    pos::SVector{2,AbstractFloat}
    q::Complex
end

struct ContourPoint
    pos::SVector{2,AbstractFloat}
    v::Complex
end

"""
Distance between two points
"""
distance(a::SVector, b::SVector) = norm(a - b)

"""
Distance between point a and image of point b
"""
distanceimg(a::SVector, b::SVector) = norm(a - SVector(b.x, -b.y))

"""
Calculates potential coefficient between contour point cont and a charge
"""
function pcoeff(cont::SVector, charge::LineCharge)
    return (1/(2π*ϵ)) * ln(distanceimg(cont, charge.pos) / distance(cont, charge.pos))
end
pcoeff(charge::LineCharge, cont::SVector) = pcoeff(cont::SVector, charge::LineCharge)

function check_accuracy(check_pts::Vector{ContourPoint}, charges::Vector{<:Charge};
    verbose=false,
    tol=1e-5)

    errs_percent = []

    for pt in check_pts
        v = 0
        for ch in charges
            v += ch.q * pcoeff(pt.pos, ch)
        end



        err = abs(pt.v - v) / abs(pt.v)
        push!(errs_percent, err*100)

        if verbose
            println("Error at point $(pt.pos): $(100*err)%")
        end

        if err > tol
            println("ERROR EXCEEDED LIMIT!")
        end
    end

    return errs_percent

end

function solve_csm(pts::Vector{ContourPoint}, charges::Vector{<:Charge};
    check_pts::Vector{ContourPoint})
    ϕ = [pt.v for pt in pts]

    P = zeros(length(pts), length(charges))
    for i ∈ eachindex(pts), j ∈ eachindex(charges)
        P[i, j] = pcoeff(pts[i].pos, charges[j])
    end

    Q = (P^-1) * ϕ

    # println(check_pts)
    for i in eachindex(Q)
        charges[i].q = Q[i]
    end

    errs_percent = check_accuracy(check_pts, charges, verbose=false)


    return charges, errs_percent
end
