import ChargeSimulationMethod as csm
using StaticArrays
using LinearAlgebra
using Test

@testset "ChargeSimulationMethod.jl" begin

    @testset "Hayt 5.21" begin
        q1 = csm.LineCharge(SVector(0, 1), 30e-9)
        q2 = csm.LineCharge(SVector(0, 2), 30e-9)

        P = SVector(1, 2)
        V = q1.q * csm.pcoeff(P, q1) +
            q2.q * csm.pcoeff(P, q2)
        E = q1.q * csm.fcoeff(P, q1) +
            q2.q * csm.fcoeff(P, q2)
        @test V ≈ 1200 atol=3
        @test E ≈ SVector(723, -18.9) atol=3
    end

    @testset "Electric Field Analysis Problem 9.2" begin
        q1 = csm.LineCharge(SVector(0, 5), 1e-9)

        P = SVector(0, 2)
        V = q1.q * csm.pcoeff(P, q1)
        E = q1.q * csm.fcoeff(P, q1)
        @test V ≈ 15.23 atol=1e-3
        @test E ≈ SVector(0, -8.56) atol=1e-3
    end

end
