import ChargeSimulationMethod as csm
using StaticArrays
using LinearAlgebra
using Test

@testset "ChargeSimulationMethod.jl" begin

    # Tests `pcoeff` and `fcoeff` functions for line charges
    #
    # The goal of this test is to check whether the functions for calculating
    # the potential coefficient `pcoeff` and electric field coefficient
    # `fcoeff` are working properly. Problem 5.21 from the 9th edition of
    # Engineering Electromagnetics by Hayt and Buck is used as a numerical test
    # case.
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


    # Tests `pcoeff` and `fcoeff` functions for line charges
    #
    # The goal of this test is to check whether the functions for calculating
    # the potential coefficient `pcoeff` and electric field coefficient
    # `fcoeff` are working properly. Problem 9.2 from the 1st edition of
    # Electric Field Analysis by Sivaji Chakravorti is used as a numerical test
    # case.
    @testset "Electric Field Analysis Problem 9.2" begin
        q1 = csm.LineCharge(SVector(0, 5), 1e-9)

        P = SVector(0, 2)
        V = q1.q * csm.pcoeff(P, q1)
        E = q1.q * csm.fcoeff(P, q1)
        @test V ≈ 15.23 atol=1e-3
        @test E ≈ SVector(0, -8.56) atol=1e-3
    end

    # Tests generation of line charges and determination of its values
    #
    # The goal of this test is to check whether the line charge is allocated at
    # the expected position and if its charge can be correctly determined. The
    # values from Problem 9.2 from the 1st edition of Electric Field Analysis
    # by Sivaji Chakravorti are used once more, but the problem is solved in
    # reverse order, that is, from V determine Q
    @testset "Electric Field Analysis Problem 9.2 Reversed" begin
        V = 15.23
        R = 5
        P = SVector(0, 2)
        conds = [csm.SolidConductor(SVector(0,7), R, V)]
        sim = csm.SimulationGeneral(2/5, 1, -90, [P])

        res, charges, _ = csm.calc_general(conds, sim)
        E = res[1][2]

        @test charges[1].q ≈ 1e-9 atol=1e-3
        @test E ≈ SVector(0, -8.56) atol=1e-3
    end

end
