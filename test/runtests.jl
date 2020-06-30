using NetMSA
using Test

@testset "NetMSA.jl" begin
    # Write your tests here.
    S1 = "abcbcdem";
    S2 = "acbcfg";
    S3 = "abchimn";
    S4 = "abcbcjkm";

    L = [S1, S2, S3, S4];
    M = [
        ['a' 'a' 'a' 'a']
        ['b' 'c' 'b' 'b']
        ['c' 'b' 'c' 'c']
        ['b' 'c' 'h' 'b']
        ['c' 'f' 'i' 'c']
        ['d' 'g' 'm' 'j']
        ['e' missing 'n' 'k']
        ['m' missing missing 'm']
        ]

    @test isequal(NetMSA.createPeerMatrix(L), M);

    @test NetMSA.getposition('b', 2, M).indexes == [1, 3, 4]

    @test NetMSA.mostfrequent(M[1, :]) == (4, 'a');

    @test NetMSA.mostfrequent(M[2, :]) == (3, 'b');

    @test NetMSA.full(M[1, :]) == true;

    @test NetMSA.aligned(M[2, :]) == false;

    @test NetMSA.aligned(M[8, :]) == true;

    @test NetMSA.weight(M[1, :], w1=0.25, w2=0.5, w3=1.0) == 1.0;

    @test NetMSA.weight(M[2, :], w1=0.25, w2=0.5, w3=1.0) == 0.1875;

    @test NetMSA.weight(M[4, :], w1=0.25, w2=0.5, w3=1.0) == 0.125;

    @test NetMSA.weight(M[6, :], w1=0.25, w2=0.5, w3=1.0) == 0.0;

    @test NetMSA.weight(M[8, :], w1=0.25, w2=0.5, w3=1.0) == 0.25;

    @test sum(NetMSA.weight.(eachrow(M[2:end, :]))) == 0.875;

    @test NetMSA.objective(M, 2) == 2.625;

    @test_throws ArgumentError NetMSA.objective(M, 2, endindex=9)

end
