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

end
