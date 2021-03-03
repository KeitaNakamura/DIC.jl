using DIC
using Test

@testset "DIC searching" begin
    image = DIC.testimage("buffalo")
    subset = image[100:300, 300:500]
    indices, C = coarse_search(subset, image)
    @test C ≈ 1
    @test indices == CartesianIndices((100:300, 300:500))
end
