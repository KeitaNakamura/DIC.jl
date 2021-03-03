```@meta
DocTestSetup = :(using DIC)
```

# DIC Searching

## Coarse search

Coarse search using zero-mean normalized cross-correlation can be performed by using `coarse_search` method.
First, you need to `load` your image file.

```@example 1
using DIC # hide
image = DIC.testimage("buffalo")
save("buffalo.tif", image) # hide
```

![](buffalo.tif)

Extracting a part of image can be done using sub-array.

```@example 1
subset = image[100:300, 300:500]
save("buffalo_part1.tif", subset) # hide
```

![](buffalo_part1.tif)

Searching this `subset` in `image` returns corresponding `CartesianIndices`.

```@example 1
indices, C = coarse_search(subset, image)
image[indices]
save("buffalo_searched.tif", ans) # hide
```

![](buffalo_searched.tif)

# Functions

```@index
Order = [:function]
Pages = ["DIC Searching.md"]
```

```@autodocs
Modules = [DIC]
Order   = [:function]
Pages   = ["searching.jl"]
```
