```@meta
DocTestSetup = :(using DIC)
```

# DIC Searching

## Coarse search

In coarse searching, the following *zero-mean normalized cross-correlation* is used.
```math
C = \frac{\sum_x \sum_y (A(x,y) - \bar{A}) (B(x,y) - \bar{B})}{\sqrt{\sum_x \sum_y (A(x,y) - \bar{A})^2 \sum_x \sum_y (B(x,y) - \bar{B})^2}}
```
With this correlation value, the subset having the highest value of ``C (0 \leq C \leq 1)`` will be searched.

First, you need to `load` your image file.

```@example 1
using DIC # hide
image = DIC.testimage("buffalo") # use `load(filename)` for your own image
save("buffalo.tif", image) # hide
```

![](buffalo.tif)

To extract a part of image, just create sub-array.

```@example 1
subset = image[100:300, 300:500]
save("buffalo_part1.tif", subset) # hide
```

![](buffalo_part1.tif)

Searching this `subset` in `image` returns corresponding `CartesianIndices`.

```@example 1
indices, C = coarse_search(subset, image)
image[indices] # this should be the same as `image[100:300, 300:500]`
save("buffalo_searched.tif", ans) # hide
```

![](buffalo_searched.tif)

`C` is the correlation value defined in the range from 0 to 1.
In above example, the result `C` should be `1` since the subset is exactly a part of the original image.
You should use a subset after a certain event in practice.

## Fine search

The coarse search only considers the rigid translation of the subset with 1 pixel resolution.
However, in fine search, searching is performed with sub-pixel resolution using linear interpolation, and the deformation of subset is also taken into account.
See [this study](https://link.springer.com/article/10.1007%2FBF02321405) for more detail.

```@example 2
using DIC # hide
image = DIC.testimage("buffalo")
subset = image[100:300, 300:500]
center, C = fine_search(subset, image, CartesianIndices((101:301, 301:501)))
```

Note that, in above example, using `coarse_search` is enough to get an exact solution.

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
