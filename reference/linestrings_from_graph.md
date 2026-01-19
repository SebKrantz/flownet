# Convert Graph to Linestrings

Convert a graph data frame with node coordinates to an sf object with
LINESTRING geometries.

## Usage

``` r
linestrings_from_graph(graph_df, crs = 4326)
```

## Arguments

- graph_df:

  A data frame representing a graph with columns: `FX`, `FY`, `TX`, `TY`
  (starting and ending node coordinates), and optionally other columns
  to preserve.

- crs:

  Numeric or character (default: 4326). Coordinate reference system to
  assign to the output sf object.

## Value

An sf data frame with LINESTRING geometry, containing all columns from
`graph_df` except `FX`, `FY`, `TX`, and `TY`. Each row represents an
edge as a LINESTRING connecting the from-node (`FX`, `FY`) to the
to-node (`TX`, `TY`).

## Details

This function is the inverse operation of
[`linestrings_to_graph`](https://sebkrantz.github.io/flowr/reference/linestrings_to_graph.md).
It:

- Creates LINESTRING geometries from node coordinates (`FX`, `FY`, `TX`,
  `TY`)

- Removes the coordinate columns from the output

- Preserves all other columns from the input graph data frame

- Returns an sf object suitable for spatial operations and visualization

## See also

[linestrings_to_graph](https://sebkrantz.github.io/flowr/reference/linestrings_to_graph.md)
[flowr-package](https://sebkrantz.github.io/flowr/reference/flowr-package.md)

## Examples

``` r
library(flowr)
library(sf)

# Convert segments data frame to sf LINESTRING object
segments_sf <- linestrings_from_graph(africa_segments)
class(segments_sf)
#> [1] "sf"         "data.frame"
head(segments_sf)
#> Simple feature collection with 6 features and 3 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -17.44759 ymin: 14.40844 xmax: -16.63849 ymax: 21.28569
#> Geodetic CRS:  WGS 84
#>   passes  gravity gravity_rd                       geometry
#> 1     59 93.20669  67.218731 LINESTRING (-17.44759 14.69...
#> 2      6 59.98582  45.491608 LINESTRING (-17.44759 14.69...
#> 3     48 61.68941  42.121072 LINESTRING (-17.04861 14.70...
#> 4     11 31.51728  25.097660 LINESTRING (-17.04861 14.70...
#> 5     63  4.09779   2.981159 LINESTRING (-17.03345 20.93...
#> 6     10 20.13822  14.737626 LINESTRING (-16.96602 14.72...

if (FALSE) { # \dontrun{
# Plot segments colored by route importance
plot(segments_sf["passes"])
} # }
```
