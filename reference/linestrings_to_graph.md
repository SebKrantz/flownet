# Convert Linestring to Graph

Convert Linestring to Graph

## Usage

``` r
linestrings_to_graph(
  lines,
  digits = 6,
  keep.cols = is.atomic,
  compute.length = TRUE
)
```

## Arguments

- lines:

  An sf data frame of LINESTRING geometries.

- digits:

  Numeric rounding applied to coordinates (to ensure that matching
  points across different linestrings is not impaired by numeric
  precision issues). Set to `NA/Inf/FALSE` to disable.

- keep.cols:

  Character vector of column names to keep from the input data frame.

- compute.length:

  Applies
  [`st_length()`](https://r-spatial.github.io/sf/reference/geos_measures.html)
  to and saves it as an additional column named `".length"`.

## Value

A data.frame representing the graph with columns:

- `edge` - Edge identifier

- `from` - Starting node ID

- `FX` - Starting node X-coordinate (longitude)

- `FY` - Starting node Y-coordinate (latitude)

- `to` - Ending node ID

- `TX` - Ending node X-coordinate (longitude)

- `TY` - Ending node Y-coordinate (latitude)

## See also

[simplify_network](https://sebkrantz.github.io/flowr/reference/simplify_network.md)
[flowr-package](https://sebkrantz.github.io/flowr/reference/flowr-package.md)

## Examples

``` r
library(flowr)
library(sf)

# Load existing network edges (exclude proposed new links)
africa_net <- africa_network[!africa_network$add, ]

# Convert network LINESTRING geometries to graph
graph <- linestrings_to_graph(africa_net)
head(graph)
#>   edge from        FX       FY to        TX       TY      .length from to
#> 1    1    1 -17.44671 14.69281  2 -17.04453 14.70297 64809.24 [m]    1  2
#> 2    2    2 -17.04453 14.70297  3 -16.63849 15.11222 72572.59 [m]    2  3
#> 3    3    2 -17.04453 14.70297  4 -16.46332 14.75651 74979.30 [m]    2  4
#> 4    4    2 -17.04453 14.70297  5 -16.42054 14.34236 95828.61 [m]    2  5
#> 5    5    3 -16.63849 15.11222  4 -16.46332 14.75651 46890.83 [m]    3  4
#> 6    6    3 -16.63849 15.11222 10 -16.24391 15.63569 71956.96 [m]    3 10
#>   from_ctry to_ctry        FX       FY        TX       TY sp_distance distance
#> 1       SEN     SEN -17.44671 14.69281 -17.04453 14.70297    43272.22  65988.0
#> 2       SEN     SEN -17.04453 14.70297 -16.63849 15.11222    63043.01  87290.5
#> 3       SEN     SEN -17.04453 14.70297 -16.46332 14.75651    62786.56  89469.5
#> 4       SEN     SEN -17.04453 14.70297 -16.42054 14.34236    78226.09  99998.0
#> 5       SEN     SEN -16.63849 15.11222 -16.46332 14.75651    43802.37  61111.5
#> 6       SEN     SEN -16.63849 15.11222 -16.24391 15.63569    71956.96  73391.0
#>   duration speed_kmh   passes     gravity  gravity_rd border_dist total_dist
#> 1    46.15  85.79155 104.3333 215.3302926 157.5228265           0    65988.0
#> 2    72.25  72.49038  10.0000  20.1382178  14.7376257           0    87290.5
#> 3    61.10  87.85876   1.5000  11.6420372  10.5577017           0    89469.5
#> 4    72.90  82.30288  55.0000 122.2011739  88.0080152           0    99998.0
#> 5    72.90  50.29753   1.0000   0.5259515   0.3953356           0    61111.5
#> 6    61.20  71.95196  31.0000  23.0342955  16.8520354           0    73391.0
#>   border_time total_time duration_100kmh total_time_100kmh      rugg  pop_wpop
#> 1           0      46.15         39.5928           39.5928 25203.645 670689.38
#> 2           0      72.25         52.3743           52.3743 23656.895 191412.11
#> 3           0      61.10         53.6817           53.6817 25761.240 177027.86
#> 4           0      72.90         59.9988           59.9988 24773.338 137874.58
#> 5           0      72.90         36.6669           36.6669 11331.973  65316.04
#> 6           0      61.20         44.0346           44.0346  1745.889 102835.91
#>   pop_wpop_km2  cost_km             upgrade_cat ug_cost_km   add
#> 1    2262.9802 699437.0 Asphalt Mix Resurfacing   165533.4 FALSE
#> 2     450.2150 605210.4             Mixed Works   325804.9 FALSE
#> 3     415.0577 607229.6 Asphalt Mix Resurfacing   143711.0 FALSE
#> 4     262.0664 581289.6 Asphalt Mix Resurfacing   137571.9 FALSE
#> 5     216.0512 520634.9                 Upgrade   440804.2 FALSE
#> 6     212.9464 415457.6             Mixed Works   223654.7 FALSE

# Graph contains edge, from/to nodes, and coordinates
names(graph)
#>  [1] "edge"              "from"              "FX"               
#>  [4] "FY"                "to"                "TX"               
#>  [7] "TY"                ".length"           "from"             
#> [10] "to"                "from_ctry"         "to_ctry"          
#> [13] "FX"                "FY"                "TX"               
#> [16] "TY"                "sp_distance"       "distance"         
#> [19] "duration"          "speed_kmh"         "passes"           
#> [22] "gravity"           "gravity_rd"        "border_dist"      
#> [25] "total_dist"        "border_time"       "total_time"       
#> [28] "duration_100kmh"   "total_time_100kmh" "rugg"             
#> [31] "pop_wpop"          "pop_wpop_km2"      "cost_km"          
#> [34] "upgrade_cat"       "ug_cost_km"        "add"              
```
