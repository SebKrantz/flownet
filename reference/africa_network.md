# Trans-African Road Transport Network

A spatial dataset representing a discretized road transport network
connecting major African cities and ports. The network combines existing
road infrastructure (2,344 edges) with proposed new links (481 edges)
identified through network efficiency analysis. Each edge contains
distance, travel time, border crossing costs, terrain characteristics,
and road upgrade cost estimates.

## Usage

``` r
data(africa_network)
```

## Format

A Simple feature collection (sf object) with 2,825 LINESTRING features
and 28 fields:

- from:

  Integer. Origin node index (1 to 1,377).

- to:

  Integer. Destination node index (2 to 1,379).

- from_ctry:

  Character. Origin country ISO3 code (49 countries).

- to_ctry:

  Character. Destination country ISO3 code (49 countries).

- FX:

  Numeric. Origin node longitude.

- FY:

  Numeric. Origin node latitude.

- TX:

  Numeric. Destination node longitude.

- TY:

  Numeric. Destination node latitude.

- sp_distance:

  Numeric. Spherical (great-circle) distance in meters.

- distance:

  Numeric. Road distance in meters from OSRM routing.

- duration:

  Numeric. Travel duration in minutes from OSRM routing (NA for proposed
  links).

- speed_kmh:

  Numeric. Average speed in km/h (distance/duration) (NA for proposed
  links).

- passes:

  Numeric. Number of optimal inter-city routes passing through this edge
  (NA for proposed links).

- gravity:

  Numeric. Sum of population gravity weights from routes using this edge
  (NA for proposed links).

- gravity_rd:

  Numeric. Sum of road-distance-weighted gravity from routes (NA for
  proposed links).

- border_dist:

  Numeric. Additional distance for border crossing in meters (0 for
  domestic links).

- total_dist:

  Numeric. Total distance including border crossing penalty in meters.

- border_time:

  Numeric. Additional time for border crossing in minutes.

- total_time:

  Numeric. Total travel time including border crossing in minutes.

- duration_100kmh:

  Numeric. Hypothetical travel time at 100 km/h in minutes.

- total_time_100kmh:

  Numeric. Hypothetical total time at 100 km/h including border
  penalties.

- rugg:

  Numeric. Terrain ruggedness index along the edge.

- pop_wpop:

  Numeric. Population within corridor (WorldPop data).

- pop_wpop_km2:

  Numeric. Population density per km2 along corridor.

- cost_km:

  Numeric. Estimated road construction/maintenance cost per km in USD.

- upgrade_cat:

  Character. Road upgrade category: "Nothing", "Asphalt Mix
  Resurfacing", "Mixed Works", "Upgrade", or NA.

- ug_cost_km:

  Numeric. Upgrade cost per km in USD.

- add:

  Logical. TRUE for proposed new links, FALSE for existing road network
  edges.

- geometry:

  LINESTRING. Spatial geometry in WGS 84 (EPSG:4326) coordinate
  reference system.

## Source

Road network derived from OpenStreetMap via OSRM routing. Border
crossing data from World Bank estimates. Terrain data from SRTM
elevation models. Population data from WorldPop.

Dataset constructed for: Krantz, S. (2024). Optimal Investments in
Africa's Road Network. Policy Research Working Paper 10893. World Bank.
[doi:10.1596/1813-9450-10893](https://doi.org/10.1596/1813-9450-10893) .
Replication materials:
<https://github.com/SebKrantz/OptimalAfricanRoads>.

## Details

The network was constructed through the following process:

1.  Computing optimal routes between all city pairs within 2,000km using
    OSRM

2.  Filtering routes using network efficiency criteria (alpha = 45
    degrees, EU-grade efficiency)

3.  Intersecting and aggregating overlapping route segments

4.  Contracting the network to reduce complexity while preserving
    connectivity

5.  Identifying proposed new links that would improve network route
    efficiency

6.  Adding border crossing costs based on country pairs

7.  Computing terrain, population, and road cost attributes

The `gravity` and `gravity_rd` fields measure edge importance based on
the population gravity model: routes between larger, closer cities
contribute more weight to edges they traverse.

The bounding box spans continental Africa from approximately 34S to 37N
latitude and 17W to 49E longitude.

## See also

[`africa_cities_ports`](https://sebkrantz.github.io/flownet/reference/africa_cities_ports.md),
[`africa_segments`](https://sebkrantz.github.io/flownet/reference/africa_segments.md),
[`africa_trade`](https://sebkrantz.github.io/flownet/reference/africa_trade.md),
[flownet-package](https://sebkrantz.github.io/flownet/reference/flownet-package.md)

## Examples

``` r
library(sf)
data(africa_network)
head(africa_network)
#> Simple feature collection with 6 features and 28 fields
#> Geometry type: LINESTRING
#> Dimension:     XY
#> Bounding box:  xmin: -17.44671 ymin: 14.34236 xmax: -16.24391 ymax: 15.63569
#> Geodetic CRS:  WGS 84
#>   from to from_ctry to_ctry        FX       FY        TX       TY sp_distance
#> 1    1  2       SEN     SEN -17.44671 14.69281 -17.04453 14.70297    43272.22
#> 2    2  3       SEN     SEN -17.04453 14.70297 -16.63849 15.11222    63043.01
#> 3    2  4       SEN     SEN -17.04453 14.70297 -16.46332 14.75651    62786.56
#> 4    2  5       SEN     SEN -17.04453 14.70297 -16.42054 14.34236    78226.09
#> 5    3  4       SEN     SEN -16.63849 15.11222 -16.46332 14.75651    43802.37
#> 6    3 10       SEN     SEN -16.63849 15.11222 -16.24391 15.63569    71956.96
#>   distance duration speed_kmh   passes     gravity  gravity_rd border_dist
#> 1  65988.0    46.15  85.79155 104.3333 215.3302926 157.5228265           0
#> 2  87290.5    72.25  72.49038  10.0000  20.1382178  14.7376257           0
#> 3  89469.5    61.10  87.85876   1.5000  11.6420372  10.5577017           0
#> 4  99998.0    72.90  82.30288  55.0000 122.2011739  88.0080152           0
#> 5  61111.5    72.90  50.29753   1.0000   0.5259515   0.3953356           0
#> 6  73391.0    61.20  71.95196  31.0000  23.0342955  16.8520354           0
#>   total_dist border_time total_time duration_100kmh total_time_100kmh      rugg
#> 1    65988.0           0      46.15         39.5928           39.5928 25203.645
#> 2    87290.5           0      72.25         52.3743           52.3743 23656.895
#> 3    89469.5           0      61.10         53.6817           53.6817 25761.240
#> 4    99998.0           0      72.90         59.9988           59.9988 24773.338
#> 5    61111.5           0      72.90         36.6669           36.6669 11331.973
#> 6    73391.0           0      61.20         44.0346           44.0346  1745.889
#>    pop_wpop pop_wpop_km2  cost_km             upgrade_cat ug_cost_km   add
#> 1 670689.38    2262.9802 699437.0 Asphalt Mix Resurfacing   165533.4 FALSE
#> 2 191412.11     450.2150 605210.4             Mixed Works   325804.9 FALSE
#> 3 177027.86     415.0577 607229.6 Asphalt Mix Resurfacing   143711.0 FALSE
#> 4 137874.58     262.0664 581289.6 Asphalt Mix Resurfacing   137571.9 FALSE
#> 5  65316.04     216.0512 520634.9                 Upgrade   440804.2 FALSE
#> 6 102835.91     212.9464 415457.6             Mixed Works   223654.7 FALSE
#>                         geometry
#> 1 LINESTRING (-17.44671 14.69...
#> 2 LINESTRING (-17.04453 14.70...
#> 3 LINESTRING (-17.04453 14.70...
#> 4 LINESTRING (-17.04453 14.70...
#> 5 LINESTRING (-16.63849 15.11...
#> 6 LINESTRING (-16.63849 15.11...

# Existing vs proposed links
table(africa_network$add)
#> 
#> FALSE  TRUE 
#>  2344   481 

# Cross-border links
cross_border <- africa_network[africa_network$from_ctry != africa_network$to_ctry, ]
nrow(cross_border)
#> [1] 457

# Upgrade categories
table(africa_network$upgrade_cat, useNA = "ifany")
#> 
#> Asphalt Mix Resurfacing             Mixed Works                 Nothing 
#>                     628                     901                       5 
#>                 Upgrade                    <NA> 
#>                     810                     481 

if (FALSE) { # \dontrun{
# Plot by gravity
plot(africa_network["gravity_rd"])

# Highlight proposed new links
plot(africa_network[africa_network$add, "geometry"], col = "red", add = TRUE)
} # }
```
