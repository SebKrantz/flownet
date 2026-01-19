#' African Cities and International Ports
#'
#' @name africa_cities_ports
#' @description
#' A spatial dataset containing 453 major African cities (population > 100,000) and international ports.
#' Cities are deduplicated within 50-100km radii, with populations aggregated from nearby settlements.
#' Port cities include cargo flow data from the World Bank Global Ports dataset.
#'
#' @format
#' A Simple feature collection (sf object, also inheriting from data.table) with 453 POINT features and 12 fields:
#' \describe{
#'   \item{city_country}{Character. Unique city-country identifier (e.g., "Cairo - Egypt", "Lagos - Nigeria").}
#'   \item{city}{Character. City name.}
#'   \item{country}{Character. Country name.}
#'   \item{iso2}{Character. ISO 3166-1 alpha-2 country code.}
#'   \item{iso3}{Character. ISO 3166-1 alpha-3 country code.}
#'   \item{admin_name}{Character. Administrative region or province name.}
#'   \item{capital}{Character. Capital status: "" (none), "admin" (administrative), "minor", or "primary" (national capital).}
#'   \item{population}{Numeric. City population including nearby settlements within 30km.}
#'   \item{port_locode}{Character. UN/LOCODE port identifier (empty string for non-port cities).}
#'   \item{port_name}{Character. Official port name (empty string for non-port cities).}
#'   \item{port_status}{Character. Port status code (empty string for non-port cities).}
#'   \item{outflows}{Numeric. Outflows in TEU in Q1 of 2020 (NA for non-port cities). 51 cities have port outflow data.}
#'   \item{geometry}{POINT. Spatial geometry in WGS 84 (EPSG:4326) coordinate reference system.}
#' }
#'
#' @details
#' The dataset was constructed by:
#' \enumerate{
#'   \item Selecting cities with population > 50,000 from Simplemaps World Cities database
#'   \item Weighting by administrative importance (capital status)
#'   \item Deduplicating within 50-100km radii, keeping largest weighted city
#'   \item Aggregating populations from settlements within 30km
#'   \item Matching with World Bank international ports within 30km
#' }
#'
#' The bounding box spans from approximately 34S to 37N latitude and 17W to 49E longitude,
#' covering continental Africa.
#'
#' @usage
#' data(africa_cities_ports)
#'
#' @source
#' City data from Simplemaps World Cities Database (\url{https://simplemaps.com/data/world-cities}).
#' Port data from World Bank Global International Ports dataset (\url{https://datacatalog.worldbank.org/search/dataset/0038118}).
#'
#' Dataset constructed for: Krantz, S. (2024). Optimal Investments in Africa's Road Network.
#' Policy Research Working Paper 10893. World Bank. \doi{10.1596/1813-9450-10893}.
#' Replication materials: \url{https://github.com/SebKrantz/OptimalAfricanRoads}.
#'
#' @seealso \code{\link{africa_network}}, \code{\link{africa_trade}}, \link{flowr-package}
#'
#' @examples
#' library(sf)
#' data(africa_cities_ports)
#' head(africa_cities_ports)
#'
#' # View largest cities
#' largest <- africa_cities_ports[order(-africa_cities_ports$population), ]
#' largest[1:10, c("city", "country", "population")]
#'
#' # Filter port cities
#' ports <- africa_cities_ports[!is.na(africa_cities_ports$port_locode), ]
#' nrow(ports)  # 51 ports
#'
#' \dontrun{
#' plot(africa_cities_ports["population"])
#' }
#'
"africa_cities_ports"

#' Trans-African Road Transport Network
#'
#' @name africa_network
#' @description
#' A spatial dataset representing a discretized road transport network connecting major African cities
#' and ports. The network combines existing road infrastructure (2,344 edges) with proposed new links
#' (481 edges) identified through network efficiency analysis. Each edge contains distance, travel time,
#' border crossing costs, terrain characteristics, and road upgrade cost estimates.
#'
#' @format
#' A Simple feature collection (sf object) with 2,825 LINESTRING features and 28 fields:
#' \describe{
#'   \item{from}{Integer. Origin node index (1 to 1,377).}
#'   \item{to}{Integer. Destination node index (2 to 1,379).}
#'   \item{from_ctry}{Character. Origin country ISO3 code (49 countries).}
#'   \item{to_ctry}{Character. Destination country ISO3 code (49 countries).}
#'   \item{FX}{Numeric. Origin node longitude.}
#'   \item{FY}{Numeric. Origin node latitude.}
#'   \item{TX}{Numeric. Destination node longitude.}
#'   \item{TY}{Numeric. Destination node latitude.}
#'   \item{sp_distance}{Numeric. Spherical (great-circle) distance in meters.}
#'   \item{distance}{Numeric. Road distance in meters from OSRM routing.}
#'   \item{duration}{Numeric. Travel duration in minutes from OSRM routing (NA for proposed links).}
#'   \item{speed_kmh}{Numeric. Average speed in km/h (distance/duration) (NA for proposed links).}
#'   \item{passes}{Numeric. Number of optimal inter-city routes passing through this edge (NA for proposed links).}
#'   \item{gravity}{Numeric. Sum of population gravity weights from routes using this edge (NA for proposed links).}
#'   \item{gravity_rd}{Numeric. Sum of road-distance-weighted gravity from routes (NA for proposed links).}
#'   \item{border_dist}{Numeric. Additional distance for border crossing in meters (0 for domestic links).}
#'   \item{total_dist}{Numeric. Total distance including border crossing penalty in meters.}
#'   \item{border_time}{Numeric. Additional time for border crossing in minutes.}
#'   \item{total_time}{Numeric. Total travel time including border crossing in minutes.}
#'   \item{duration_100kmh}{Numeric. Hypothetical travel time at 100 km/h in minutes.}
#'   \item{total_time_100kmh}{Numeric. Hypothetical total time at 100 km/h including border penalties.}
#'   \item{rugg}{Numeric. Terrain ruggedness index along the edge.}
#'   \item{pop_wpop}{Numeric. Population within corridor (WorldPop data).}
#'   \item{pop_wpop_km2}{Numeric. Population density per km2 along corridor.}
#'   \item{cost_km}{Numeric. Estimated road construction/maintenance cost per km in USD.}
#'   \item{upgrade_cat}{Character. Road upgrade category: "Nothing", "Asphalt Mix Resurfacing", "Mixed Works", "Upgrade", or NA.}
#'   \item{ug_cost_km}{Numeric. Upgrade cost per km in USD.}
#'   \item{add}{Logical. TRUE for proposed new links, FALSE for existing road network edges.}
#'   \item{geometry}{LINESTRING. Spatial geometry in WGS 84 (EPSG:4326) coordinate reference system.}
#' }
#'
#' @details
#' The network was constructed through the following process:
#' \enumerate{
#'   \item Computing optimal routes between all city pairs within 2,000km using OSRM
#'   \item Filtering routes using network efficiency criteria (alpha = 45 degrees, EU-grade efficiency)
#'   \item Intersecting and aggregating overlapping route segments
#'   \item Contracting the network to reduce complexity while preserving connectivity
#'   \item Identifying proposed new links that would improve network route efficiency
#'   \item Adding border crossing costs based on country pairs
#'   \item Computing terrain, population, and road cost attributes
#' }
#'
#' The \code{gravity} and \code{gravity_rd} fields measure edge importance based on the population
#' gravity model: routes between larger, closer cities contribute more weight to edges they traverse.
#'
#' The bounding box spans continental Africa from approximately 34S to 37N latitude
#' and 17W to 49E longitude.
#'
#' @usage
#' data(africa_network)
#'
#' @source
#' Road network derived from OpenStreetMap via OSRM routing.
#' Border crossing data from World Bank estimates.
#' Terrain data from SRTM elevation models.
#' Population data from WorldPop.
#'
#' Dataset constructed for: Krantz, S. (2024). Optimal Investments in Africa's Road Network.
#' Policy Research Working Paper 10893. World Bank. \doi{10.1596/1813-9450-10893}.
#' Replication materials: \url{https://github.com/SebKrantz/OptimalAfricanRoads}.
#'
#' @seealso \code{\link{africa_cities_ports}}, \code{\link{africa_segments}},
#'   \code{\link{africa_trade}}, \link{flowr-package}
#'
#' @examples
#' library(sf)
#' data(africa_network)
#' head(africa_network)
#'
#' # Existing vs proposed links
#' table(africa_network$add)
#'
#' # Cross-border links
#' cross_border <- africa_network[africa_network$from_ctry != africa_network$to_ctry, ]
#' nrow(cross_border)
#'
#' # Upgrade categories
#' table(africa_network$upgrade_cat, useNA = "ifany")
#'
#' \dontrun{
#' # Plot by gravity
#' plot(africa_network["gravity_rd"])
#'
#' # Highlight proposed new links
#' plot(africa_network[africa_network$add, "geometry"], col = "red", add = TRUE)
#' }
#'
"africa_network"

#' Raw Network Segments for Trans-African Transport Network
#'
#' @name africa_segments
#' @description
#' A dataset containing 14,358 raw network segments representing intersected road routes
#' between African cities. Each segment is defined by start and end coordinates with
#' aggregate importance metrics. This dataset is provided to demonstrate how package
#' functions like \code{\link[=consolidate_graph]{consolidate_graph()}} and
#' \code{\link[=simplify_network]{simplify_network()}} can process messy segment data
#' into clean analytical networks like \code{\link{africa_network}}.
#'
#' @format
#' A data frame with 14,358 rows and 7 columns:
#' \describe{
#'   \item{FX}{Numeric. Start point longitude (range: -17.4 to 49.2).}
#'   \item{FY}{Numeric. Start point latitude (range: -34.2 to 37.2).}
#'   \item{TX}{Numeric. End point longitude (range: -17.0 to 49.1).}
#'   \item{TY}{Numeric. End point latitude (range: -34.2 to 37.2).}
#'   \item{passes}{Integer. Number of optimal inter-city routes passing through this segment.
#'     Range: 1 to 1,615, median: 46.}
#'   \item{gravity}{Numeric. Sum of population gravity weights from routes using this segment.
#'     Computed as sum of (pop_origin * pop_destination / spherical_distance_km) / 1e9.}
#'   \item{gravity_rd}{Numeric. Sum of road-distance-weighted gravity from routes.
#'     Computed as sum of (pop_origin * pop_destination / road_distance_m) / 1e9.}
#' }
#'
#' @details
#' This dataset represents an intermediate stage in network construction, after routes have been
#' intersected but before network simplification. The segments have been simplified using
#' \code{\link[=linestrings_from_graph]{linestrings_from_graph()}} to retain only start and end coordinates.
#'
#' The segments can be used to demonstrate the flowr network processing workflow:
#' \enumerate{
#'   \item Convert segments to an sf LINESTRING object using \code{\link[=linestrings_from_graph]{linestrings_from_graph()}}
#'   \item Apply \code{\link[=consolidate_graph]{consolidate_graph()}} to merge nearby nodes
#'   \item Apply \code{\link[=simplify_network]{simplify_network()}} to remove intermediate nodes
#' }
#'
#' The \code{passes} field indicates how many optimal city-to-city routes use each segment,
#' serving as a measure of segment importance in the network. Higher values indicate
#' segments that are critical for efficient inter-city connectivity.
#'
#' @usage
#' data(africa_segments)
#'
#' @source
#' Derived from OpenStreetMap routing data via OSRM, processed through route intersection
#' and aggregation.
#'
#' Dataset constructed for: Krantz, S. (2024). Optimal Investments in Africa's Road Network.
#' Policy Research Working Paper 10893. World Bank. \doi{10.1596/1813-9450-10893}.
#' Replication materials: \url{https://github.com/SebKrantz/OptimalAfricanRoads}.
#'
#' @seealso \code{\link{africa_network}}, \code{\link[=consolidate_graph]{consolidate_graph()}},
#'   \code{\link[=simplify_network]{simplify_network()}}, \code{\link[=linestrings_from_graph]{linestrings_from_graph()}},
#'   \link{flowr-package}
#'
#' @examples
#' data(africa_segments)
#' head(africa_segments)
#'
#' # Summary statistics
#' summary(africa_segments[, c("passes", "gravity", "gravity_rd")])
#'
#' # Segments with highest traffic
#' africa_segments[order(-africa_segments$passes), ][1:10, ]
#'
#' \dontrun{
#' # Convert to sf and plot
#' library(sf)
#' segments_sf <- linestrings_from_graph(africa_segments)
#' plot(segments_sf["passes"])
#' }
#'
"africa_segments"

#' Intra-African Trade Flows by HS Section
#'
#' @name africa_trade
#' @description
#' A dataset containing bilateral trade flows between 47 African countries, aggregated by
#' HS (Harmonized System) section. Values represent annual averages over 2012-2022 from
#' the CEPII BACI database (HS96 nomenclature).
#'
#' @format
#' A data.table with 27,721 rows and 8 columns:
#' \describe{
#'   \item{iso3_o}{Factor. Exporter (origin) country ISO 3166-1 alpha-3 code (47 countries).}
#'   \item{iso3_d}{Factor. Importer (destination) country ISO 3166-1 alpha-3 code (47 countries).}
#'   \item{section_code}{Integer. HS section code (1 to 21).}
#'   \item{section_name}{Factor. HS section description (21 categories, e.g.,
#'     "Live animals and animal products", "Mineral products", "Machinery and mechanical appliances...").}
#'   \item{hs2_codes}{Factor. Comma-separated HS 2-digit codes within the section
#'     (e.g., "84, 85" for machinery).}
#'   \item{value}{Numeric. Trade value in thousands of USD (current prices).}
#'   \item{value_kd}{Numeric. Trade value in thousands of constant 2015 USD.}
#'   \item{quantity}{Numeric. Trade quantity in metric tons.}
#' }
#'
#' @details
#' The dataset provides bilateral trade flows aggregated from HS 6-digit product codes
#' (via HS 2-digit) to 21 HS sections. Trade values and quantities are annual averages
#' computed over the 2012-2022 period.
#'
#' HS sections cover broad product categories:
#' \itemize{
#'   \item Sections 1-5: Animal and vegetable products
#'   \item Sections 6-7: Chemical and plastic products
#'   \item Sections 8-14: Raw materials and manufactured goods
#'   \item Sections 15-16: Base metals and machinery
#'   \item Sections 17-21: Transport, instruments, and miscellaneous
#' }
#'
#' Note: Some country pairs may have sparse trade relationships. Very small values
#' indicate limited trade below typical reporting thresholds.
#'
#' @usage
#' data(africa_trade)
#'
#' @source
#' CEPII BACI Database (HS96 nomenclature), Version 202401b, released 2024-04-09.
#' Available at \url{http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37}.
#'
#' Reference: Gaulier, G. and Zignago, S. (2010). BACI: International Trade Database
#' at the Product-Level. The 1994-2007 Version. CEPII Working Paper, N 2010-23.
#'
#' @seealso \code{\link{africa_cities_ports}}, \code{\link{africa_network}}, \link{flowr-package}
#'
#' @examples
#' data(africa_trade)
#' head(africa_trade)
#'
#' # Number of trading pairs
#' length(unique(paste(africa_trade$iso3_o, africa_trade$iso3_d)))
#'
#' # Total trade by section
#' aggregate(value ~ section_name, data = africa_trade, FUN = sum)
#'
#' # Largest bilateral flows
#' africa_trade[order(-africa_trade$value), ][1:10, ]
#'
#' # Trade between specific countries
#' subset(africa_trade, iso3_o == "ZAF" & iso3_d == "NGA")
#'
"africa_trade"
