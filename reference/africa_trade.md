# Intra-African Trade Flows by HS Section

A dataset containing bilateral trade flows between 47 African countries,
aggregated by HS (Harmonized System) section. Values represent annual
averages over 2012-2022 from the CEPII BACI database (HS96
nomenclature).

## Usage

``` r
data(africa_trade)
```

## Format

A data.table with 27,721 rows and 8 columns:

- iso3_o:

  Factor. Exporter (origin) country ISO 3166-1 alpha-3 code (47
  countries).

- iso3_d:

  Factor. Importer (destination) country ISO 3166-1 alpha-3 code (47
  countries).

- section_code:

  Integer. HS section code (1 to 21).

- section_name:

  Factor. HS section description (21 categories, e.g., "Live animals and
  animal products", "Mineral products", "Machinery and mechanical
  appliances...").

- hs2_codes:

  Factor. Comma-separated HS 2-digit codes within the section (e.g.,
  "84, 85" for machinery).

- value:

  Numeric. Trade value in thousands of USD (current prices).

- value_kd:

  Numeric. Trade value in thousands of constant 2015 USD.

- quantity:

  Numeric. Trade quantity in metric tons.

## Source

CEPII BACI Database (HS96 nomenclature), Version 202401b, released
2024-04-09. Available at
<https://www.cepii.fr/DATA_DOWNLOAD/baci/doc/baci_webpage.html>.

Reference: Gaulier, G. and Zignago, S. (2010). BACI: International Trade
Database at the Product-Level. The 1994-2007 Version. CEPII Working
Paper, N 2010-23.

## Details

The dataset provides bilateral trade flows aggregated from HS 6-digit
product codes (via HS 2-digit) to 21 HS sections. Trade values and
quantities are annual averages computed over the 2012-2022 period.

HS sections cover broad product categories:

- Sections 1-5: Animal and vegetable products

- Sections 6-7: Chemical and plastic products

- Sections 8-14: Raw materials and manufactured goods

- Sections 15-16: Base metals and machinery

- Sections 17-21: Transport, instruments, and miscellaneous

Note: Some country pairs may have sparse trade relationships. Very small
values indicate limited trade below typical reporting thresholds.

## See also

[`africa_cities_ports`](https://sebkrantz.github.io/flownet/reference/africa_cities_ports.md),
[`africa_network`](https://sebkrantz.github.io/flownet/reference/africa_network.md),
[flownet-package](https://sebkrantz.github.io/flownet/reference/flownet-package.md)

## Examples

``` r
data(africa_trade)
head(africa_trade)
#>   iso3_o iso3_d section_code
#> 1    AGO    BDI           16
#> 2    AGO    BEN            1
#> 3    AGO    BEN            4
#> 4    AGO    BEN            5
#> 5    AGO    BEN            6
#> 6    AGO    BEN            7
#>                                                                                                                                                                                                  section_name
#> 1 Machinery and mechanical appliances, electrical equipment, parts thereof, sound recorders and reproducers, television image and souch recorders and reproducers, and parts and accessories of such articles
#> 2                                                                                                                                                                            Live animals and animal products
#> 3                                                                                                           Prepared foodstuffs, beverages, spirits and vinegar, tobacco and manufactured tobacco substitutes
#> 4                                                                                                                                                                                            Mineral products
#> 5                                                                                                                                                               Product of the chemicals or allied industries
#> 6                                                                                                                                                  Plastics and articles thereof, rubber and articles thereof
#>        hs2_codes     value     value_kd  quantity
#> 1             84    5.1040    5.1040000    0.9870
#> 2              3 6810.4528 6527.0068247 8777.0133
#> 3     19, 20, 22    9.5135    9.0984152   20.7370
#> 4             27 3583.7897 3536.4303649 8369.7665
#> 5 30, 32, 33, 34    0.9435    0.9319342    0.0275
#> 6         39, 40   71.2590   71.9540899   31.2510

# Number of trading pairs
length(unique(paste(africa_trade$iso3_o, africa_trade$iso3_d)))
#> [1] 1971

# Total trade by section
aggregate(value ~ section_name, data = africa_trade, FUN = sum)
#>                                                                                                                                                                                                    section_name
#> 1                                                                                                 Animal or vegetable fats and oils and their cleavage products; prepared edible fats; animal or vegetble waxes
#> 2                                                                                                                                                            Arms and ammunition, parts and accessories thereof
#> 3                                                                                               Articles of stone, plaster, cement, asbestos, mica, or similar materials, ceramic products, glass and glassware
#> 4                                                                                                                                                                        Base metals and articles of base metal
#> 5  Footwear, headgear, umbrellas, sun umbrellas, walking-sticks, seat-sticks, whips, riding-crops and parts thereof, prepared feathers and articles made therewith, artificial flowers, articles of human hair.
#> 6                                                                                                                                                                              Live animals and animal products
#> 7   Machinery and mechanical appliances, electrical equipment, parts thereof, sound recorders and reproducers, television image and souch recorders and reproducers, and parts and accessories of such articles
#> 8                                                                                                                                                                                              Mineral products
#> 9                                                                                                                                                                           Miscellaneous manufactured articles
#> 10                                              Natural or cultured pearls, precious or semi-precious stones, precious metals, metal clad with precious metal, and articles thereof, imitation jewellery, coins
#> 11                Optical, photographic, cinematographic, measuring, checking, precision, medical or surgical instruments and apparatus, clocks and watches, musical instruments, parts and accessories thereof
#> 12                                                                                                                                                   Plastics and articles thereof, rubber and articles thereof
#> 13                                                                                                            Prepared foodstuffs, beverages, spirits and vinegar, tobacco and manufactured tobacco substitutes
#> 14                                                                                                                                                                Product of the chemicals or allied industries
#> 15                                                             Pulp of wood or of other fibrous cellulosic material, recovered (waste and scrap) paper or paperboard, paper and paperboard and articles thereof
#> 16                          Raw hides and skins, leather, furskins and articles thereof, saddlery and harness, travel goods, handbags and similar containers, articles of animal gut (other than silk-worm gut)
#> 17                                                                                                                                                                                 Textile and textile articles
#> 18                                                                                                                                                                                           Vegetable products
#> 19                                                                                                                                               Vehicles, aircraft, vessels and associated transport equipment
#> 20                                            Wood and articles of wood, wood charcoal, cork and articles of cork, manufacturers of straw, of esparto or of other plaiting materials, basketwork and wickerwork
#> 21                                                                                                                                                                Works of art, collectors' pieces and antiques
#>          value
#> 1   1613945.96
#> 2    175980.52
#> 3   1323219.72
#> 4   9615324.31
#> 5    492184.56
#> 6   3014069.79
#> 7   7144013.10
#> 8  25996492.65
#> 9    728451.52
#> 10  7376528.86
#> 11   654594.12
#> 12  3019133.56
#> 13  7824476.75
#> 14  9326746.35
#> 15  1753353.32
#> 16   171264.71
#> 17  2415768.17
#> 18  4772233.71
#> 19 12884142.52
#> 20   744624.84
#> 21    63022.15

# Largest bilateral flows
africa_trade[order(-africa_trade$value), ][1:10, ]
#>       iso3_o iso3_d section_code
#> 19294    NGA    ZAF            5
#> 18678    NGA    CIV            5
#> 6977     DZA    TUN            5
#> 664      AGO    ZAF            5
#> 26129    ZAF    MOZ            5
#> 5477     COD    TZA           15
#> 10431    GHA    ZAF           14
#> 25654    ZAF    BWA            5
#> 16026    MOZ    ZAF            5
#> 18764    NGA    DZA            5
#>                                                                                                                                                          section_name
#> 19294                                                                                                                                                Mineral products
#> 18678                                                                                                                                                Mineral products
#> 6977                                                                                                                                                 Mineral products
#> 664                                                                                                                                                  Mineral products
#> 26129                                                                                                                                                Mineral products
#> 5477                                                                                                                           Base metals and articles of base metal
#> 10431 Natural or cultured pearls, precious or semi-precious stones, precious metals, metal clad with precious metal, and articles thereof, imitation jewellery, coins
#> 25654                                                                                                                                                Mineral products
#> 16026                                                                                                                                                Mineral products
#> 18764                                                                                                                                                Mineral products
#>                                hs2_codes     value  value_kd     quantity
#> 19294                         25, 26, 27 3128312.6 3041533.2  5567383.378
#> 18678                         25, 26, 27 1572964.5 1524679.8  2718137.515
#> 6977                              25, 27 1393337.9 1359928.5  1477835.669
#> 664                           25, 26, 27 1204953.6 1198880.0  2473175.036
#> 26129                         25, 26, 27 1150394.2 1064112.0  6965476.414
#> 5477  72, 73, 74, 76, 78, 79, 81, 82, 83 1105493.8 1026106.0   187368.177
#> 10431                                 71 1093580.9 1106132.6      569.086
#> 25654                         25, 26, 27 1014182.6  979699.1  1943577.323
#> 16026                         25, 26, 27  892348.6  861015.7 15235699.239
#> 18764                             26, 27  866681.5  907334.7     6418.847

# Trade between specific countries
subset(africa_trade, iso3_o == "ZAF" & iso3_d == "NGA")
#>       iso3_o iso3_d section_code
#> 26230    ZAF    NGA            1
#> 26231    ZAF    NGA            2
#> 26232    ZAF    NGA            3
#> 26233    ZAF    NGA            4
#> 26234    ZAF    NGA            5
#> 26235    ZAF    NGA            6
#> 26236    ZAF    NGA            7
#> 26237    ZAF    NGA            8
#> 26238    ZAF    NGA            9
#> 26239    ZAF    NGA           10
#> 26240    ZAF    NGA           11
#> 26241    ZAF    NGA           12
#> 26242    ZAF    NGA           13
#> 26243    ZAF    NGA           14
#> 26244    ZAF    NGA           15
#> 26245    ZAF    NGA           16
#> 26246    ZAF    NGA           17
#> 26247    ZAF    NGA           18
#> 26248    ZAF    NGA           19
#> 26249    ZAF    NGA           20
#> 26250    ZAF    NGA           21
#>                                                                                                                                                                                                       section_name
#> 26230                                                                                                                                                                             Live animals and animal products
#> 26231                                                                                                                                                                                           Vegetable products
#> 26232                                                                                                Animal or vegetable fats and oils and their cleavage products; prepared edible fats; animal or vegetble waxes
#> 26233                                                                                                            Prepared foodstuffs, beverages, spirits and vinegar, tobacco and manufactured tobacco substitutes
#> 26234                                                                                                                                                                                             Mineral products
#> 26235                                                                                                                                                                Product of the chemicals or allied industries
#> 26236                                                                                                                                                   Plastics and articles thereof, rubber and articles thereof
#> 26237                          Raw hides and skins, leather, furskins and articles thereof, saddlery and harness, travel goods, handbags and similar containers, articles of animal gut (other than silk-worm gut)
#> 26238                                            Wood and articles of wood, wood charcoal, cork and articles of cork, manufacturers of straw, of esparto or of other plaiting materials, basketwork and wickerwork
#> 26239                                                             Pulp of wood or of other fibrous cellulosic material, recovered (waste and scrap) paper or paperboard, paper and paperboard and articles thereof
#> 26240                                                                                                                                                                                 Textile and textile articles
#> 26241 Footwear, headgear, umbrellas, sun umbrellas, walking-sticks, seat-sticks, whips, riding-crops and parts thereof, prepared feathers and articles made therewith, artificial flowers, articles of human hair.
#> 26242                                                                                              Articles of stone, plaster, cement, asbestos, mica, or similar materials, ceramic products, glass and glassware
#> 26243                                              Natural or cultured pearls, precious or semi-precious stones, precious metals, metal clad with precious metal, and articles thereof, imitation jewellery, coins
#> 26244                                                                                                                                                                       Base metals and articles of base metal
#> 26245  Machinery and mechanical appliances, electrical equipment, parts thereof, sound recorders and reproducers, television image and souch recorders and reproducers, and parts and accessories of such articles
#> 26246                                                                                                                                               Vehicles, aircraft, vessels and associated transport equipment
#> 26247                Optical, photographic, cinematographic, measuring, checking, precision, medical or surgical instruments and apparatus, clocks and watches, musical instruments, parts and accessories thereof
#> 26248                                                                                                                                                           Arms and ammunition, parts and accessories thereof
#> 26249                                                                                                                                                                          Miscellaneous manufactured articles
#> 26250                                                                                                                                                                Works of art, collectors' pieces and antiques
#>                                                    hs2_codes       value
#> 26230                                          1, 2, 3, 4, 5   4047.6297
#> 26231                         6, 7, 8, 9, 10, 11, 12, 13, 14  44466.0138
#> 26232                                                     15    306.9128
#> 26233                     16, 17, 18, 19, 20, 21, 22, 23, 24  51346.5875
#> 26234                                             25, 26, 27  18430.5662
#> 26235             28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38  70410.8983
#> 26236                                                 39, 40 107293.9331
#> 26237                                             41, 42, 43    367.3535
#> 26238                                             44, 45, 46   1321.2536
#> 26239                                             47, 48, 49  22110.2727
#> 26240 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63   9011.9387
#> 26241                                         64, 65, 66, 67    973.3085
#> 26242                                             68, 69, 70   4164.7806
#> 26243                                                     71    175.0221
#> 26244             72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83  75819.1181
#> 26245                                                 84, 85  85302.7382
#> 26246                                         86, 87, 88, 89 143485.1386
#> 26247                                             90, 91, 92  11087.3447
#> 26248                                                     93    852.2147
#> 26249                                             94, 95, 96   8806.5570
#> 26250                                                     97    149.3303
#>          value_kd     quantity
#> 26230   3873.1014 1.352117e+03
#> 26231  42116.8754 5.019097e+04
#> 26232    279.5778 2.736282e+02
#> 26233  49500.2415 2.533764e+04
#> 26234  17881.1629 1.625416e+05
#> 26235  67480.4950 3.381519e+04
#> 26236 102448.6112 8.191542e+04
#> 26237    356.0907 4.876490e+02
#> 26238   1300.7074 1.766760e+03
#> 26239  21560.8662 1.317702e+04
#> 26240   8653.2046 4.482534e+03
#> 26241    958.3432 2.306605e+02
#> 26242   3994.2057 5.822344e+03
#> 26243    172.3183 1.363955e+01
#> 26244  74971.9623 7.522371e+04
#> 26245  82879.6271 2.325024e+04
#> 26246 143010.5935 1.061381e+04
#> 26247  10743.7826 2.134918e+03
#> 26248    823.5848 1.939155e+02
#> 26249   8580.9290 2.617359e+03
#> 26250    146.0542 7.546455e+00
```
