#!/usr/bin/env python3
"""
Iso-chromatic energy-quality assessment of rotary microwave (MW) vs hot-air (HA)
roasting of Chiang Rai Robusta coffee  --  fully reproducible analysis script.

All raw data are embedded below (no external files needed). Running this script
regenerates every table and figure of the manuscript into ./mw_roasting_outputs/.

Design : 2 roasters (HA, MW) x 3 roast degrees (Light, Medium, Dark), n = 3 batches.
Batch  : 100 g green Robusta per batch (all conditions).
Basis  : caffeine, TPC, DPPH and proximate values refer to drip brews (90 C water).
GC-MS  : HS-SPME-GC-MS (one injection per condition); the table embedded at the end
         is the cross-sample aligned, curated feature table (see gcms_curation.py).
Requires: numpy, pandas, scipy, statsmodels, matplotlib  (all open-source)
"""
# =============================================================================
#  CONSTANTS AND ASSUMPTIONS
# =============================================================================
M0  = 100.0    # g green coffee per batch
CP  = 1.5      # kJ kg-1 K-1, mean specific heat of coffee beans (sensitivity 1.3-1.7)
HFG = 2257.0   # kJ kg-1, latent heat of water vaporisation (roast loss treated as water: upper bound)
EF  = 0.4999   # kg CO2e kWh-1, Thailand grid, TGO CFO emission factor (grid mix 2016-2018; in force from Apr 2022); sensitivity 0.4750 (TGO 2026 update) and 0.40-0.50

# =============================================================================
#  RAW DATA  (replicate index i = batch i; identical index across variables = same batch)
# =============================================================================
RAW = {
 "Tend": {"HA-L":[180,180,180], "HA-M":[190,190,190], "HA-D":[210,210,210],
          "MW-L":[168.8,169.4,169.1], "MW-M":[192.8,191.1,191.8], "MW-D":[201.8,202.3,201.0]},   # final bean T, C
 "T0":   {"HA-L":[25.30,25.70,24.90], "HA-M":[25.80,26.50,25.70], "HA-D":[25.38,25.40,25.57],
          "MW-L":[25.73,25.88,24.89], "MW-M":[25.50,26.40,25.20], "MW-D":[25.70,26.80,25.80]},     # initial bean T, C
 "kWh":  {"HA-L":[0.276,0.290,0.284], "HA-M":[0.334,0.348,0.327], "HA-D":[0.456,0.466,0.446],
          "MW-L":[0.168,0.169,0.158], "MW-M":[0.194,0.191,0.198], "MW-D":[0.268,0.273,0.255]},      # metered electricity
 "mass": {"HA-L":[84.12,84.32,84.55], "HA-M":[82.94,81.99,83.01], "HA-D":[81.65,81.69,81.87],
          "MW-L":[88.80,87.80,89.20], "MW-M":[86.03,86.50,85.55], "MW-D":[84.39,85.26,84.12]},      # roasted mass, g
 "time": {"HA-L":[10,11,10], "HA-M":[13,12,13], "HA-D":[16,15,17],
          "MW-L":[12,11,13], "MW-M":[16,15,16], "MW-D":[21,20,21]},                                 # min
 "L":    {"HA-L":[30.69,29.78,30.64], "HA-M":[27.43,28.49,27.81], "HA-D":[21.31,22.65,22.05],
          "MW-L":[29.83,30.29,30.21], "MW-M":[28.23,27.23,25.72], "MW-D":[21.32,22.49,22.70]},
 "a":    {"HA-L":[11.26,11.12,10.71], "HA-M":[9.14,9.59,9.86], "HA-D":[7.63,7.07,7.04],
          "MW-L":[12.26,10.71,11.61], "MW-M":[9.00,9.44,8.69], "MW-D":[7.37,7.80,7.54]},
 "b":    {"HA-L":[17.16,17.98,18.45], "HA-M":[15.08,15.12,14.55], "HA-D":[11.45,11.17,10.03],
          "MW-L":[18.53,18.59,17.59], "MW-M":[14.06,14.44,14.73], "MW-D":[10.61,10.73,10.09]},
 "caffeine": {"HA-L":[1179.20,1168.20,1175.10], "HA-M":[1341.85,1325.80,1338.75], "HA-D":[1449.80,1430.50,1438.42],
              "MW-L":[1160.14,1172.14,1166.28], "MW-M":[1236.12,1247.52,1241.10], "MW-D":[1579.85,1565.32,1573.35]},  # mg kg-1 brew
 "DPPH": {"HA-L":[1.89,1.81,1.82], "HA-M":[1.66,1.60,1.56], "HA-D":[1.59,1.50,1.51],
          "MW-L":[1.92,1.99,1.93], "MW-M":[1.36,1.32,1.35], "MW-D":[2.36,2.30,2.33]},               # mg VCE g-1 brew
 "TPC":  {"HA-L":[0.86,0.82,0.81], "HA-M":[1.38,1.31,1.33], "HA-D":[1.39,1.49,1.42],
          "MW-L":[1.24,1.21,1.29], "MW-M":[1.18,1.11,1.13], "MW-D":[1.89,1.80,1.82]},               # mg GAE g-1 brew
 "carb": {"HA-L":[1.14,1.22,1.18], "HA-M":[1.42,1.56,1.49], "HA-D":[1.16,1.24,1.15],
          "MW-L":[1.51,1.46,1.48], "MW-M":[1.46,1.40,1.35], "MW-D":[1.61,1.66,1.56]},               # g/100 g brew
 "protein": {"HA-L":[0.35,0.43,0.38], "HA-M":[0.23,0.21,0.20], "HA-D":[0.44,0.52,0.48],
             "MW-L":[0.51,0.59,0.56], "MW-M":[0.64,0.61,0.63], "MW-D":[0.78,0.71,0.80]},
 "fat":  {"HA-L":[0.05,0.06,0.07], "HA-M":[0.065,0.065,0.085], "HA-D":[0.08,0.09,0.10],
          "MW-L":[0.08,0.06,0.07], "MW-M":[0.11,0.09,0.10], "MW-D":[0.11,0.10,0.11]},
 "ash":  {"HA-L":[0.43,0.39,0.45], "HA-M":[0.46,0.41,0.43], "HA-D":[0.47,0.41,0.45],
          "MW-L":[0.29,0.34,0.31], "MW-M":[0.30,0.34,0.29], "MW-D":[0.41,0.37,0.38]},
 "moist":{"HA-L":[96.94,98.92,97.91], "HA-M":[97.82,96.65,98.85], "HA-D":[97.20,97.76,98.40],
          "MW-L":[97.85,98.88,95.90], "MW-M":[97.81,95.94,98.80], "MW-D":[97.80,95.40,98.10]},
}

# =============================================================================
#  CURATED, CROSS-SAMPLE-ALIGNED HS-SPME-GC-MS FEATURE TABLE (areas = deconvoluted peak area, counts)
# =============================================================================
VOLATILES = """rt|name|class|formula|mf_max|HA-L|HA-M|HA-D|MW-L|MW-M|MW-D
2.797|Unidentified (RT 2.80 min)|Unidentified|–|81.9|3346024.4|0.0|0.0|0.0|0.0|2606124.4
2.889|Unidentified (RT 2.89 min)|Unidentified|–|94.2|1661903.1|0.0|497840.2|0.0|0.0|1383332.7
2.954|Unidentified (RT 2.95 min)|Unidentified|–|92.9|0.0|0.0|0.0|0.0|438302.3|0.0
3.36|Unidentified (RT 3.36 min)|Unidentified|–|90.4|0.0|0.0|0.0|0.0|0.0|75098.6
3.486|Unidentified (RT 3.49 min)|Unidentified|–|83.2|0.0|0.0|0.0|0.0|0.0|134892.3
3.904|Unidentified (RT 3.90 min)|Unidentified|–|91.6|150287.9|0.0|0.0|131741.8|0.0|0.0
3.904|Unidentified (RT 3.90 min)|Unidentified|–|88.9|0.0|325526.9|0.0|0.0|0.0|0.0
3.921|Unidentified (RT 3.92 min)|Unidentified|–|91.5|0.0|0.0|0.0|0.0|151519.3|0.0
3.937|1,4-Dioxin, 2,3-dihydro-5-methyl-|Others|C5H8O2|93.3|1594467.5|324124.6|0.0|0.0|0.0|0.0
3.965|Unidentified (RT 3.96 min)|Unidentified|–|84.2|0.0|0.0|0.0|0.0|0.0|413711.8
3.985|Unidentified (RT 3.99 min)|Unidentified|–|86.0|0.0|0.0|0.0|582467.6|219038.7|0.0
4.03|Unidentified (RT 4.03 min)|Unidentified|–|81.3|505033.2|511610.7|0.0|393388.1|0.0|0.0
4.061|Unidentified (RT 4.06 min)|Unidentified|–|87.9|1266531.6|0.0|0.0|107026.9|0.0|0.0
4.094|Unidentified (RT 4.09 min)|Unidentified|–|93.4|0.0|0.0|0.0|729483.8|0.0|0.0
4.306|Unidentified (RT 4.31 min)|Unidentified|–|87.2|454256.2|0.0|0.0|0.0|0.0|0.0
4.318|Unidentified (RT 4.32 min)|Unidentified|–|82.7|0.0|0.0|0.0|190620.5|0.0|0.0
4.383|Unidentified (RT 4.38 min)|Unidentified|–|88.2|0.0|254478.4|0.0|0.0|0.0|0.0
4.508|Unidentified (RT 4.51 min)|Unidentified|–|92.7|806676.0|0.0|523683.9|0.0|0.0|0.0
4.555|Unidentified (RT 4.55 min)|Unidentified|–|92.5|389606.4|0.0|0.0|0.0|0.0|0.0
4.59|Unidentified (RT 4.59 min)|Unidentified|–|93.9|0.0|0.0|0.0|0.0|457933.5|0.0
4.604|Unidentified (RT 4.60 min)|Unidentified|–|95.4|459381.4|327681.2|343574.6|0.0|0.0|0.0
4.616|Unidentified (RT 4.62 min)|Unidentified|–|80.5|0.0|0.0|912456.5|0.0|767773.3|0.0
5.657|N-(Pyridin-2-yl)formamide|Pyridines|C6H6N2O|83.7|0.0|0.0|432466.6|0.0|0.0|0.0
5.668|cis-3-Methylenebicyclo[3.3.0]octan-2-one|Ketones|C9H12O|89.3|0.0|711218.1|666283.3|0.0|0.0|0.0
5.683|3-methyl-5-propan-2-yl-pyridazine|Other N-compounds|C8H12N2|98.0|0.0|0.0|0.0|2636268.7|861013.8|0.0
5.694|D-Limonene|Terpenes|C10H16|88.4|0.0|5349976.1|2983754.6|11016590.0|17706294.0|14595770.6
5.747|Cyanogen|Other N-compounds|C2N2|82.3|0.0|0.0|0.0|1482237.0|0.0|5303664.5
5.775|D-Limonene|Terpenes|C10H16|90.5|9905698.2|0.0|0.0|0.0|0.0|0.0
5.777|Pyridine|Pyridines|C5H5N|97.2|3792261.2|5675021.9|6213447.8|6128124.2|4466813.2|17260829.8
5.786|2-Propynenitrile|Other N-compounds|C3HN|89.0|0.0|0.0|0.0|0.0|883973.7|0.0
5.837|1H-pyrrole-3-carboxylic acid tert-butyl ester|Pyrroles|C9H13NO2|82.7|0.0|0.0|227551.6|0.0|0.0|0.0
5.938|3-(exo)-Methyl-(exo)-tricyclo[3.2.1.0(2,4)]octane|Hydrocarbons|C9H14|84.1|0.0|0.0|0.0|0.0|234991.5|0.0
6.16|1H-Pyrrole, 2-methyl-|Pyrroles|C5H7N|83.7|0.0|0.0|285701.7|353942.3|557388.0|0.0
6.307|2-Azabicyclo[2.2.1]heptane|Other N-compounds|C6H11N|87.3|0.0|0.0|0.0|0.0|0.0|1367845.0
6.38|(2E)-5-hydroxy-2-methyl-2-pentenenitrile|Other N-compounds|C6H9NO|87.0|0.0|0.0|0.0|0.0|0.0|290567.2
6.589|Thiazole, 2-methyl-|Sulfur compounds|C4H5NS|81.6|0.0|0.0|0.0|0.0|0.0|104717.4
6.974|3(2H)-Furanone, dihydro-2-methyl-|Furans|C5H8O2|90.8|1612369.7|1801221.9|0.0|1561274.7|1549892.4|408568.4
6.995|2-Iminomalonamide|Other N-compounds|C3H5N3O2|92.5|0.0|0.0|850745.3|0.0|0.0|183196.3
7.011|Malonaldehyde|Aldehydes|C3H4O2|81.9|0.0|464249.0|0.0|0.0|0.0|0.0
7.138|(2S)-1-(2-furanyl)-2-hydroxy-1-propanone|Furans|C7H8O3|80.0|373694.2|0.0|0.0|0.0|0.0|255728.9
7.155|Pyrazine, methyl-|Pyrazines|C5H6N2|97.2|7648223.9|1393548.1|6633275.8|7377100.6|6706655.2|4267930.5
7.164|4-Methylene-1,3-dioxaspiro[4.4]nonan-2-one|Ketones|C8H10O3|80.0|0.0|0.0|0.0|0.0|1872926.4|0.0
7.165|1-Phenoxy-3-vinylcyclopropane|Others|C11H12O|86.2|1860885.6|1380574.2|0.0|1704393.4|0.0|0.0
7.288|2-Cyclohexene-1-carboxylic acid|Acids|C7H10O2|84.8|0.0|0.0|0.0|229315.4|257117.7|169920.0
7.416|2-Furanpropanoic acid, tetrahydro-.alpha.-methyl-, methyl ester|Furans|C9H16O3|87.0|0.0|0.0|463170.6|0.0|338786.5|0.0
7.534|1H-Imidazole-4-carboxaldehyde|Aldehydes|C4H4N2O|82.0|0.0|0.0|402538.6|0.0|0.0|0.0
7.535|2,2-dimethyl-3-(3-methylpent-4-enyl)oxirane|Others|C10H18O|81.7|0.0|0.0|0.0|0.0|0.0|622788.9
7.77|6-Heptyl-5,6-dihydro-2H-pyran-2-one|Pyranones & lactones|C12H20O2|80.9|0.0|0.0|0.0|0.0|147079.9|0.0
7.782|2-Propanone, 1-hydroxy-|Ketones|C3H6O2|94.0|1895090.4|751069.3|805938.0|786137.0|1002209.7|271398.6
7.947|5-Methoxy-2-pentyne|Others|C6H10O|88.1|115727.8|359484.4|236353.2|142676.0|308515.9|0.0
8.376|1,4,7,7-tetramethyl-3-bicyclo[2.2.1]heptanol|Alcohols|C11H20O|89.3|0.0|0.0|1069182.8|0.0|1069524.8|0.0
8.394|Dimethylpyrazine isomer (C6H8N2)|Pyrazines|C6H8N2|91.9|46066608.3|11181939.3|12559836.2|17873382.8|15476738.4|14727416.8
8.41|Pyrazine|Pyrazines|C4H4N2|85.9|0.0|986218.6|513064.4|0.0|0.0|0.0
8.427|(S)-1-(2-Pyridyl)-2-methylpropylamine|Other N-compounds|C9H14N2|96.8|2965729.9|1860417.6|2170362.8|2692578.2|0.0|1634424.7
8.562|1-Propanol, 2-methoxy-|Alcohols|C4H10O2|87.3|228623.2|265096.3|154799.0|0.0|0.0|0.0
8.596|2-Butenoic acid, 2-methyl-, propyl ester, (E)-|Esters|C8H14O2|81.9|108736.0|0.0|0.0|0.0|133220.5|0.0
8.824|1-(hydroxymethyl)-7,7-dimethyl-norbornan-2-ol|Alcohols|C10H18O2|88.4|0.0|0.0|424838.0|0.0|568592.1|0.0
8.832|Pyrazine, 2,3-dimethyl-|Pyrazines|C6H8N2|89.2|2210982.8|1065462.4|0.0|1531497.7|1063495.2|1052501.1
8.842|3-Methyl-1H-1,2,4-triazol-5-yl Azide|Alcohols|C3H4N6|90.5|502554.5|0.0|0.0|0.0|375767.7|0.0
8.898|2-Hexanol, 3,4-dimethyl-|Alcohols|C8H18O|81.3|0.0|442507.4|190495.3|0.0|210818.7|0.0
8.916|(rac)-2-Methylbut-3-ene-1,2-diol|Alcohols|C5H10O2|95.0|166832.2|360254.9|418534.9|372006.9|319813.3|466503.9
9.053|(1S)-1-(1,3-dioxolan-2-yl)ethanamine|Other N-compounds|C5H11NO2|86.6|0.0|381461.7|390782.4|0.0|0.0|330600.0
9.267|(1SR,2RS,4SR,5RS)-Bicyclo[3.1.0]hexane-2,4-diol|Alcohols|C6H10O2|82.8|0.0|0.0|257703.0|0.0|336286.3|0.0
9.283|4-Methyl-1-hexen-3-ol|Alcohols|C7H14O|84.3|485437.7|448775.5|534100.7|0.0|172093.8|141466.7
9.288|2,2-dimethyl-4-oxidanyl-hexan-3-one|Ketones|C8H16O2|86.6|371574.4|0.0|366381.8|189328.6|0.0|0.0
9.489|5-Oxo-endo-tricyclo[5.2.1.0(2,6)]dec-8-en-2-carbonitrile|Other N-compounds|C11H11NO|81.9|0.0|0.0|0.0|221278.9|0.0|0.0
9.492|Pyrazine, 2-ethyl-6-methyl-|Pyrazines|C7H10N2|95.1|0.0|0.0|2405463.1|4166737.9|3807458.4|1921469.0
9.498|1,3,3a,4,5,6,7,8,9,9a-Decahydro-4,9-epoxybenzo[f]isobenzouran-1,3-dione|Ketones|C12H12O4|92.2|0.0|1118037.9|0.0|0.0|1187749.5|0.0
9.503|Pyrazine, 2-ethyl-6-methyl-|Pyrazines|C7H10N2|95.2|3709158.2|3296329.9|0.0|0.0|3375159.4|0.0
9.697|2-oxatertracyclo[4,3,1,1(4,8),o(1,4)]undecane|Others|C10H16O|94.3|160046.6|165667.3|295735.0|1019030.9|269957.5|0.0
9.977|2-(2-Propen-1-ylidene)aminoacetonitrile|Other N-compounds|C5H6N2|81.3|349140.5|290136.0|0.0|0.0|0.0|0.0
10.013|(S)-1-(2-Pyridyl)pentylamine|Other N-compounds|C10H16N2|87.6|0.0|0.0|0.0|296095.3|0.0|0.0
10.03|Pyrazine, 2-ethyl-3-methyl-|Pyrazines|C7H10N2|88.4|0.0|0.0|4013737.2|0.0|0.0|3819096.8
10.033|3a,4,5,6-Tetrahydro-2(3H)-pentalenone|Ketones|C8H10O|94.9|3847841.6|4084098.8|2793109.8|4350559.1|3516274.3|2033541.2
10.23|(4E)-1-hepta-4,6-dienol|Alcohols|C7H12O|89.7|137910.5|0.0|145498.5|123088.1|0.0|172990.5
10.367|Pyridine, 2-ethyl-5-methyl-|Pyridines|C8H11N|91.9|0.0|0.0|182965.0|3231028.4|3057733.1|2765673.0
10.744|(E)-4-hydroxy-2-hexenal|Aldehydes|C6H10O2|80.2|0.0|207470.5|0.0|85164.1|110847.0|122593.2
10.766|Disulfide, dipentyl|Sulfur compounds|C10H22S2|80.6|0.0|0.0|0.0|0.0|0.0|85873.8
10.863|Acetic acid|Acids|C2H4O2|98.7|3266160.0|3002458.7|2095641.5|3160131.7|2425854.2|617322.7
11.181|Pyrazine, 2,6-diethyl-|Pyrazines|C8H12N2|89.5|0.0|0.0|0.0|0.0|0.0|1107987.5
11.255|methyl(1-propylbutyl)amine|Other N-compounds|C8H19N|82.6|0.0|0.0|0.0|0.0|0.0|732934.8
11.262|Oxiraneethanol, 2,3-dimethyl-, cis-(+)-|Alcohols|C6H12O2|94.2|0.0|0.0|6066677.6|0.0|0.0|1975892.4
11.265|2-Propanone, 1-(acetyloxy)-|Ketones|C5H8O3|93.8|5174026.5|8531622.9|7791182.0|6070634.2|6685906.5|2655843.3
11.29|[(1S)-1-cyclopent-2-enyl]methanol|Alcohols|C6H10O|82.0|0.0|209527.6|0.0|0.0|0.0|0.0
11.299|Furfural|Furans|C5H4O2|97.1|9003969.9|6209521.6|3994778.6|3911484.9|3748671.3|728633.4
11.459|1-(3-methyl-1H-pyrazol-4-yl)ethanone|Ketones|C6H8N2O|82.2|0.0|0.0|0.0|0.0|0.0|104796.0
11.572|Furan, 2-[(methylthio)methyl]-|Sulfur compounds|C6H8OS|86.3|115504.3|344447.9|634651.9|479471.5|726793.1|1856684.9
11.807|[(E)-1-butyl-2-hex-2-enylidene]cyanamide|Other N-compounds|C11H18N2|87.4|192405.8|0.0|0.0|0.0|191895.0|143879.9
11.819|(3aR,7aS)-5-methyl-2,3,3a,4,7,7a-hexahydroinden-1-one|Ketones|C10H14O|91.0|0.0|434874.3|425047.9|487449.0|0.0|0.0
11.838|Furfuryl formate|Furans|C6H6O3|86.1|1716260.4|1999334.6|1449588.9|2016713.1|0.0|0.0
11.844|1H-Imidazol-1-ylacetic acid|Acids|C5H6N2O2|91.4|1060926.6|988218.3|536618.0|345331.0|735783.9|161618.4
11.845|2-Furylmethyl formate|Others|C6H6O3|85.0|466270.7|550054.2|0.0|0.0|0.0|0.0
11.901|{1-[2-(Methylamino)(phenyl)methyl]cyclopropyl}methanol|Alcohols|C12H17NO|82.2|457526.1|200915.9|410606.4|0.0|165746.3|0.0
12.005|2-Cyclopenten-1-one, 2,3,4-trimethyl-|Ketones|C8H12O|90.5|0.0|0.0|147356.1|106348.7|0.0|250580.0
12.22|Oct-7-en-1-ol|Alcohols|C8H16O|83.6|0.0|0.0|0.0|0.0|722720.0|651443.8
12.223|Ethanone, 1-(2-furanyl)-|Furans|C6H6O2|92.3|999308.9|0.0|2119733.8|0.0|2091775.0|2665037.7
12.237|(1E)-1-butenyl acetate|Esters|C6H10O2|83.0|0.0|0.0|238982.7|0.0|189657.2|299119.4
12.245|2,5-Hexanedione|Ketones|C6H10O2|94.5|214256.6|553551.3|416571.8|138697.0|0.0|0.0
12.304|4-Ethyl-2,5,6-Trimethylpyrimidine|Other N-compounds|C9H14N2|96.0|183254.0|169199.2|131730.4|160913.8|220925.1|121175.5
12.451|Bicyclo[3.3.1]nonane|Hydrocarbons|C9H16|96.9|0.0|542043.6|359889.1|0.0|416560.9|332579.8
12.457|1-(2-Furyl)-2-propanone|Ketones|C7H8O2|94.3|253136.9|0.0|482651.5|470662.0|550859.2|0.0
12.46|3,3-dimethyl-4-oxopentanenitrile|Other N-compounds|C7H11NO|82.5|0.0|0.0|0.0|124958.1|0.0|0.0
12.615|(R)-(+)-3-Heptene-6-olide|Pyranones & lactones|C7H10O2|80.5|0.0|0.0|197531.4|0.0|0.0|178311.0
12.675|Propanoic acid, anhydride|Acids|C6H10O3|88.0|0.0|1649937.7|0.0|0.0|0.0|0.0
12.732|2-Furanmethanol, acetate|Furans|C7H8O3|96.7|6007411.4|13257936.0|14684087.7|12471538.2|15222686.9|10998445.2
12.76|1H-Imidazole-4-carboxaldehyde|Aldehydes|C4H4N2O|82.2|0.0|0.0|147295.3|0.0|120208.2|188805.6
12.949|1,6-Octadien-3-ol, 3,7-dimethyl-, formate|Alcohols|C11H18O2|83.8|234018.1|0.0|0.0|246025.3|133290.5|0.0
13.115|(2Z,4E)-3,6-Dimethylhepta-2,4-dienamide|Other N-compounds|C9H15NO|82.4|0.0|402008.3|0.0|0.0|317084.5|0.0
13.116|2-Cyclopenten-1-one, 2,3-dimethyl-|Ketones|C7H10O|89.0|0.0|259178.0|275044.1|197712.3|343017.7|862258.8
13.629|1,1-(Ethylenedioxy)-3-hydroxy-5-(methoxymethyl)-5-methylcycloheptane|Others|C12H22O4|86.0|99353.6|152103.6|0.0|0.0|138718.5|202128.9
13.738|2-methyl-5-propyl-furan|Furans|C8H12O|95.9|0.0|739889.4|626158.6|0.0|0.0|2311134.4
13.746|(S)-(+)-4-sec-Butylpyrazole|Other N-compounds|C7H12N2|98.7|427012.9|730487.3|790007.8|646847.2|685855.4|654445.3
13.789|Hydrogen isocyanate|Other N-compounds|CHNO|90.8|401274.1|673819.2|0.0|0.0|0.0|0.0
13.809|Propanoic acid, 2-hydroxy-, ethyl ester|Esters|C5H10O3|84.9|1071676.7|0.0|959642.6|1351272.8|903886.7|0.0
13.846|6-Oxatricyclo[3.2.1.02,7]octan-8-ol, 8-(1,1-dimethylethyl)-, (1.alpha.,2.beta.,5.alpha.,7.beta.,8R*)-|Alcohols|C11H18O2|80.8|165251.0|0.0|0.0|0.0|0.0|0.0
13.846|3-Methoxycyclohex-2-en-1-one|Ketones|C7H10O2|85.0|221317.6|0.0|224118.2|0.0|236062.4|0.0
13.857|2,3,6-Heptanetrione|Ketones|C7H10O3|90.2|0.0|323964.7|223733.6|122058.0|295460.0|182392.5
14.057|2-Furanmethanol, propanoate|Furans|C8H10O3|83.0|0.0|462083.2|727141.7|0.0|540100.6|601095.0
14.097|5-methyl-1H-benzimidazole|Other N-compounds|C8H8N2|88.6|0.0|0.0|0.0|0.0|0.0|81700.2
14.175|(3S)-1-hepten-3-ol|Alcohols|C7H14O|88.5|111698.9|333702.9|320369.9|309529.8|325254.6|178697.6
14.197|Ethanone, 1-(1H-pyrrol-3-yl)-|Pyrroles|C6H7NO|83.7|0.0|0.0|0.0|0.0|0.0|164388.0
14.258|Propanoic acid, 2-methyl-, anhydride|Acids|C8H14O3|91.1|0.0|0.0|0.0|0.0|119776.1|118882.5
14.285|Furan, 2,2'-methylenebis-|Furans|C9H8O2|95.1|147437.0|784178.6|1246034.3|747826.2|1112001.6|2896520.5
14.374|2-Cyclopenten-1-one, 3-isopropyl-4-methyl-|Ketones|C9H14O|87.8|0.0|0.0|0.0|0.0|0.0|95524.4
14.515|1-Ethyl-1H-pyrrole-2-carbaldehyde|Pyrroles|C7H9NO|92.8|580905.0|457623.8|511270.4|767668.7|1128211.0|199867.4
14.547|5-Hydroxymethyl-2-furaldehyde|Aldehydes|C6H6O3|83.0|0.0|173344.2|160247.7|253886.9|130842.8|138356.2
14.714|1-(3-methyl-1H-pyrazol-4-yl)ethanone|Ketones|C6H8N2O|87.0|0.0|0.0|0.0|0.0|374307.3|424641.7
14.719|2-keto-5,5-dimethyl-tetrahydrofuran-3-carbonitrile|Furans|C7H9NO2|90.6|0.0|0.0|139440.7|0.0|0.0|0.0
14.784|5H-5-Methyl-6,7-dihydrocyclopentapyrazine|Pyrazines|C8H10N2|85.3|496217.4|553679.3|480975.3|528031.7|424691.1|488627.4
14.825|(4S)-4-methyl-2-oxolanone|Ketones|C5H8O2|85.4|0.0|0.0|0.0|132719.0|0.0|0.0
14.838|1-Methyl-2-pyrrolecarboxaldehyde|Pyrroles|C6H7NO|97.1|1331924.1|1572818.7|1478312.7|1953041.0|2069402.7|3245629.0
15.023|1-(2-Pyrazinyl)ethanone|Ketones|C6H6N2O|94.2|299795.5|287773.9|245394.0|287780.5|237652.8|82535.9
15.048|3-(Hydroxymethyl)-2-cyclopentenone|Ketones|C6H8O2|81.0|0.0|0.0|0.0|0.0|0.0|96508.1
15.113|1-Octene-3,4-dione|Ketones|C8H12O2|90.0|0.0|0.0|0.0|0.0|0.0|1494289.3
15.116|Butyrolactone|Pyranones & lactones|C4H6O2|95.6|2210749.9|3166968.5|3409225.5|2878021.5|3216612.6|4762378.2
15.221|4-isopropylcyclohex-3-en-1-one|Ketones|C9H14O|86.8|0.0|0.0|448919.7|0.0|0.0|0.0
15.223|2-Butanone, 4-(2-furanyl)-|Furans|C8H10O2|84.2|372074.6|543009.4|312666.7|420100.5|578364.1|672669.7
15.229|1-Oxaspiro[4.5]dec-2-ene|Others|C9H14O|90.7|0.0|242350.6|263572.0|130420.5|0.0|168250.0
15.449|N,1-Dimethyl-2-cyanopyrrole|Pyrroles|C7H8N2|96.9|379254.2|404183.8|482148.0|480805.8|565162.4|889915.6
15.476|(5E)-4-Methyl-1,5-heptadiene|Hydrocarbons|C8H14|89.2|0.0|536434.4|260821.9|0.0|249296.3|388704.4
15.558|1,3-Dimethyl-2-ethylpyrrole|Pyrroles|C8H13N|96.5|271031.8|621736.8|580157.5|515128.8|580371.6|667419.3
15.578|2-Furanmethanol (furfuryl alcohol)|Furans|C5H6O2|96.5|19779354.9|26854392.2|27643185.9|21222505.0|21740275.2|12637412.7
15.631|Pentanoic acid|Acids|C5H10O2|83.1|6204850.7|4303592.7|3681728.5|3897625.5|3233350.1|2733217.8
15.636|2,2-dimethylpropanoic acid tert-butyl ester|Esters|C9H18O2|83.4|324645.3|0.0|0.0|0.0|0.0|0.0
15.801|Furan, 2-(2-furanylmethyl)-5-methyl-|Furans|C10H10O2|90.3|271312.2|415749.6|831271.3|390599.6|565014.3|1895903.1
15.942|Benzene, pentamethyl-|Hydrocarbons|C11H16|83.1|109215.3|0.0|150147.7|0.0|149706.9|163218.6
15.951|Benzene, [(2-methyl-2-propenyl)oxy]-|Others|C10H12O|92.7|0.0|0.0|185522.2|93511.0|0.0|0.0
16.076|2-Cyclopenten-1-one, 3-ethyl-2-hydroxy-|Ketones|C7H10O2|88.2|362583.2|1127987.0|1348835.7|747729.8|897124.6|934460.8
16.117|Pentanenitrile, 2-(aminomethylene)-3,3-dimethyl-, (Z)-|Other N-compounds|C8H14N2|80.5|102257.5|178016.9|219404.5|166716.2|222290.4|270187.3
16.124|2-[(1S,2S)-2-methylcyclopropyl]aniline|Other N-compounds|C10H13N|86.9|0.0|0.0|0.0|177453.6|189893.3|186008.0
16.157|Furfuryl isovalerate|Furans|C10H14O3|92.9|188059.2|289883.6|441122.9|167312.6|276769.1|546306.2
16.258|1-(5-Methyl-2-pyrazinyl)-1-ethanone|Ketones|C7H8N2O|91.2|770267.3|666019.5|547376.4|639388.0|1230862.9|142619.9
16.451|1-(3,5-Dihydroxyphenyl)ethanone|Ketones|C8H8O3|83.3|0.0|0.0|0.0|0.0|0.0|154584.8
16.492|1-(6-Methyl-2-pyrazinyl)-1-ethanone|Ketones|C7H8N2O|96.7|1148665.2|1057794.9|935672.0|1196409.3|1018575.6|190199.2
16.701|2-butyl-5-methyl-hexa-1,4-diene|Hydrocarbons|C11H20|84.4|0.0|0.0|0.0|0.0|0.0|75003.4
16.872|4-Octene, 2-methyl-, (Z)-|Hydrocarbons|C9H18|86.2|0.0|140754.4|147086.7|171070.8|278952.6|200247.9
16.927|Pyrazinecarboxamide|Pyrazines|C5H5N3O|89.0|1708391.1|2154284.2|1844182.2|1674517.9|1685525.2|1756813.2
17.066|Thiophene, 2-methyl-5-propyl-|Sulfur compounds|C8H12S|80.7|0.0|575719.0|0.0|0.0|0.0|199719.2
17.071|1,7-dioxaspiro[5.5]undec-4-en-3-ol|Alcohols|C9H14O3|83.8|0.0|320898.1|0.0|0.0|311063.0|0.0
17.088|4-Pyrimidinol, 5-methoxy-|Alcohols|C5H6N2O2|90.5|0.0|243872.0|205108.6|256081.7|251009.5|280410.3
17.183|1-Vinyl-8-oxabicyclo[3.2.1]oct-6-en-3-one|Ketones|C9H10O2|85.0|0.0|0.0|0.0|117242.7|142909.4|78222.7
17.222|2,4-Dimethylcyclopent-4-ene-1,3-dione|Ketones|C7H8O2|82.3|0.0|176215.4|729205.7|157608.4|316687.7|1219931.9
17.33|(2R,4S)-4-Methyl-2-phenyl-tetrahydro-pyran|Pyranones & lactones|C12H16O|81.6|0.0|0.0|0.0|0.0|0.0|78633.2
17.372|3,4-Dimethylmaleic Anhydride|Others|C6H6O3|84.1|0.0|478536.4|304365.4|309493.9|374386.3|0.0
17.386|2-Cyclohexen-1-one, 2-methyl-|Ketones|C7H10O|85.7|150543.6|0.0|0.0|0.0|0.0|0.0
17.511|2,5-Furandione, 3-ethyl-4-methyl-|Furans|C7H8O3|87.1|121240.5|454607.5|656699.1|144365.6|453067.4|369698.1
17.524|Acetic acid, (dihydro-4,5,5-trimethyl-2-oxo-3(2H)-furanylidene)-, methyl ester, (E)-|Furans|C10H14O4|88.5|0.0|0.0|0.0|231636.9|0.0|0.0
17.525|2,5-Furandione, 3-ethyl-4-methyl-|Furans|C7H8O3|83.1|0.0|0.0|0.0|384100.4|0.0|0.0
17.652|1-Vinyl-8-oxabicyclo[3.2.1]oct-6-en-3-one|Ketones|C9H10O2|82.2|117274.1|243312.4|361945.7|459565.9|0.0|207194.2
17.654|1,3-syn-3-Amino-4-methyl-1-phenylpentan-1-ol|Alcohols|C12H19NO|82.8|0.0|0.0|157702.2|110383.8|141733.8|161936.6
17.803|6-Methylhexahydrocycloprop[a]pentalen-3a,6-diol|Alcohols|C10H16O2|88.0|305478.4|288459.8|262335.8|186216.9|244381.0|588627.1
17.852|2-(2-methylpropyl)-1-cyclopentanone|Ketones|C9H16O|83.6|0.0|132778.9|0.0|0.0|0.0|0.0
17.862|2(3H)-Furanone|Furans|C4H4O2|84.0|228941.4|262833.5|195619.2|116088.9|0.0|0.0
18.054|5,5-Dimethylhex-3-en-2-one|Ketones|C8H14O|83.9|100396.9|200111.5|192738.4|119466.4|137096.2|127382.2
18.106|Methyl salicylate|Others|C8H8O3|84.9|134384.2|231916.8|141750.3|137667.7|142114.3|0.0
18.106|Pyrido[4,3-d]pyrimidin-5(6H)-one|Ketones|C7H5N3O|88.5|0.0|0.0|0.0|0.0|190920.4|91681.6
18.136|2 - ( propane - 2,3 - dione ) - furan|Furans|C7H6O3|90.5|461974.2|724308.9|543789.3|303834.4|461152.8|318400.8
18.184|4-acetyl-1,2-dihydropyrazol-3-one|Ketones|C5H6N2O2|92.9|103631.3|140536.7|0.0|180032.9|118516.4|164632.9
18.194|5-ethyl-4-methoxy-1H-imidazole|Other N-compounds|C6H10N2O|90.3|0.0|0.0|0.0|106052.6|0.0|0.0
18.329|Methyl nicotinate|Other N-compounds|C7H7NO2|92.1|0.0|146775.7|153438.1|0.0|153385.0|502836.5
18.416|2-Butenoic acid, 3-methyl-|Acids|C5H8O2|96.5|2144800.1|1496923.8|999934.8|1323477.4|998808.8|200582.1
18.445|5-acetyl-1H-pyrimidin-2-one|Ketones|C6H6N2O2|81.8|0.0|0.0|0.0|0.0|0.0|318686.7
18.481|2-Cyclopenten-1-one, 2-hydroxy-3,4-dimethyl-|Ketones|C7H10O2|89.5|97185.9|202769.7|395934.1|194363.0|179526.9|323829.6
18.638|2-tert-butyl-3-methyl-thiophene|Sulfur compounds|C9H14S|85.0|0.0|159097.9|222816.0|0.0|126858.7|112714.8
18.739|Benzoic acid, 2-hydroxy-, ethyl ester|Esters|C9H10O3|93.6|446602.1|456551.6|228409.7|0.0|207729.4|0.0
18.786|(2-Hydroxyphenyl)acetonitrile|Other N-compounds|C8H7NO|82.5|0.0|0.0|0.0|0.0|0.0|77193.3
18.868|5-(3-methylphenoxy)methyl-2-amino-1,3,4-thiadiazoles|Sulfur compounds|C10H11N3OS|80.0|0.0|0.0|0.0|0.0|0.0|287106.3
18.954|Cyclobutanone, 3,3-dimethyl-|Ketones|C6H10O|81.6|0.0|0.0|0.0|0.0|0.0|310827.0
19.055|1H-Pyrrole, 1-(2-furanylmethyl)-|Pyrroles|C9H9NO|94.5|1173692.1|1853003.3|1901478.8|1567517.7|1817065.6|1980870.1
19.28|2-Cyclopenten-1-one, 2-hydroxy-3-methyl-|Ketones|C6H8O2|97.2|537154.0|938472.4|962952.1|637862.8|808390.9|1049354.5
19.349|Ethanone, 1-(3,4-dihydro-2H-pyran-2-yl)-|Pyranones & lactones|C7H10O2|83.9|0.0|176943.6|0.0|95110.9|0.0|193613.1
19.447|2,4-Azetidinedione, 3,3-bis(1-methylethyl)-|Ketones|C9H15NO2|82.7|0.0|164031.3|226546.1|109379.7|225905.7|255100.8
19.66|Furan, 3-phenyl-|Furans|C10H8O|96.0|159275.8|180372.1|194198.1|177896.1|221099.9|178054.9
19.803|Phenol, 2-methoxy-|Phenols|C7H8O2|97.9|749662.9|2239138.8|3227560.5|1832565.2|2738281.9|9285534.3
19.824|6-Amino-1,3,5-triazine-2,4(1H,3H)-dione|Ketones|C3H4N4O2|94.3|243597.3|159945.0|0.0|114698.0|0.0|0.0
20.147|3(2H)-Indolizinone, 1,5,8,8a-tetrahydro-5,7-dimethyl-, cis-(.+-.)-|Ketones|C10H15NO|80.4|0.0|0.0|0.0|0.0|103034.6|0.0
20.168|Propanedinitrile, [2-methyl-2-(2-propenyl)cyclohexyl]-, cis-|Other N-compounds|C13H18N2|94.0|0.0|0.0|0.0|0.0|0.0|90600.1
20.361|1-Ethyl-2,2-dimethylcyclohexane|Hydrocarbons|C10H20|80.2|0.0|0.0|0.0|0.0|0.0|78912.6
20.567|2-Cyclopenten-1-one, 3-ethyl-2-hydroxy-|Ketones|C7H10O2|92.0|219689.8|513995.9|648860.2|320224.3|452271.0|907819.4
20.839|Phenylethyl Alcohol|Alcohols|C8H10O|97.1|834805.2|826232.3|831553.6|541808.7|607616.5|482373.3
20.889|3-Buten-2-one, 4-(2-furanyl)-|Furans|C8H8O2|90.7|0.0|0.0|153630.2|0.0|108503.8|115084.4
21.266|1(3-Methoxybenzyl)-2,3,4,5-Tetrahydro-1H-2-benzazepine|Other N-compounds|C18H21NO|83.4|0.0|161086.3|0.0|0.0|0.0|0.0
21.272|Benzeneacetaldehyde, .alpha.-ethylidene-|Aldehydes|C10H10O|95.4|342245.0|474222.6|257424.0|322404.4|199934.3|0.0
21.479|2-Thienylmethanol|Sulfur compounds|C5H6OS|91.5|109420.5|199863.2|201461.2|137037.2|114462.9|73073.4
21.628|4-Methoxytricyclo[5.2.1.0(2,6)]deca-3,8-diene-3-carbonitrile|Other N-compounds|C12H13NO|81.6|0.0|0.0|0.0|0.0|0.0|155326.6
21.728|Phenol, 4-methoxy-3-methyl-|Phenols|C8H10O2|93.1|0.0|0.0|0.0|0.0|0.0|72830.0
21.798|2,2-Dimethyl-1-(5-methylfuran-2-yl)propan-1-one|Furans|C10H14O2|80.4|0.0|183410.6|142800.1|0.0|0.0|0.0
22.106|2-Propylbuta-2,3-dien-1-ol|Alcohols|C7H12O|82.8|0.0|0.0|0.0|0.0|106421.9|0.0
22.11|Ethanone, 1-(1H-pyrrol-2-yl)-|Pyrroles|C6H7NO|94.0|936886.7|2450920.4|3903491.9|1751750.6|3053270.3|3381708.0
22.123|Maltol|Pyranones & lactones|C6H6O3|91.7|953346.7|0.0|0.0|950034.2|0.0|0.0
22.19|Furan, 2,2'-[oxybis(methylene)]bis-|Furans|C10H10O3|97.0|0.0|398698.8|789036.1|303170.1|560474.0|1050965.8
22.613|o-Cresol|Phenols|C7H8O|87.8|0.0|189487.5|231909.7|178395.8|0.0|167106.5
22.625|1,2,3,7b-Tetrahydro-1,1a-diazacyclopropa[a]naphthalene|Other N-compounds|C9H10N2|89.9|152288.0|228651.6|217644.0|164338.1|208045.1|148608.2
22.657|Phenol|Phenols|C6H6O|95.5|218201.6|673331.1|907198.5|443279.4|687022.3|2449896.5
22.833|(E)-4-(5-methyl-2-furanyl)-3-buten-2-one|Furans|C9H10O2|82.5|0.0|285279.7|146451.1|153188.5|275367.0|96779.5
23.109|Guaiacol, 4-ethyl-|Phenols|C9H12O2|96.7|572902.4|1270744.7|2100582.2|1069495.5|1514899.3|6337424.8
23.123|4-endo-Amino-2-oxabicyclo[3.3.0]oct-7-en-3-one|Ketones|C7H9NO2|96.5|1327956.5|767445.9|1899603.1|0.0|0.0|0.0
23.131|1H-Pyrrole-2-carboxaldehyde|Pyrroles|C5H5NO|97.9|826028.8|1060503.2|0.0|952353.4|947157.2|247056.8
23.224|(1.alpha.,2.beta.,4.beta.,5.alpha.)-6-methyl-6-azatricyclo[3.2.0.0(2,4)]heptane-7-one|Ketones|C7H9NO|82.3|0.0|0.0|0.0|0.0|0.0|91205.3
23.298|Furaneol|Furans|C6H8O3|95.3|803366.0|418866.2|285638.4|293504.0|399414.1|91286.6
23.561|Octanoic acid|Acids|C8H16O2|86.6|132860.3|159119.0|175435.6|0.0|130521.7|0.0
23.697|4-p-Tolyl-5H-furan-2-one|Furans|C11H10O2|83.7|0.0|0.0|0.0|0.0|119659.6|96607.7
23.709|2-(cyanamino)-5-methyl-6-ethylpyridine|Pyridines|C8H11N3|80.9|113346.0|226304.6|249856.0|113359.1|161288.2|84934.6
23.825|2-Pyrrolidinone|Pyrroles|C4H7NO|90.2|169841.3|0.0|142025.2|138070.3|138774.3|120780.8
23.85|2-(2-furanyl)-5-methylpyrazine|Pyrazines|C9H8N2O|81.9|113434.7|150929.6|0.0|134519.4|131759.6|92460.2
24.087|2(3H)-Furanone, 5-acetyldihydro-|Furans|C6H8O3|88.4|0.0|158410.3|143697.4|87852.6|114597.4|0.0
24.116|Carbamic acid, methyl-, 3-methylphenyl ester|Esters|C9H11NO2|83.7|0.0|0.0|0.0|0.0|0.0|83094.8
24.17|4-(Furan-2-yl)aniline|Furans|C10H9NO|90.9|0.0|0.0|0.0|0.0|0.0|155629.2
24.268|p-Cresol|Phenols|C7H8O|86.5|0.0|150036.4|168877.1|109332.2|144745.1|112911.0
24.513|Anhydro - sugar|Others|C5H8O4|91.8|0.0|155656.4|301516.2|151026.8|251239.1|1254126.7
24.753|2-amino-1-(3-pyridinyl)ethanol|Pyridines|C7H10N2O|95.7|173369.5|255049.3|286193.2|244750.3|272546.7|128011.4
25.183|1-ethyl-2H-isoquinolin-3-one|Ketones|C11H11NO|82.8|0.0|0.0|0.0|0.0|0.0|73575.6
25.518|Nonanoic acid|Acids|C9H18O2|97.3|612495.4|1036106.8|1035988.5|976111.7|825790.5|519424.0
25.78|Phenol, 4-ethyl-|Phenols|C8H10O|86.3|0.0|0.0|0.0|0.0|0.0|89445.5
25.78|(Z)-4-(3-Decenyl)phenol|Phenols|C16H24O|81.7|0.0|0.0|0.0|0.0|0.0|62813.7
26.259|2-methoxy-4-vinyl-phenol|Phenols|C9H10O2|96.4|4789648.1|6346490.9|6719104.8|3758978.5|3912927.9|3159982.0
26.37|5-Acetoxymethyl-2-furaldehyde|Aldehydes|C8H8O4|81.3|107437.9|149797.2|0.0|0.0|0.0|0.0
26.848|Ethanone, 1-(2,5-dihydroxyphenyl)-|Ketones|C8H8O3|92.5|230936.2|726010.2|844841.7|393202.0|622881.4|489369.0
26.848|Benzene, 1,4-dimethoxy-2-methyl-|Others|C9H12O2|92.7|0.0|728379.5|0.0|312248.6|296389.3|713746.3
27.16|1-(2-furanylmethyl)-2-pyrrolecarboxaldehyde|Pyrroles|C10H9NO2|95.0|103486.0|144063.2|175226.0|158096.2|172371.8|90235.8
27.573|4H-Pyran-4-one, 2,3-dihydro-3,5-dihydroxy-6-methyl-|Pyranones & lactones|C6H8O4|92.0|125511.6|0.0|0.0|0.0|0.0|0.0
27.798|2(3H)-Furanone, 3-acetyldihydro-3-methyl-|Furans|C7H10O3|80.0|0.0|0.0|0.0|0.0|0.0|81181.2
28.268|Furan, 2,2'-(oxydimethylene)di-|Furans|C10H10O3|84.3|0.0|0.0|168768.5|0.0|0.0|79253.0
29.647|Benzaldehyde, 4-methyl-|Aldehydes|C8H8O|88.9|100695.2|0.0|0.0|0.0|0.0|0.0
30.234|3-Pyridinol|Pyridines|C5H5NO|95.7|0.0|134414.0|185563.9|123071.3|152329.3|191217.6
30.443|Indole|Pyrroles|C8H7N|97.0|111287.7|0.0|148615.3|0.0|0.0|218927.2
31.192|1H-Indole, 4-methyl-|Pyrroles|C9H9N|95.9|0.0|0.0|0.0|0.0|0.0|80382.6
31.579|5-Hydroxymethylfurfural|Furans|C6H6O3|93.0|103775.7|0.0|0.0|0.0|0.0|0.0
37.997|Ethylene oxide heptamer|Others|C14H28O7|83.9|0.0|0.0|0.0|0.0|111946.3|0.0
40.571|Hexaethylene glycol|Alcohols|C12H26O7|81.3|0.0|0.0|0.0|0.0|105568.4|0.0
41.672|1H-Purine-2,6-dione, 3,7-dihydro-1,3,7-trimethyl-|Ketones|C8H10N4O2|91.7|138667.1|208503.0|291070.0|233096.4|252246.7|255660.7
41.934|1,4,7,10,13,16-Hexaoxacyclooctadecane|Others|C12H24O6|86.0|0.0|0.0|0.0|0.0|0.0|66376.1
42.194|Glutamic acid, 2-methyl-N-(phenylmethylene)-, dimethyl ester|Esters|C15H19NO4|82.8|0.0|0.0|0.0|0.0|0.0|568501.6
42.475|2-[2-[2-(2-Ethenoxyethoxy)ethoxy]ethoxy]ethanol|Alcohols|C10H20O5|80.4|91759.9|0.0|0.0|0.0|0.0|0.0
43.991|3-[2-Methoxy-6-(methoxymethoxy)phenyl]-2-propylpropen-2-al|Aldehydes|C15H20O4|82.1|0.0|0.0|0.0|0.0|0.0|979955.4
44.407|1,1-Tetramethylene-3-(methylthio)-5-hydroxy-1,3-dihydroisobenzofuran|Sulfur compounds|C13H16O2S|87.3|0.0|0.0|0.0|0.0|0.0|484955.4
46.373|2-methyl-3-phthalimido-propionitrile|Other N-compounds|C12H10N2O2|80.9|102083.1|0.0|0.0|0.0|0.0|0.0
47.695|17-Methyltetraspiro[3.0.3.0.3.0.3.1]heptadecane|Hydrocarbons|C18H26|88.3|0.0|0.0|0.0|182340.1|0.0|0.0
48.465|Octaethylene glycol|Alcohols|C16H34O9|83.7|0.0|338505.0|0.0|0.0|0.0|0.0
"""

# =============================================================================
#  ANALYSIS
# =============================================================================
import io, os, numpy as np, pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list

OUT = "mw_roasting_outputs"; os.makedirs(OUT, exist_ok=True)
CONDS = ["HA-L", "HA-M", "HA-D", "MW-L", "MW-M", "MW-D"]
DEG = {"L": "Light", "M": "Medium", "D": "Dark"}
COL = {"HA": "#eb6834", "MW": "#2a78d6"}          # validated CVD-safe pair
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 8, "pdf.fonttype": 42, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
                     "savefig.dpi": 900, "legend.frameon": False})
# ---------------- Elsevier figure export ----------------
# Final widths: single column 90 mm, 1.5 column 140 mm, double column 190 mm.
# Each figure is resized so that its final (tight) width equals the target width,
# which keeps all text at its true point size (>= 7 pt), then exported as
# RGB TIFF (LZW, 1000 dpi), vector PDF and a 300 dpi PNG preview.
import io, logging
from PIL import Image
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
def save_elsevier(fig, name, width_mm, out=None, dpi=1000):
    out = out or OUT
    target = width_mm / 25.4
    for _ in range(12):
        fig.canvas.draw()
        w = fig.get_tightbbox(fig.canvas.get_renderer()).width + 0.04
        if 0 <= w - target < 0.002 * target: break
        target_adj = target * 1.0005
        s = target_adj / w
        fig.set_size_inches(fig.get_figwidth() * s, fig.get_figheight() * s)
        fig.tight_layout()
    kw = dict(bbox_inches="tight", pad_inches=0.02)
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=dpi, **kw); buf.seek(0)
    im = Image.open(buf).convert("RGB")
    need = int(round(width_mm / 25.4 * dpi))                   # Elsevier minimum width in pixels
    if im.width < need:                                         # pad with white margin, no resampling
        canvas = Image.new("RGB", (need, im.height), "white"); canvas.paste(im, ((need - im.width) // 2, 0)); im = canvas
    im.save(f"{out}/{name}.tif", dpi=(dpi, dpi), compression="tiff_lzw")
    fig.savefig(f"{out}/{name}.pdf", **kw)
    fig.savefig(f"{out}/{name}.png", dpi=300, **kw)


# ---------------- batch-level data frame ----------------
rows = []
for c in CONDS:
    for i in range(3):
        r = {"cond": c, "roaster": c[:2], "degree": c[-1], "rep": i + 1}
        for v in RAW: r[v] = RAW[v][c][i]
        rows.append(r)
df = pd.DataFrame(rows)
df["RL"]    = (M0 - df.mass) / M0 * 100                       # roast (mass) loss, %
df["SEC_g"] = df.kWh / (M0 / 1000)                            # kWh kg-1 green
df["SEC_r"] = df.kWh / (df.mass / 1000)                       # kWh kg-1 roasted
df["P"]     = df.kWh / (df.time / 60)                         # mean electrical power, kW
df["HR"]    = (df.Tend - df.T0) / df.time                     # mean bean heating rate, C min-1
df["Qs"]    = (M0 / 1000) * CP * (df.Tend - df.T0)            # sensible heat, kJ
df["Ql"]    = ((M0 - df.mass) / 1000) * HFG                   # latent heat (upper bound), kJ
df["eta"]   = (df.Qs + df.Ql) / (df.kWh * 3600) * 100         # apparent 1st-law efficiency, %
df["CO2_r"] = df.SEC_r * EF                                   # kg CO2 kg-1 roasted
df.to_csv(f"{OUT}/batch_level_data.csv", index=False)

VARS = ["time", "Tend", "HR", "kWh", "P", "mass", "RL", "SEC_g", "SEC_r", "eta", "CO2_r",
        "L", "a", "b", "caffeine", "DPPH", "TPC", "carb", "protein", "fat", "ash", "moist"]

def cld(g, v):
    """Tukey HSD + compact letter display by the insert-and-absorb algorithm
    (Piepho, 2004); conditions share a letter if and only if they do not differ
    significantly (alpha = 0.05). Letters are ordered so that 'a' marks the highest mean."""
    t = pairwise_tukeyhsd(g[v], g.cond, alpha=0.05)
    res = pd.DataFrame(t._results_table.data[1:], columns=t._results_table.data[0])
    conds = sorted(g.cond.unique())
    sets = [set(conds)]
    for r in res.itertuples():
        if not r.reject: continue
        new = []
        for L in sets:
            if r.group1 in L and r.group2 in L: new += [L - {r.group1}, L - {r.group2}]
            else: new.append(L)
        new = [L for L in new if not any(L < M for M in new)]
        sets = []
        for L in new:
            if L not in sets: sets.append(L)
    means = g.groupby("cond")[v].mean().sort_values(ascending=False)
    rank = {c: i for i, c in enumerate(means.index)}
    sets.sort(key=lambda L: min(rank[c] for c in L))
    lab = {c: "" for c in conds}
    for k, L in enumerate(sets):
        for c in L: lab[c] += "abcdefgh"[k]
    return lab

summ, anov = [], []
for v in VARS:
    g = df[["cond", "roaster", "degree", v]].rename(columns={v: "y"})
    a = anova_lm(smf.ols("y ~ C(roaster) * C(degree)", g).fit(), typ=2); ss = a.sum_sq
    rec = {"var": v}
    for i, k in enumerate(["Roaster", "Degree", "RxD"]):
        rec[f"{k}_F"] = a.F.iloc[i]; rec[f"{k}_p"] = a["PR(>F)"].iloc[i]
        rec[f"{k}_eta2p"] = ss.iloc[i] / (ss.iloc[i] + ss.iloc[-1])
    anov.append(rec)
    lab = cld(g.rename(columns={"y": v}), v)
    for c in CONDS:
        x = df.loc[df.cond == c, v]
        summ.append({"var": v, "cond": c, "mean": x.mean(), "sd": x.std(ddof=1), "letter": lab[c]})
S = pd.DataFrame(summ); A = pd.DataFrame(anov)
S.to_csv(f"{OUT}/Table_means_sd_tukey.csv", index=False); A.to_csv(f"{OUT}/Table_two_way_anova.csv", index=False)
mean = S.pivot(index="cond", columns="var", values="mean")

# ---------------- TOST equivalence of colour (Welch) ----------------
BOUNDS = {"L": 2.0, "a": 1.5, "b": 1.5}
tost = []
for d in "LMD":
    for v, delta in BOUNDS.items():
        x, y = np.array(RAW[v]["MW-" + d]), np.array(RAW[v]["HA-" + d])
        diff = x.mean() - y.mean(); vx, vy = x.var(ddof=1) / 3, y.var(ddof=1) / 3
        se = np.sqrt(vx + vy); dfw = (vx + vy) ** 2 / (vx ** 2 / 2 + vy ** 2 / 2)
        p = max(stats.t.sf((diff + delta) / se, dfw), stats.t.cdf((diff - delta) / se, dfw))
        h = stats.t.ppf(0.95, dfw) * se
        tost.append({"degree": d, "var": v, "diff_MW_minus_HA": diff, "CI90_lo": diff - h,
                     "CI90_hi": diff + h, "bound": delta, "p_TOST": p})
T = pd.DataFrame(tost); T.to_csv(f"{OUT}/Table_TOST_colour.csv", index=False)
rng = np.random.default_rng(42); dE = []
for d in "LMD":
    point = np.sqrt(sum((mean.loc["MW-" + d, v] - mean.loc["HA-" + d, v]) ** 2 for v in "Lab"))
    bs = [np.sqrt(sum((rng.choice(RAW[v]["MW-" + d], 3).mean() - rng.choice(RAW[v]["HA-" + d], 3).mean()) ** 2
                      for v in "Lab")) for _ in range(5000)]
    dE.append({"degree": d, "dE_ab": point, "boot95_lo": np.percentile(bs, 2.5), "boot95_hi": np.percentile(bs, 97.5)})
pd.DataFrame(dE).to_csv(f"{OUT}/Table_dE_MW_vs_HA.csv", index=False)

# ---------------- ANCOVA: colour vs roast loss (decoupling test) ----------------
anc = smf.ols("L ~ RL + C(roaster)", df).fit()
anc_int = smf.ols("L ~ RL * C(roaster)", df).fit()
with open(f"{OUT}/ANCOVA_L_vs_RL.txt", "w") as f:
    f.write(anc.summary().as_text() + "\n\n" + anc_int.summary().as_text())

# ---------------- volatile block ----------------
FT = pd.read_csv(io.StringIO(VOLATILES), sep="|")
Av = FT[CONDS]; rel = Av / Av.sum() * 100
H = lambda x: -((x[x > 0] / x.sum()) * np.log(x[x > 0] / x.sum())).sum()
div = pd.DataFrame({"richness": (Av > 0).sum(), "shannon_H": Av.apply(H)})
div["pielou_J"] = div.shannon_H / np.log(div.richness)
div.to_csv(f"{OUT}/Table_volatile_diversity.csv")
cls_pct = rel.groupby(FT["class"]).sum()[CONDS]; cls_pct.to_csv(f"{OUT}/Table_class_relative_area.csv")
overlap = []
for d in "LMD":
    a_, b_ = Av["HA-" + d] > 0, Av["MW-" + d] > 0
    overlap.append({"degree": d, "HA": a_.sum(), "MW": b_.sum(), "shared": (a_ & b_).sum(),
                    "HA_only": (a_ & ~b_).sum(), "MW_only": (~a_ & b_).sum(), "Jaccard": (a_ & b_).sum() / (a_ | b_).sum()})
pd.DataFrame(overlap).to_csv(f"{OUT}/Table_volatile_overlap.csv", index=False)

# Thermal-severity volatile index (TSVI): products of trigonelline / phenolic-acid pyrolysis
# relative to early-stage sugar-derived products.
num = cls_pct.loc[["Pyridines", "Phenols"]].sum(); den = cls_pct.loc[["Furans", "Acids", "Aldehydes"]].sum()
TSVI = (num / den).rename("TSVI"); TSVI.to_csv(f"{OUT}/Table_TSVI.csv")

# ---------------- multi-block fusion (MFA-weighted PCA on condition means) ----------------
blocks = {
    "Process":   mean.loc[CONDS, ["SEC_r", "RL", "HR", "time"]],
    "Colour":    mean.loc[CONDS, ["L", "a", "b"]],
    "Bioactive": mean.loc[CONDS, ["caffeine", "TPC", "DPPH"]],
    "Volatile":  pd.concat([cls_pct.loc[["Pyridines", "Phenols", "Furans", "Pyrazines", "Acids", "Sulfur compounds"]].T,
                            div[["richness"]]], axis=1).loc[CONDS],
}
Z, owner = [], []
for bname, B in blocks.items():
    Bz = (B - B.mean()) / B.std(ddof=0)
    s1 = np.linalg.svd(Bz.values, compute_uv=False)[0]
    Z.append(Bz / s1); owner += [bname] * B.shape[1]
Zc = pd.concat(Z, axis=1)
U, s, Vt = np.linalg.svd(Zc.values - Zc.values.mean(0), full_matrices=False)
expl = s ** 2 / (s ** 2).sum() * 100
scores = pd.DataFrame(U[:, :2] * s[:2], index=CONDS, columns=["Dim1", "Dim2"])
corr = pd.DataFrame({k: [np.corrcoef(Zc[c], scores[k])[0, 1] for c in Zc.columns] for k in ["Dim1", "Dim2"]},
                    index=Zc.columns); corr["block"] = owner
scores.to_csv(f"{OUT}/MFA_scores.csv"); corr.to_csv(f"{OUT}/MFA_variable_correlations.csv")

# ---------------- entropy-weighted TOPSIS + Monte-Carlo weight robustness ----------------
crit = pd.DataFrame({"SEC_r": mean.loc[CONDS, "SEC_r"], "RL": mean.loc[CONDS, "RL"],
                     "TPC": mean.loc[CONDS, "TPC"], "DPPH": mean.loc[CONDS, "DPPH"],
                     "Shannon_H": div.loc[CONDS, "shannon_H"]})
benefit = np.array([False, False, True, True, True])
def topsis(X, w):
    N = X / np.sqrt((X ** 2).sum(0)); V = N * w
    best = np.where(benefit, V.max(0), V.min(0)); worst = np.where(benefit, V.min(0), V.max(0))
    dp = np.sqrt(((V - best) ** 2).sum(1)); dn = np.sqrt(((V - worst) ** 2).sum(1))
    return dn / (dp + dn)
X = crit.values
P = X / X.sum(0); E = -(P * np.log(P)).sum(0) / np.log(len(X)); w_ent = (1 - E) / (1 - E).sum()
cc = topsis(X, w_ent)
NDIR = 100000                                                  # Dirichlet weight vectors
mc = np.array([topsis(X, w_) for w_ in rng.dirichlet(np.ones(5), NDIR)])
ranks = (-mc).argsort(1).argsort(1) + 1
rank_acc = pd.DataFrame({f"P(rank={k})": (ranks == k).mean(0) for k in range(1, 7)}, index=CONDS)
topsis_tab = crit.assign(closeness_entropy=cc, rank_entropy=(-cc).argsort().argsort() + 1).join(rank_acc)
pd.Series(w_ent, index=crit.columns, name="entropy_weight").to_csv(f"{OUT}/TOPSIS_entropy_weights.csv")
topsis_tab.to_csv(f"{OUT}/Table_TOPSIS.csv")

# ---------------- sensitivity: efficiency (cp) and carbon (EF) ----------------
sens = []
for cp in (1.3, 1.5, 1.7):
    e = ((M0 / 1000) * cp * (df.Tend - df.T0) + df.Ql) / (df.kWh * 3600) * 100
    for c, v in e.groupby(df.cond).mean().items(): sens.append({"param": "cp", "value": cp, "cond": c, "eta": v})
for ef in (0.40, 0.4750, 0.4999, 0.50):
    for c, v in (df.SEC_r * ef).groupby(df.cond).mean().items(): sens.append({"param": "EF", "value": ef, "cond": c, "CO2_r": v})
pd.DataFrame(sens).to_csv(f"{OUT}/Table_sensitivity.csv", index=False)


# =============================================================================
#  ROBUSTNESS AND UNCERTAINTY ANALYSES (added in revision)
# =============================================================================
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.oneway import anova_oneway
rob = []
for v in VARS:
    g = df[["cond", "roaster", "degree", v]].rename(columns={v: "y"})
    fit = smf.ols("y ~ C(roaster) * C(degree)", g).fit()
    sw = stats.shapiro(fit.resid).pvalue
    lev = stats.levene(*[x["y"].values for _, x in g.groupby("cond")], center="median").pvalue
    welch = anova_oneway(g["y"], g["cond"], use_var="unequal").pvalue          # Welch one-way across 6 conditions
    # aligned-rank-transform style check for the roaster main effect: rank-based two-way ANOVA
    g = g.assign(r=g["y"].rank())
    ar = anova_lm(smf.ols("r ~ C(roaster) * C(degree)", g).fit(), typ=2)
    rob.append({"var": v, "Shapiro_resid_p": sw, "Levene_p": lev, "Welch_oneway_p": welch,
                "rankANOVA_roaster_p": ar["PR(>F)"].iloc[0], "rankANOVA_degree_p": ar["PR(>F)"].iloc[1],
                "rankANOVA_RxD_p": ar["PR(>F)"].iloc[2]})
ROB = pd.DataFrame(rob)
for k in ["Roaster", "Degree", "RxD"]:
    ROB[f"{k}_p_BH"] = multipletests(A[f"{k}_p"].values, method="fdr_bh")[1]
ROB.to_csv(f"{OUT}/Table_S_robustness_assumptions_FDR.csv", index=False)

# ANCOVA restricted to the overlapping roast-loss window and leave-one-degree-out check
lo, hi = df.loc[df.roaster == "HA", "RL"].min(), df.loc[df.roaster == "MW", "RL"].max()
overlap = df[(df.RL >= lo) & (df.RL <= hi)]
anc_light_out = smf.ols("L ~ RL + C(roaster)", df[df.degree != "L"]).fit()
with open(f"{OUT}/ANCOVA_overlap.txt", "w") as f:
    f.write(f"HA RL range {df.loc[df.roaster=='HA','RL'].min():.2f}-{df.loc[df.roaster=='HA','RL'].max():.2f}; "
            f"MW RL range {df.loc[df.roaster=='MW','RL'].min():.2f}-{df.loc[df.roaster=='MW','RL'].max():.2f}\n")
    f.write(f"Overlap window {lo:.2f}-{hi:.2f} %, batches in window: {len(overlap)} ({overlap.cond.tolist()})\n")
    f.write(f"Within window mean L*: " + str(overlap.groupby('roaster').L.mean().round(2).to_dict()) + "\n")
    f.write(f"Within window mean RL: " + str(overlap.groupby('roaster').RL.mean().round(2).to_dict()) + "\n")
    f.write("ANCOVA excluding light roasts: MW coef %.2f, 95%% CI %s, p=%.4g\n" % (
        anc_light_out.params['C(roaster)[T.MW]'], anc_light_out.conf_int().loc['C(roaster)[T.MW]'].round(2).tolist(),
        anc_light_out.pvalues['C(roaster)[T.MW]']))

# Parametric (t-based Monte-Carlo, 10^6 draws) uncertainty of ΔE*ab, replacing the n = 3 bootstrap as primary
rngp = np.random.default_rng(7); dEp = []
for d in "LMD":
    comp = []
    for v in "Lab":
        x, y = np.array(RAW[v]["MW-" + d]), np.array(RAW[v]["HA-" + d])
        comp.append((x.mean() - y.mean()) + stats.t.rvs(2, size=1_000_000, random_state=rngp) * x.std(ddof=1) / np.sqrt(3)
                    - stats.t.rvs(2, size=1_000_000, random_state=rngp) * y.std(ddof=1) / np.sqrt(3))
    de = np.sqrt(sum(c ** 2 for c in comp))
    dEp.append({"degree": d, "dE_point": np.sqrt(sum((np.mean(RAW[v]['MW-'+d]) - np.mean(RAW[v]['HA-'+d])) ** 2 for v in 'Lab')),
                "MC_median": np.median(de), "MC_95_hi": np.percentile(de, 95)})
pd.DataFrame(dEp).to_csv(f"{OUT}/Table_S_dE_parametric.csv", index=False)

# Energy-meter resolution (±0.001 kWh, rectangular) and mass (±0.01 g) propagated to SEC_r
res_kwh, res_g = 0.001 / np.sqrt(3), 0.01 / np.sqrt(3)
unc = []
for c in CONDS:
    k, m_ = mean.loc[c, "kWh"], mean.loc[c, "mass"]
    sec = k / (m_ / 1000)
    u_meas = sec * np.sqrt((res_kwh / k) ** 2 + (res_g / m_) ** 2)
    u_rep = S.query("var=='SEC_r' and cond==@c")["sd"].iloc[0] / np.sqrt(3)
    unc.append({"cond": c, "SEC_r": sec, "u_meter_rel_%": u_meas / sec * 100, "u_repeat_rel_%": u_rep / sec * 100,
                "u_combined_rel_%": np.sqrt(u_meas ** 2 + u_rep ** 2) / sec * 100})
pd.DataFrame(unc).to_csv(f"{OUT}/Table_S_SEC_uncertainty.csv", index=False)

# TOPSIS criteria-set sensitivity: (i) baseline, (ii) roast loss treated as benefit, (iii) roast loss removed
def topsis_set(Xm, ben, w):
    N = Xm / np.sqrt((Xm ** 2).sum(0)); V = N * w
    b_ = np.where(ben, V.max(0), V.min(0)); w_ = np.where(ben, V.min(0), V.max(0))
    dp = np.sqrt(((V - b_) ** 2).sum(1)); dn = np.sqrt(((V - w_) ** 2).sum(1)); return dn / (dp + dn)
def entw(Xm):
    P_ = Xm / Xm.sum(0); E_ = -(P_ * np.log(P_)).sum(0) / np.log(len(Xm)); return (1 - E_) / (1 - E_).sum()
sets = {"baseline (RL cost)": (X, benefit),
        "RL as benefit": (X, np.array([False, True, True, True, True])),
        "RL removed": (X[:, [0, 2, 3, 4]], np.array([False, True, True, True]))}
ts = []
for name, (Xm, ben) in sets.items():
    cc_ = topsis_set(Xm, ben, entw(Xm)); rk = (-cc_).argsort().argsort() + 1
    mc_ = np.array([topsis_set(Xm, ben, w_) for w_ in rng.dirichlet(np.ones(Xm.shape[1]), NDIR)])
    first = (mc_.argmax(1)[:, None] == np.arange(6)).mean(0)
    if name.startswith("baseline"):                       # baseline = the same simulation as Table 6
        first = (mc.argmax(1)[:, None] == np.arange(6)).mean(0)
    for i, c in enumerate(CONDS):
        ts.append({"criteria_set": name, "cond": c, "C_entropy": cc_[i], "rank_entropy": rk[i], "P_rank1": first[i]})
pd.DataFrame(ts).to_csv(f"{OUT}/Table_S_TOPSIS_criteria_sensitivity.csv", index=False)
print("\n=== ROBUSTNESS ==="); print(ROB.round(4).to_string()); print(open(f"{OUT}/ANCOVA_overlap.txt").read())
print(pd.DataFrame(dEp).round(2).to_string()); print(pd.DataFrame(unc).round(2).to_string())
print(pd.DataFrame(ts).pivot(index="cond", columns="criteria_set", values="rank_entropy"))
print(pd.DataFrame(ts).pivot(index="cond", columns="criteria_set", values="P_rank1").round(3))

# =============================================================================
#  FIGURES
# =============================================================================
xpos = {"L": 0, "M": 1, "D": 2}; W = 0.36
def bars(ax, v, ylabel, fmt="{:.2f}", letters=True, ylim0=True):
    for r, off in (("HA", -W / 2 - 0.01), ("MW", W / 2 + 0.01)):
        m = [mean.loc[f"{r}-{d}", v] for d in "LMD"]
        sd = [S.query("var==@v and cond==@r+'-'+@d")["sd"].iloc[0] for d in "LMD"]
        xs = [xpos[d] + off for d in "LMD"]
        ax.bar(xs, m, W, color=COL[r], edgecolor="white", linewidth=1.5, label="Hot-air (HA)" if r == "HA" else "Microwave (MW)", zorder=3)
        ax.errorbar(xs, m, yerr=sd, fmt="none", ecolor=INK, elinewidth=0.9, capsize=2.5, zorder=4)
        if letters:
            for x, mm, s_, d in zip(xs, m, sd, "LMD"):
                lt = S.query("var==@v and cond==@r+'-'+@d")["letter"].iloc[0]
                ax.text(x, mm + s_ + 0.02 * ax.get_ylim()[1], lt, ha="center", va="bottom", fontsize=8, color=INK2)
    ax.set_xticks(range(3)); ax.set_xticklabels([DEG[d] for d in "LMD"]); ax.set_ylabel(ylabel)
    ax.yaxis.grid(True, color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)

def panel(ax, s): ax.text(-0.16, 1.04, s, transform=ax.transAxes, fontsize=9, fontweight="bold", color=INK)

# Fig. 1 – process & energy
fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6))
for ax, v, yl, tag in zip(axs.flat, ["time", "HR", "SEC_r", "eta"],
        ["Roasting time (min)", "Mean bean heating rate (°C min$^{-1}$)",
         "SEC (kWh kg$^{-1}$ roasted coffee)", "Apparent energy efficiency (%)"], "abcd"):
    ax.set_ylim(0, df.groupby("cond")[v].mean().max() * (1.38 if v == "SEC_r" else 1.25)); bars(ax, v, yl); panel(ax, f"({tag})")
for d in "LMD":
    sav = (1 - mean.loc["MW-" + d, "SEC_r"] / mean.loc["HA-" + d, "SEC_r"]) * 100
    axs[1, 0].text(xpos[d], mean.loc["HA-" + d, "SEC_r"] + 1.05, f"MW −{sav:.0f}%", ha="center", fontsize=7.5, color=INK, fontweight="bold")
axs[0, 0].legend(loc="upper left", fontsize=8)
fig.tight_layout(); save_elsevier(fig, "Fig1_process_energy", 190); plt.close(fig)

# Fig. 2 – iso-chromatic decoupling
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
for r, mk in (("HA", "o"), ("MW", "s")):
    g = df[df.roaster == r]; ax.scatter(g.RL, g.L, s=34, marker=mk, color=COL[r], edgecolor="white", linewidth=1.2, zorder=3,
                                        label="Hot-air (HA)" if r == "HA" else "Microwave (MW)")
    sl = stats.linregress(g.RL, g.L); xx = np.linspace(g.RL.min() - .3, g.RL.max() + .3, 20)
    ax.plot(xx, sl.intercept + sl.slope * xx, color=COL[r], lw=2, zorder=2)
for d in "LMD":
    sub = mean.loc[["HA-" + d, "MW-" + d]]
    ax.annotate("", xy=(sub.RL.iloc[1], sub.L.iloc[1]), xytext=(sub.RL.iloc[0], sub.L.iloc[0]),
                arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8, ls="--"))
    ax.text(sub.RL.mean(), sub.L.mean() + 0.9, DEG[d], ha="center", fontsize=8, color=INK2)
ax.set_xlabel("Roast loss (%)"); ax.set_ylabel("Lightness, L*"); ax.legend(loc="lower left", fontsize=8)
ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True); panel(ax, "(a)")
ax = axs[1]; y = 0; yt, yl = [], []
for d in "DML":
    for v in ["b", "a", "L"]:
        t = T.query("degree==@d and var==@v").iloc[0]
        eq = t.p_TOST < 0.05
        ax.plot([-t.bound, -t.bound], [y - .35, y + .35], color=INK2, lw=0.8); ax.plot([t.bound, t.bound], [y - .35, y + .35], color=INK2, lw=0.8)
        ax.plot([t.CI90_lo, t.CI90_hi], [y, y], color=COL["MW"], lw=2.2)
        ax.scatter([t.diff_MW_minus_HA], [y], s=30, color=COL["MW"] if eq else "white", edgecolor=COL["MW"], zorder=3, linewidth=1.4)
        yt.append(y); yl.append(f"{DEG[d]} · {v}*  (p={t.p_TOST:.3f})"); y += 1
    y += 0.6
ax.axvline(0, color=GRID, lw=1, zorder=0); ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=7.5)
ax.set_xlabel("MW − HA difference (90% CI)"); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.set_xlim(-3.3, 3.3); panel(ax, "(b)")
ax.text(0.5, -0.30, "Filled = equivalent (TOST p < 0.05); vertical ticks = equivalence bounds", transform=ax.transAxes, ha="center", fontsize=7, color=INK2)
fig.tight_layout(); save_elsevier(fig, "Fig2_colour_equivalence_offset", 190); plt.close(fig)

# Fig. 3 – bioactives (brew basis)
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.8))
for ax, v, yl, tag in zip(axs, ["caffeine", "TPC", "DPPH"],
        ["Caffeine (mg kg$^{-1}$ brew)", "TPC (mg GAE g$^{-1}$ brew)", "DPPH (mg VCE g$^{-1}$ brew)"], "abc"):
    ax.set_ylim(0, df.groupby("cond")[v].mean().max() * 1.22); bars(ax, v, yl); panel(ax, f"({tag})")
axs[0].legend(loc="upper left", fontsize=7.5)
fig.tight_layout(); save_elsevier(fig, "Fig3_bioactives", 190); plt.close(fig)

# Fig. 4 – volatile class composition + log2 fold change
order_cls = cls_pct.mean(axis=1).sort_values(ascending=False).index.tolist()
main = [c for c in order_cls if c not in ("Others",)][:9]
plot_pct = cls_pct.loc[main].copy(); plot_pct.loc["Other classes"] = 100 - plot_pct.sum()
seq = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948", "#86b6ef", "#c3c2b7"]
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.6), gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]; bottom = np.zeros(6)
for k, (cname, rowv) in enumerate(plot_pct.iterrows()):
    ax.bar(range(6), rowv.values, 0.7, bottom=bottom, color=seq[k], edgecolor="white", linewidth=1.2, label=cname)
    bottom += rowv.values
ax.set_xticks(range(6)); ax.set_xticklabels(CONDS, rotation=45, ha="right", fontsize=8); ax.set_ylabel("Relative peak area (%)")
ax.set_ylim(0, 100); ax.legend(bbox_to_anchor=(1.0, 1.0), loc="upper left", fontsize=7); panel(ax, "(a)")
ax = axs[1]
fc = np.log2((cls_pct.loc[main, [f"MW-{d}" for d in "LMD"]].values + 0.05) / (cls_pct.loc[main, [f"HA-{d}" for d in "LMD"]].values + 0.05))
im = ax.imshow(fc, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_xticks(range(3)); ax.set_xticklabels([DEG[d] for d in "LMD"]); ax.set_yticks(range(len(main))); ax.set_yticklabels(main, fontsize=7.5)
for i in range(fc.shape[0]):
    for j in range(3): ax.text(j, i, f"{fc[i, j]:+.1f}", ha="center", va="center", fontsize=7, color="white" if abs(fc[i, j]) > 1.6 else INK)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.set_label("log$_2$(MW/HA)", fontsize=8); cb.outline.set_visible(False)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.text(-0.75, 1.04, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold")
fig.tight_layout(); save_elsevier(fig, "Fig4_volatile_classes", 190); plt.close(fig)

# Fig. 5 – clustered heatmap of 30 most abundant aligned features (row z-score of log area)
top = rel.assign(t=rel.sum(axis=1)).sort_values("t", ascending=False).head(30).drop(columns="t")
Lg = np.log10(top + 0.01); Zr = Lg.sub(Lg.mean(axis=1), axis=0).div(Lg.std(axis=1, ddof=0) + 1e-9, axis=0)
o = leaves_list(linkage(Zr.values, "average", metric="euclidean"))
Zr = Zr.iloc[o]; labels = [f"{FT.loc[i, 'name'][:38]}  [{FT.loc[i, 'rt']:.2f}]" for i in Zr.index]
fig, ax = plt.subplots(figsize=(6.4, 7.6))
im = ax.imshow(Zr.values, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
ax.set_xticks(range(6)); ax.set_xticklabels(CONDS); ax.xaxis.tick_top()
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.axvline(2.5, color="white", lw=3)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("Row z-score of log$_{10}$ relative area", fontsize=8); cb.outline.set_visible(False)
fig.tight_layout(); save_elsevier(fig, "Fig5_volatile_heatmap", 140); plt.close(fig)

# Fig. 6 – MFA fusion (scores + correlation circle) and TOPSIS robustness
fig = plt.figure(figsize=(7.48, 7.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1.15, 0.85])
axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
ax = axs[0]
for c in CONDS:
    ax.scatter(*scores.loc[c], s=80, marker="o" if c[:2] == "HA" else "s", color=COL[c[:2]], edgecolor="white", linewidth=1.5, zorder=4)
    up = c[:2] == "HA"
    ax.text(scores.loc[c, "Dim1"], scores.loc[c, "Dim2"] + (0.14 if up else -0.14), c, ha="center", va="bottom" if up else "top", fontsize=7.5, color=INK, fontweight="bold")
for r in ("HA", "MW"):
    pts = scores.loc[[f"{r}-{d}" for d in "LMD"]].values
    ax.plot(pts[:, 0], pts[:, 1], color=COL[r], lw=1.2, ls="--", zorder=2)
ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
lim = np.abs(scores.values).max() * 1.3; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
ax.set_xlabel(f"Dim 1 ({expl[0]:.1f}%)"); ax.set_ylabel(f"Dim 2 ({expl[1]:.1f}%)"); panel(ax, "(a)")
ax = axs[1]
bcol = {"Process": INK2, "Colour": "#c98500", "Bioactive": "#199e70", "Volatile": "#6250d6"}
t_ = np.linspace(0, 2 * np.pi, 200); ax.plot(np.cos(t_), np.sin(t_), color=GRID, lw=1)
LBL = {"Sulfur compounds": "S-compounds", "SEC_r": "SEC", "RL": "Roast loss", "HR": "Heating rate", "L": "L*", "a": "a*", "b": "b*", "time": "Time", "richness": "Richness", "caffeine": "Caffeine"}
pos = {}
for var, r in corr.iterrows():
    ax.annotate("", xy=(r.Dim1, r.Dim2), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=bcol[r.block], lw=0.9))
    ang = np.arctan2(r.Dim2, r.Dim1); rr = np.hypot(r.Dim1, r.Dim2) + 0.24
    pos[var] = np.array([rr * np.cos(ang), rr * np.sin(ang)])
for _ in range(300):                                  # simple label repulsion
    keys = list(pos)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            d = pos[keys[i]] - pos[keys[j]]
            wx = 0.045 * (len(LBL.get(keys[i], keys[i])) + len(LBL.get(keys[j], keys[j]))) / 2 + 0.04
            if abs(d[0]) < wx and abs(d[1]) < 0.11:
                push = np.array([0, 0.012 if d[1] >= 0 else -0.012])
                pos[keys[i]] += push; pos[keys[j]] -= push
for var, r in corr.iterrows():
    lab = LBL.get(var, var)
    ax.text(*pos[var], lab, fontsize=7, color=bcol[r.block], ha="center", va="center", zorder=5,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85))
ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
ax.set_xlim(-1.55, 1.55); ax.set_ylim(-1.45, 1.45); ax.set_aspect("equal")
ax.set_xlabel("Correlation with Dim 1"); ax.set_ylabel("Correlation with Dim 2")
ax.legend(handles=[Patch(color=v, label=k) for k, v in bcol.items()], fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=4, handlelength=1); panel(ax, "(b)")
ax = axs[2]; ordc = topsis_tab.sort_values("closeness_entropy").index
left = np.zeros(len(ordc)); shades = ["#0d366b", "#256abf", "#6da7ec", "#b7d3f6", "#e4e3df", "#f0efec"]
for k in range(1, 7):
    vals = rank_acc.loc[ordc, f"P(rank={k})"].values * 100
    ax.barh(range(len(ordc)), vals, left=left, color=shades[k - 1], edgecolor="white", linewidth=1, label=f"Rank {k}"); left += vals
for i, c in enumerate(ordc):
    ax.text(102, i, f"C*={topsis_tab.loc[c, 'closeness_entropy']:.2f}", va="center", fontsize=7.5, color=INK)
ax.set_yticks(range(len(ordc))); ax.set_yticklabels(ordc); ax.set_xlabel("Rank acceptability (%)")
ax.set_xlim(0, 112); ax.set_xticks([0, 25, 50, 75, 100])
ax.legend(ncol=6, fontsize=7, bbox_to_anchor=(0.45, -0.22), loc="upper center", handlelength=1); ax.text(-0.09, 1.04, "(c)", transform=ax.transAxes, fontsize=9, fontweight="bold", color=INK)
fig.tight_layout(); save_elsevier(fig, "Fig6_fusion_topsis", 190); plt.close(fig)

# Fig. S1 – TIC overlay (optional; only if Agilent tic_back.csv files are supplied)
print("\n=== KEY RESULTS ===")
print(S.assign(v=lambda d: d.apply(lambda r: f"{r['mean']:.3g}±{r.sd:.2g}{r.letter}", axis=1)).pivot(index="var", columns="cond", values="v").loc[VARS, CONDS].to_string())
print(A.round(4).to_string()); print(T.round(3).to_string()); print(pd.DataFrame(dE).round(2).to_string())
print(anc.params.round(3).to_string(), "\nroaster p =", round(anc.pvalues["C(roaster)[T.MW]"], 5), "| interaction p =", round(anc_int.pvalues["RL:C(roaster)[T.MW]"], 3))
print(div.round(3).to_string()); print(pd.DataFrame(overlap).round(2).to_string()); print(TSVI.round(3).to_string())
print("MFA explained:", expl[:3].round(1)); print(scores.round(2).to_string()); print(corr.round(2).to_string())
print("entropy weights:", dict(zip(crit.columns, w_ent.round(3)))); print(topsis_tab.round(3).to_string())
print(f"\nAll outputs written to ./{OUT}/")
