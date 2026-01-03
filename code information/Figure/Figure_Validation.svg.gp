set terminal svg size 1600,1000 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure_Validation.svg'

set multiplot layout 2,3 title 'Model Validation: Predicted vs Actual (24h)' font 'Arial Bold,16'

set title 'Bacteroides' font 'Arial Bold,12'
set xlabel 'Actual (24h)' font 'Arial,10'
set ylabel 'Predicted (24h)' font 'Arial,10'
set key off
set grid
set xrange [0:0.7327]
set yrange [0:0.7327]

$data_0g_0 << EOD
0.169362 0.421051
0.106555 0.438830
0.183745 0.429116
0.085313 0.410558
0.079801 0.393732
0.067883 0.415421
0.359088 0.567390
0.319270 0.555656
0.311064 0.539999
EOD

$data_1g_0 << EOD
0.226038 0.476345
0.234775 0.482077
0.161177 0.500418
0.104579 0.444514
0.114321 0.436553
0.090335 0.438887
0.420107 0.592385
0.396320 0.629168
0.355050 0.587033
EOD

$data_5g_0 << EOD
0.238238 0.539957
0.219866 0.498445
0.196630 0.500388
0.124112 0.449505
0.136845 0.445621
0.089471 0.457310
0.531169 0.635564
0.519468 0.666056
0.475728 0.657506
EOD

plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \
     $data_0g_0 using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \
     $data_1g_0 using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \
     $data_5g_0 using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'

set title 'Staphylococcus' font 'Arial Bold,12'
set xlabel 'Actual (24h)' font 'Arial,10'
set ylabel 'Predicted (24h)' font 'Arial,10'
set key off
set grid
set xrange [0:0.3255]
set yrange [0:0.3255]

$data_0g_74 << EOD
0.196202 0.191710
0.259546 0.170698
0.235794 0.182576
0.226223 0.121656
0.218245 0.120456
0.149605 0.129540
0.100313 0.069166
0.072643 0.073025
0.059465 0.067875
EOD

$data_1g_74 << EOD
0.199589 0.126112
0.208327 0.098780
0.220932 0.129087
0.145572 0.091656
0.169530 0.082238
0.162184 0.091914
0.074394 0.056463
0.093383 0.045699
0.072555 0.036492
EOD

$data_5g_74 << EOD
0.175962 0.092830
0.172856 0.101217
0.205014 0.109529
0.125296 0.084115
0.159142 0.108697
0.144786 0.062863
0.067077 0.056593
0.062775 0.041949
0.045617 0.040254
EOD

plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \
     $data_0g_74 using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \
     $data_1g_74 using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \
     $data_5g_74 using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'

set title 'Pseudomonas' font 'Arial Bold,12'
set xlabel 'Actual (24h)' font 'Arial,10'
set ylabel 'Predicted (24h)' font 'Arial,10'
set key off
set grid
set xrange [0:0.3554]
set yrange [0:0.3554]

$data_0g_76 << EOD
0.030225 0.017211
0.027507 0.029826
0.020494 0.018864
0.000149 0.000000
0.000000 0.000000
0.000000 0.000000
0.251305 0.135068
0.237442 0.128342
0.264787 0.142132
EOD

$data_1g_76 << EOD
0.018396 0.029625
0.013317 0.016157
0.018086 0.036192
0.000000 0.000000
0.000205 0.000000
0.000000 0.000000
0.202671 0.139334
0.191282 0.111254
0.233530 0.143551
EOD

$data_5g_76 << EOD
0.011369 0.010508
0.018630 0.015061
0.014653 0.011202
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.148211 0.100435
0.152434 0.087437
0.186372 0.072982
EOD

plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \
     $data_0g_76 using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \
     $data_1g_76 using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \
     $data_5g_76 using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'

set title 'Prevotella' font 'Arial Bold,12'
set xlabel 'Actual (24h)' font 'Arial,10'
set ylabel 'Predicted (24h)' font 'Arial,10'
set key off
set grid
set xrange [0:0.3845]
set yrange [0:0.3845]

$data_0g_1 << EOD
0.001739 0.000000
0.004959 0.000000
0.001126 0.000000
0.210360 0.117267
0.209774 0.119843
0.240398 0.119260
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

$data_1g_1 << EOD
0.002447 0.000000
0.001646 0.000000
0.004483 0.000000
0.269122 0.129104
0.246748 0.150699
0.280355 0.137545
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

$data_5g_1 << EOD
0.003825 0.002967
0.002448 0.004803
0.003053 0.003432
0.289199 0.166176
0.265956 0.186694
0.349514 0.195992
0.000078 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \
     $data_0g_1 using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \
     $data_1g_1 using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \
     $data_5g_1 using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'

set title 'Escherichia_Shigella' font 'Arial Bold,12'
set xlabel 'Actual (24h)' font 'Arial,10'
set ylabel 'Predicted (24h)' font 'Arial,10'
set key off
set grid
set xrange [0:0.1645]
set yrange [0:0.1645]

$data_0g_48 << EOD
0.144919 0.047003
0.149589 0.058497
0.116092 0.067171
0.022604 0.014283
0.021158 0.026298
0.025981 0.027461
0.051401 0.030747
0.062390 0.029100
0.078970 0.039317
EOD

$data_1g_48 << EOD
0.091742 0.036944
0.082336 0.046210
0.098811 0.039543
0.013173 0.021226
0.013007 0.017912
0.015670 0.009422
0.027208 0.022581
0.038948 0.029542
0.042701 0.023520
EOD

$data_5g_48 << EOD
0.090769 0.023708
0.083316 0.037318
0.053808 0.033151
0.012065 0.018316
0.021600 0.008733
0.011105 0.004660
0.020871 0.016063
0.012877 0.018934
0.028704 0.018327
EOD

plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \
     $data_0g_48 using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \
     $data_1g_48 using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \
     $data_5g_48 using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'

set title 'Blautia' font 'Arial Bold,12'
set xlabel 'Actual (24h)' font 'Arial,10'
set ylabel 'Predicted (24h)' font 'Arial,10'
set key off
set grid
set xrange [0:0.1064]
set yrange [0:0.1064]

$data_0g_4 << EOD
0.052881 0.030369
0.064787 0.033331
0.070525 0.044151
0.025732 0.018539
0.026140 0.022793
0.033783 0.017292
0.030560 0.035191
0.044363 0.024214
0.032812 0.027511
EOD

$data_1g_4 << EOD
0.065056 0.034385
0.049603 0.026406
0.074761 0.028933
0.055105 0.020009
0.027930 0.020724
0.024190 0.025792
0.041744 0.024183
0.031395 0.014324
0.054039 0.020980
EOD

$data_5g_4 << EOD
0.070034 0.042976
0.096718 0.063613
0.081607 0.037615
0.057427 0.028482
0.032574 0.023133
0.028940 0.034739
0.043230 0.029143
0.030626 0.030558
0.029151 0.029704
EOD

plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \
     $data_0g_4 using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \
     $data_1g_4 using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \
     $data_5g_4 using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'

unset multiplot
