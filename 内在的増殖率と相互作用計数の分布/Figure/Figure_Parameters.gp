set terminal svg size 1400,600 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure_Parameters.svg'

set multiplot layout 1,2 title 'gLV Parameter Distributions' font 'Arial Bold,16'

# Panel 1: Intrinsic Growth Rates
set title 'Intrinsic Growth Rates (r_i)' font 'Arial Bold,12'
set xlabel 'Growth Rate' font 'Arial,10'
set ylabel 'Frequency' font 'Arial,10'
set grid
set style fill solid 0.7
set boxwidth 0.008196
set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'gray' lw 1 dt 2
set label 1 'Mean: -0.0338\nSD: 0.0601\nN: 120' at graph 0.95,0.95 right font 'Arial,9'

$data_r << EOD
-0.159393 1
-0.150286 2
-0.141180 1
-0.132074 2
-0.122968 3
-0.113862 2
-0.104756 7
-0.095649 7
-0.086543 6
-0.077437 7
-0.068331 3
-0.059225 5
-0.050119 3
-0.041013 4
-0.031906 8
-0.022800 11
-0.013694 7
-0.004588 0
0.004518 17
0.013624 4
0.022731 1
0.031837 1
0.040943 2
0.050049 4
0.059155 3
0.068261 4
0.077367 3
0.086474 1
0.095580 0
0.104686 1
EOD
plot $data_r using 1:2 with boxes lc rgb '#E64B35' notitle

# Panel 2: Interaction Coefficients
unset arrow
unset label
set title 'Interaction Coefficients (A_{ij})' font 'Arial Bold,12'
set xlabel 'Coefficient Value' font 'Arial,10'
set ylabel 'Frequency' font 'Arial,10'
set boxwidth 0.009164
set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'red' lw 2 dt 2
set label 1 'Mean: -0.0003\nSD: 0.0151\nNear zero: 53.0%%' at graph 0.95,0.95 right font 'Arial,9'

$data_a << EOD
-0.265066 1
-0.254883 1
-0.244701 1
-0.234518 0
-0.224336 0
-0.214153 1
-0.203971 0
-0.193788 3
-0.183606 0
-0.173423 0
-0.163241 4
-0.153058 3
-0.142876 6
-0.132693 4
-0.122511 7
-0.112328 4
-0.102146 7
-0.091963 9
-0.081781 11
-0.071599 28
-0.061416 27
-0.051234 54
-0.041051 64
-0.030869 122
-0.020686 284
-0.010504 817
-0.000321 11421
0.009861 946
0.020044 255
0.030226 132
0.040409 49
0.050591 46
0.060774 22
0.070956 20
0.081139 20
0.091321 7
0.101504 3
0.111686 6
0.121869 1
0.132051 5
0.142233 3
0.152416 2
0.162598 1
0.172781 0
0.182963 0
0.193146 0
0.203328 0
0.213511 2
0.223693 0
0.233876 1
EOD
plot $data_a using 1:2 with boxes lc rgb '#4DBBD5' notitle

unset multiplot
