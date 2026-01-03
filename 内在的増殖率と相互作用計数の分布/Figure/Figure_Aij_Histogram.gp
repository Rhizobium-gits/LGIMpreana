set terminal svg size 800,600 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure_Aij_Histogram.svg'

set title 'Distribution of Interaction Coefficients (A_{ij})' font 'Arial Bold,14'
set xlabel 'Coefficient Value' font 'Arial,12'
set ylabel 'Frequency' font 'Arial,12'
set grid
set style fill solid 0.7
set boxwidth 0.009164
set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'red' lw 2 dt 2
set label 1 'Mean: -0.0003\nSD: 0.0151\nNear zero: 53.0%%' at graph 0.95,0.95 right font 'Arial,10'

$data << EOD
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

plot $data using 1:2 with boxes lc rgb '#4DBBD5' notitle
