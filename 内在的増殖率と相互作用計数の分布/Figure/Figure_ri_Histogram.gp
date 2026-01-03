set terminal svg size 800,600 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure_ri_Histogram.svg'

set title 'Distribution of Intrinsic Growth Rates (r_i)' font 'Arial Bold,14'
set xlabel 'Growth Rate' font 'Arial,12'
set ylabel 'Frequency' font 'Arial,12'
set grid
set style fill solid 0.7
set boxwidth 0.004917
set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'red' lw 2 dt 2
set label 1 'Mean: -0.0338\nSD: 0.0601\nN: 120' at graph 0.95,0.95 right font 'Arial,10'

$data << EOD
-0.161214 1
-0.155750 0
-0.150286 2
-0.144823 0
-0.139359 1
-0.133895 0
-0.128432 2
-0.122968 2
-0.117504 3
-0.112041 0
-0.106577 6
-0.101113 2
-0.095649 3
-0.090186 5
-0.084722 4
-0.079258 4
-0.073795 4
-0.068331 2
-0.062867 3
-0.057404 2
-0.051940 2
-0.046476 1
-0.041013 4
-0.035549 3
-0.030085 5
-0.024621 6
-0.019158 6
-0.013694 4
-0.008230 2
-0.002767 0
0.002697 16
0.008161 2
0.013624 1
0.019088 3
0.024552 0
0.030015 1
0.035479 1
0.040943 1
0.046407 2
0.051870 2
0.057334 2
0.062798 2
0.068261 3
0.073725 1
0.079189 2
0.084652 0
0.090116 1
0.095580 0
0.101043 0
0.106507 1
EOD

plot $data using 1:2 with boxes lc rgb '#E64B35' notitle
