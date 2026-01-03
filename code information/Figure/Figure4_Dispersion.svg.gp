set terminal svg size 1600,500 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure4_Dispersion.svg'

set multiplot layout 1,3 title 'Beta Dispersion by Gravity Condition' font 'Arial Bold,16'

set title 'Donor 1' font 'Arial Bold,14'
set ylabel 'Distance to Centroid' font 'Arial,11'
set style fill solid 0.6
set boxwidth 0.6
set grid y

$box_D1_0g << EOD
1 0.331040
1 0.332648
1 0.332361
1 0.327701
1 0.325627
1 0.321456
1 0.290103
1 0.269444
1 0.276572
EOD

$box_D1_1_6g << EOD
2 0.318079
2 0.315903
2 0.316590
2 0.315004
2 0.317189
2 0.317577
2 0.287512
2 0.292356
2 0.289236
EOD

$box_D1_1g << EOD
3 0.305833
3 0.303963
3 0.305955
3 0.302931
3 0.303902
3 0.304207
3 0.273899
3 0.259478
3 0.276180
EOD

$box_D1_1g_s << EOD
4 0.345393
4 0.355998
4 0.332952
4 0.327025
4 0.327473
4 0.280206
4 0.229668
4 0.251335
4 0.225779
EOD

$box_D1_5g << EOD
5 0.299824
5 0.299252
5 0.300946
5 0.301124
5 0.299403
5 0.299459
5 0.297263
5 0.300363
5 0.292247
EOD

set xtics ('0g' 1, '1/6g' 2, '1g' 3, '1g-S' 4, '5g' 5) font 'Arial,9'

plot $box_D1_0g using 1:2 notitle with boxplot lc rgb '#E64B35', \
     $box_D1_1_6g using 1:2 notitle with boxplot lc rgb '#4DBBD5', \
     $box_D1_1g using 1:2 notitle with boxplot lc rgb '#00A087', \
     $box_D1_1g_s using 1:2 notitle with boxplot lc rgb '#3C5488', \
     $box_D1_5g using 1:2 notitle with boxplot lc rgb '#F39B7F'

set title 'Donor 2' font 'Arial Bold,14'
set ylabel 'Distance to Centroid' font 'Arial,11'
set style fill solid 0.6
set boxwidth 0.6
set grid y

$box_D2_0g << EOD
1 0.323390
1 0.325763
1 0.327131
1 0.310824
1 0.345995
1 0.323543
1 0.394800
1 0.381663
1 0.385290
EOD

$box_D2_1_6g << EOD
2 0.339449
2 0.315914
2 0.316544
2 0.345092
2 0.344923
2 0.345605
2 0.386718
2 0.371656
2 0.400691
EOD

$box_D2_1g << EOD
3 0.329287
3 0.359092
3 0.364613
3 0.358551
3 0.361152
3 0.327223
3 0.409926
3 0.391824
3 0.382591
EOD

$box_D2_1g_s << EOD
4 0.303464
4 0.308518
4 0.302243
4 0.322492
4 0.353776
4 0.346776
4 0.423762
4 0.409550
4 0.398311
EOD

$box_D2_5g << EOD
5 0.337117
5 0.346483
5 0.355915
5 0.362179
5 0.363925
5 0.343576
5 0.406995
5 0.358551
5 0.356620
EOD

set xtics ('0g' 1, '1/6g' 2, '1g' 3, '1g-S' 4, '5g' 5) font 'Arial,9'

plot $box_D2_0g using 1:2 notitle with boxplot lc rgb '#E64B35', \
     $box_D2_1_6g using 1:2 notitle with boxplot lc rgb '#4DBBD5', \
     $box_D2_1g using 1:2 notitle with boxplot lc rgb '#00A087', \
     $box_D2_1g_s using 1:2 notitle with boxplot lc rgb '#3C5488', \
     $box_D2_5g using 1:2 notitle with boxplot lc rgb '#F39B7F'

set title 'Donor 3' font 'Arial Bold,14'
set ylabel 'Distance to Centroid' font 'Arial,11'
set style fill solid 0.6
set boxwidth 0.6
set grid y

$box_D3_0g << EOD
1 0.034605
1 0.034605
1 0.034605
1 0.034605
1 0.034605
1 0.034605
1 0.034605
1 0.034605
1 0.034605
EOD

$box_D3_1_6g << EOD
2 0.145089
2 0.016638
2 0.135593
2 0.016638
2 0.016638
2 0.016638
2 0.016638
2 0.016638
2 0.016638
EOD

$box_D3_1g << EOD
3 0.125202
3 0.004363
3 0.122567
3 0.130589
3 0.130282
3 0.121817
3 0.004363
3 0.004363
3 0.004363
EOD

$box_D3_1g_s << EOD
4 0.057503
4 0.057503
4 0.057503
4 0.057503
4 0.057503
4 0.057503
4 0.057503
4 0.057503
4 0.450575
EOD

$box_D3_5g << EOD
5 0.000608
5 0.123694
5 0.128543
5 0.121991
5 0.000608
5 0.170293
5 0.000608
5 0.000608
5 0.000608
EOD

set xtics ('0g' 1, '1/6g' 2, '1g' 3, '1g-S' 4, '5g' 5) font 'Arial,9'

plot $box_D3_0g using 1:2 notitle with boxplot lc rgb '#E64B35', \
     $box_D3_1_6g using 1:2 notitle with boxplot lc rgb '#4DBBD5', \
     $box_D3_1g using 1:2 notitle with boxplot lc rgb '#00A087', \
     $box_D3_1g_s using 1:2 notitle with boxplot lc rgb '#3C5488', \
     $box_D3_5g using 1:2 notitle with boxplot lc rgb '#F39B7F'

unset multiplot
