set terminal svg size 1600,500 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure1_PCoA.svg'

set multiplot layout 1,3 title 'PCoA of Gut Microbiome Communities' font 'Arial Bold,16'

# Donor 1 panel
set title 'Donor 1' font 'Arial Bold,14'
set xlabel 'PC1 (51.2%% variance)' font 'Arial,11'
set ylabel 'PC2 (1.1%% variance)' font 'Arial,11'
set key outside right top font 'Arial,9' box
set grid

$D1_0g << EOD
-0.296435 0.000000
-0.298043 0.000000
-0.297756 0.000000
-0.293096 0.000000
-0.291022 0.000000
-0.286851 0.000000
-0.255498 0.000000
-0.234839 0.000000
-0.241967 0.000000
EOD

$D1_1_6g << EOD
-0.301442 0.000000
-0.299265 0.000000
-0.299953 0.000000
-0.298366 0.000000
-0.300551 0.000000
-0.300940 0.000000
-0.270874 0.000000
-0.275718 0.000000
-0.272599 0.000000
EOD

$D1_1g << EOD
-0.301470 0.000000
-0.299599 0.000000
-0.301592 0.000000
-0.298568 0.000000
-0.299539 0.000000
-0.299843 0.000000
-0.269535 0.000000
-0.255114 0.000000
-0.271817 0.000000
EOD

$D1_1g_s << EOD
-0.290094 0.000000
-0.300712 0.000000
-0.277638 0.000000
-0.271702 0.000000
-0.272151 0.000000
-0.224807 0.000000
-0.174153 0.000000
-0.195876 0.000000
-0.170253 0.000000
EOD

$D1_5g << EOD
-0.300431 0.000000
-0.299859 0.000000
-0.301554 0.000000
-0.301731 0.000000
-0.300010 0.000000
-0.300067 0.000000
-0.297870 0.000000
-0.300970 0.000000
-0.292855 0.000000
EOD

plot $D1_0g using 1:2 title '0g' with points pt 7 ps 1.8 lc rgb '#E64B35', \
     $D1_1_6g using 1:2 title '1/6g' with points pt 7 ps 1.8 lc rgb '#4DBBD5', \
     $D1_1g using 1:2 title '1g' with points pt 7 ps 1.8 lc rgb '#00A087', \
     $D1_1g_s using 1:2 title '1g-S' with points pt 7 ps 1.8 lc rgb '#3C5488', \
     $D1_5g using 1:2 title '5g' with points pt 7 ps 1.8 lc rgb '#F39B7F'

# Donor 2 panel
set title 'Donor 2' font 'Arial Bold,14'
set xlabel 'PC1 (51.2%% variance)' font 'Arial,11'
set ylabel 'PC2 (1.1%% variance)' font 'Arial,11'
set key outside right top font 'Arial,9' box
set grid

$D2_0g << EOD
0.357995 0.000000
0.360368 0.000000
0.361736 0.000000
0.345429 0.000000
0.380600 0.000000
0.358148 0.000000
0.429405 0.000000
0.416268 0.000000
0.419895 0.000000
EOD

$D2_1_6g << EOD
0.356087 0.000000
0.332552 0.000000
0.333181 0.000000
0.361730 0.000000
0.361561 0.000000
0.362243 0.000000
0.403356 0.000000
0.388293 0.000000
0.417329 0.000000
EOD

$D2_1g << EOD
0.333650 0.000000
0.363455 0.000000
0.368976 0.000000
0.362914 0.000000
0.365515 0.000000
0.331586 0.000000
0.414289 0.000000
0.396187 0.000000
0.386954 0.000000
EOD

$D2_1g_s << EOD
0.357846 0.000000
0.362908 0.000000
0.356623 0.000000
0.376903 0.000000
0.408228 0.000000
0.401219 0.000000
0.478283 0.000000
0.464059 0.000000
0.452810 0.000000
EOD

$D2_5g << EOD
0.336510 0.000000
0.345876 0.000000
0.355307 0.000000
0.361571 0.000000
0.363317 0.000000
0.342969 0.000000
0.406388 0.000000
0.357943 0.000000
0.356012 0.000000
EOD

plot $D2_0g using 1:2 title '0g' with points pt 7 ps 1.8 lc rgb '#E64B35', \
     $D2_1_6g using 1:2 title '1/6g' with points pt 7 ps 1.8 lc rgb '#4DBBD5', \
     $D2_1g using 1:2 title '1g' with points pt 7 ps 1.8 lc rgb '#00A087', \
     $D2_1g_s using 1:2 title '1g-S' with points pt 7 ps 1.8 lc rgb '#3C5488', \
     $D2_5g using 1:2 title '5g' with points pt 7 ps 1.8 lc rgb '#F39B7F'

# Donor 3 panel
set title 'Donor 3' font 'Arial Bold,14'
set xlabel 'PC1 (51.2%% variance)' font 'Arial,11'
set ylabel 'PC2 (1.1%% variance)' font 'Arial,11'
set key outside right top font 'Arial,9' box
set grid

$D3_0g << EOD
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

$D3_1_6g << EOD
-0.128451 0.000000
0.000000 0.000000
-0.118955 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

$D3_1g << EOD
-0.120839 0.000000
0.000000 0.000000
-0.118204 0.000000
-0.126226 0.000000
-0.125919 0.000000
-0.117454 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

$D3_1g_s << EOD
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.464422
EOD

$D3_5g << EOD
0.000000 0.000000
-0.124302 0.000000
-0.129151 0.000000
-0.122598 0.000000
0.000000 0.000000
-0.170900 0.000000
0.000000 0.000000
0.000000 0.000000
0.000000 0.000000
EOD

plot $D3_0g using 1:2 title '0g' with points pt 7 ps 1.8 lc rgb '#E64B35', \
     $D3_1_6g using 1:2 title '1/6g' with points pt 7 ps 1.8 lc rgb '#4DBBD5', \
     $D3_1g using 1:2 title '1g' with points pt 7 ps 1.8 lc rgb '#00A087', \
     $D3_1g_s using 1:2 title '1g-S' with points pt 7 ps 1.8 lc rgb '#3C5488', \
     $D3_5g using 1:2 title '5g' with points pt 7 ps 1.8 lc rgb '#F39B7F'

unset multiplot
