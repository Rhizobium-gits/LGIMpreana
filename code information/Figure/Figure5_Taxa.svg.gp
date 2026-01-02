set terminal svg size 1000,600 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure5_Taxa.svg'

set title 'Top Taxa by Gravity Condition' font 'Arial Bold,14'
set style data histogram
set style histogram clustered gap 1
set style fill solid 0.8 border -1
set boxwidth 0.9
set key outside right top font 'Arial,9'
set ylabel 'Mean Relative Abundance' font 'Arial Bold,11'
set xtics rotate by -45 right font 'Arial,9'
set grid y

$data << EOD
Taxon 0g 1_6g 1g 1g_s 5g
"Bacteroides" 0.2328 0.2673 0.2658 0.1706 0.3000
"Prevotella" 0.1073 0.1093 0.1130 0.0782 0.1273
"Staphylococcus" 0.1099 0.1003 0.0847 0.1406 0.0702
"Blautia" 0.0619 0.0587 0.0616 0.0601 0.0639
"Pseudomonas" 0.0581 0.0515 0.0462 0.0695 0.0333
"Escherichia_..." 0.0472 0.0322 0.0289 0.0538 0.0216
"Faecalibacte..." 0.0238 0.0291 0.0347 0.0188 0.0418
"Acinetobacter" 0.0305 0.0253 0.0264 0.0382 0.0182
"Streptococcus" 0.0254 0.0185 0.0169 0.0465 0.0135
"Veillonella" 0.0209 0.0180 0.0212 0.0238 0.0201
EOD

plot $data using 2:xtic(1) title 'Microgravity' lc rgb '#E64B35', \
     $data using 3:xtic(1) title 'Lunar (1/6g)' lc rgb '#4DBBD5', \
     $data using 4:xtic(1) title 'Earth (1g)' lc rgb '#00A087', \
     $data using 5:xtic(1) title 'Static' lc rgb '#3C5488', \
     $data using 6:xtic(1) title 'Hypergravity' lc rgb '#F39B7F'
