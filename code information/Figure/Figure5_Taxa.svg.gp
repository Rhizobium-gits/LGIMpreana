set terminal svg size 1200,600 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure5_Taxa.svg'

set title 'Top Taxa Abundance by Gravity (Mean ± SD)' font 'Arial Bold,14'
set style data histogram
set style histogram clustered gap 1 errorbars lw 1
set style fill solid 0.8 border -1
set boxwidth 0.9
set key outside right top font 'Arial,10' box
set ylabel 'Relative Abundance' font 'Arial Bold,11'
set xtics rotate by -45 right font 'Arial,10'
set grid y
set bars 2

$data << EOD
"Bacteroides" 0.2328 0.1517 0.2673 0.1645 0.2658 0.1519 0.3000 0.1824
"Prevotella" 0.1073 0.1587 0.1093 0.1563 0.1130 0.1638 0.1273 0.1826
"Staphyloco.." 0.1099 0.0699 0.1003 0.0669 0.0847 0.0609 0.0702 0.0565
"Blautia" 0.0619 0.0255 0.0587 0.0260 0.0616 0.0236 0.0639 0.0249
"Pseudomonas" 0.0581 0.0824 0.0515 0.0760 0.0462 0.0694 0.0333 0.0525
"Escherichi.." 0.0472 0.0393 0.0322 0.0290 0.0289 0.0258 0.0216 0.0223
"Faecalibac.." 0.0238 0.0221 0.0291 0.0286 0.0347 0.0325 0.0418 0.0348
"Acinetobac.." 0.0305 0.0352 0.0253 0.0313 0.0264 0.0339 0.0182 0.0231
EOD

plot $data using 2:3:xtic(1) title '0g' lc rgb '#E64B35', \
     $data using 4:5:xtic(1) title '1/6g' lc rgb '#4DBBD5', \
     $data using 6:7:xtic(1) title '1g' lc rgb '#00A087', \
     $data using 8:9:xtic(1) title '5g' lc rgb '#F39B7F'
