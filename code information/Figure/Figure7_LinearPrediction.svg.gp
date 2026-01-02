set terminal svg size 1200,800 enhanced font 'Arial,12' background '#FFFFFF'
set output '/tmp/microbiome_results/Figure7_LinearPrediction.svg'

set multiplot layout 2,3 title 'Linear Prediction: 48h Composition' font 'Arial Bold,14'

set title 'Microgravity' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.8
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y

$data_0g << EOD
0 0.3383 0.1500 "Staphyloco.."
1 0.1793 0.2193 "Pseudomonas"
2 0.1477 0.0954 "Escherichi.."
3 0.0945 0.0873 "Acinetobac.."
4 0.0784 0.0724 "Bacteroides"
EOD
plot $data_0g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#E64B35' notitle

set title 'Lunar (1/6g)' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.8
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y

$data_1_6g << EOD
0 0.3248 0.1486 "Staphyloco.."
1 0.1678 0.2091 "Pseudomonas"
2 0.1403 0.1053 "Bacteroides"
3 0.1085 0.0708 "Escherichi.."
4 0.0716 0.0809 "Acinetobac.."
EOD
plot $data_1_6g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#4DBBD5' notitle

set title 'Earth (1g)' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.8
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y

$data_1g << EOD
0 0.3176 0.1166 "Staphyloco.."
1 0.1793 0.0984 "Bacteroides"
2 0.1578 0.2052 "Pseudomonas"
3 0.0964 0.0707 "Escherichi.."
4 0.0916 0.0987 "Acinetobac.."
EOD
plot $data_1g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#00A087' notitle

set title 'Static' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.8
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y

$data_1g_s << EOD
0 0.4022 0.1415 "Staphyloco.."
1 0.1963 0.2529 "Pseudomonas"
2 0.1488 0.0720 "Escherichi.."
3 0.1300 0.0749 "Streptococ.."
4 0.0928 0.1045 "Acinetobac.."
EOD
plot $data_1g_s using 1:2:3:xtic(4) with boxerrorbars lc rgb '#3C5488' notitle

set title 'Hypergravity' font 'Arial Bold,12'
set style fill solid 0.7
set boxwidth 0.8
set ylabel 'Predicted Abundance'
set xtics rotate by -45 right font 'Arial,9'
set grid y

$data_5g << EOD
0 0.2844 0.1216 "Staphyloco.."
1 0.2493 0.1641 "Bacteroides"
2 0.1255 0.1584 "Pseudomonas"
3 0.0775 0.0635 "Escherichi.."
4 0.0654 0.0643 "Acinetobac.."
EOD
plot $data_5g using 1:2:3:xtic(4) with boxerrorbars lc rgb '#F39B7F' notitle

unset multiplot
