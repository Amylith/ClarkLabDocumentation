
#remove sections of promoters that overlap genes using bedtools
#NOTE: if you do this, you must use the "ExtractGenesFrommaf.R" to get promoters from maf files
#bedtools subtract -a PromoterCoords2000.txt -b ../ExonCoords.txt > PromoterCoords2000nocoding.txt
