export PATH=`pwd`/../graphlan/:$PATH

species=data/my.species_2022.txt
annot_f=annot/myannotR_2022.txt
out_prefix=species_2022

graphlan_annotate.py --annot $annot_f $species xml/$out_prefix.xml
graphlan.py xml/$out_prefix.xml img/$out_prefix.png --dpi 300 --size 3.5
graphlan.py xml/$out_prefix.xml img/$out_prefix.svg --dpi 300 --size 3.5
