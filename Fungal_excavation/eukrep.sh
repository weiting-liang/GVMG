#MAG size > 3MB
find /path/Assay/02_rmhostCHM13v2_assembly_multiBin/ -type f -name "*_mul_metabat2.bin.*.fa" -size +3M > multi_bin.list
#EukRep
for i in `cat multi_bin.list `;do id=`basename $i .fa` && EukRep -i $i -o ./multi_bin/$id.euk.fa;done
#keep 3MB
find . -type f -size -3M -exec rm {} +
