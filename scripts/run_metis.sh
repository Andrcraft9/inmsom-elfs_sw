#!/bin/bash

AREA=$1
BLOCKS=$2
PROCS=$3

python mesh_metis_graph.py ${AREA}/${BLOCKS}b_${PROCS}p_h.txt > ${AREA}/graph${BLOCKS}b
./gpmetis ${AREA}/graph${BLOCKS}b ${PROCS}
head -n 1 ${AREA}/${BLOCKS}b_${PROCS}p_h.txt >  ${AREA}/${BLOCKS}b_${PROCS}p_metis.txt
python mesh_metis_graph.py ${AREA}/${BLOCKS}b_${PROCS}p_h.txt ${AREA}/graph${BLOCKS}b.part.${PROCS} >> ${AREA}/${BLOCKS}b_${PROCS}p_metis.txt
echo "Created: "
echo  ${AREA}/${BLOCKS}b_${PROCS}p_metis.txt

