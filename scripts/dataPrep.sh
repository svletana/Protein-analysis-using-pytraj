#!/bin/bash

pdbs='2m4l_sh 1j4m_sh' # folder names
tem='250' # temperature
OLD='VAR'

for NEW in ${pdbs}
do
 sed -i "s/${OLD}/${NEW}/g" afterProd.in
 mv afterProd.in ${NEW}/temp${tem}k
 cd ${NEW}/temp${tem}k
 mkdir input
 cpptraj -i afterProd.in
 mv nowat* input/
 mv afterProd.in ../../
 cd ../../
 OLD=${NEW}
done

sed -i "s/${OLD}/VAR/g" afterProd.in
