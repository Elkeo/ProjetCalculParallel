#!/bin/bash
#zsh SU_EF.sh n
NBPROC=$1
RED='\033[0;31m'
NC='\033[0m'
rm TP_PROG_EF/metis/time.dat
rm TP_PROG_EF/scotch/time.dat

cd "TP_METIS_SCOTCH/"
echo "${RED}Actuellement dans $PWD${NC}"

echo "${RED}Run de data2tec${NC}"
./data2tec

echo "${RED}Run de m2gmetis${NC}"
m2gmetis dualformetis.dat dualformetis.dat.dgraph

echo "${RED}Run de metisdual2scotchdual${NC}"
./metisdual2scotchdual dualformetis.dat.dgraph dualforscotch.grf

echo "${RED}Run de gtst${NC}"
gtst dualforscotch.grf

for i in {1..$NBPROC}
do
   if (($i == 1))
   then
      g++ createmeshfor1proc.cpp -o createmeshfor1proc
      ./createmeshfor1proc
   else
      echo "${RED}Run de mpmetis${NC}"
      mpmetis dualformetis.dat $i
   fi
   echo "${RED}Run de la décomposition scotch${NC}"
   echo cmplt $i | gmap dualforscotch.grf - dualforscotch.map
   echo ""
   echo "${RED}Actuellement dans $PWD${NC}"
   cd "../TP_PROG_EF"
   echo "${RED}Actuellement dans $PWD${NC}"
   echo "${RED}Création du metis/mesh_for_progc.data${NC}"
   cat ../TP_METIS_SCOTCH/meshprogc.data ../TP_METIS_SCOTCH/dualformetis.dat.epart.$i > metis/mesh_for_progc.data
   echo "${RED}Création du scotch/mesh_for_progc.data${NC}"
   ./fromscotch ../TP_METIS_SCOTCH/meshprogc.data ../TP_METIS_SCOTCH/dualforscotch.map scotch/mesh_for_progc.data

   ./Preprocess $i metis/mesh_for_progc.data

   mpicc -lm -o Fem FemPar.c
   mpirun -np $i ./Fem metis/

   echo "${RED}Suppression des fichiers Data00.In${NC}"
   for j in {0..$((i-1))}
   do
      rm Data0$j.In
   done

   ./Preprocess $i scotch/mesh_for_progc.data

   mpirun -np $i ./Fem scotch/

   echo "${RED}Suppression des fichiers Data00.In${NC}"
   for j in {0..$((i-1))}
   do
      rm Data0$j.In
   done

   cd ".."
   echo "${RED}Actuellement dans $PWD${NC}"
   echo "${RED}Suppression des fichiers inutiles...${NC}"
   rm TP_PROG_EF/metis/mesh_for_progc.data
   rm TP_PROG_EF/scotch/mesh_for_progc.data
   rm TP_METIS_SCOTCH/dualformetis.dat.epart.$i
   rm TP_METIS_SCOTCH/dualformetis.dat.npart.$i
   rm TP_METIS_SCOTCH/dualforscotch.map
   rm TP_PROG_EF/Sol0*
   cd "TP_METIS_SCOTCH/"
done

