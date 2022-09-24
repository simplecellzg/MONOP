#!/bin/bash
#change ms data to lammps data
aa=`ls *.car`
echo $aa
name=`basename $aa .car`
echo $name
msi2lmp.exe $name -class I -frc cvff -i
#delete unuseful rows
sed -i '4,8d;10,13d' $name.data
bb=`sed -n -e "/Pair Coeffs/=" $name.data`
cc=`sed -n -e "/Atoms/=" $name.data`
declare -i bb=$bb
declare -i cc=$cc-1
echo $bb
echo $cc
sed -i $bb','$cc'd' $name.data
dd=`sed -n -e "/Bonds/=" $name.data`
declare -i dd=$dd
echo $dd
sed -i $dd',$d' $name.data
#copy pbs file and rename
pbs=`ls *lammps_cpu.pbs`
file=$pbs
newf=$name.pbs
#change homework name; in. file name
sed -e "s/template/$name/g" $file > $newf
sed -i "s/name/reaxc.$name/g" $newf
