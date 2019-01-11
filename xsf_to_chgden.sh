#!/bin/bash
chgden() {
   #cp $1 $2
   TMP=$(head -7 $1 | tail -1 | awk '{$1=$1;print}' | awk '{ print $1 }')
   TMP2=$(head -7 $1 | tail -1)
   sed -e "s@PRIMVEC@1@g" -e "s@ PRIMCOORD@$TMP@g" -e "s@$TMP2@Cartesian@g" $1 > temp1
   head -7 temp1 > $1.chg
   tail -n +8 temp1 | head -$TMP | awk '{ print $2, $3, $4 }' >> $1.chg
   echo "" >> $1.chg
   echo "136   136   136" >> $1.chg
   TMP3=$(echo "$(($TMP+16))")
   tail -n +$TMP3 $1>> $1.chg
   rm temp1
}
