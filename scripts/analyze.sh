#!/bin/bash

# ==============================================================================
# == Launch analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# -- Variables
# ------------------------------------------------------------------------------

home=/local1/work/ginkgodev
dir=${home}/uploads/${1}
source ${dir}/config
distMet=$distMeth
touch $dir/index.html

inFile=list
statFile=status.xml
genome=${home}/genomes/${chosen_genome}

if [ "$rmpseudoautosomal" == "1" ];
then
  genome=${genome}/pseudoautosomal
else
  genome=${genome}/original
  #genome=${genome}/win0.2/p0.6
  #genome=${genome}/win0.2/p0.5635
fi

# ------------------------------------------------------------------------------
# -- Error Check and Reformat User Files
# ------------------------------------------------------------------------------

if [ "$f" == "0" ]; then
  touch ${dir}/ploidyDummy.txt
  facs=ploidyDummy.txt
else 
  # In case upload file with \r instead of \n (Mac, Windows)
  tr '\r' '\n' < ${dir}/${facs} > ${dir}/quickTemp
  mv ${dir}/quickTemp ${dir}/${facs}
  # 
  sed "s/.bed//g" ${dir}/${facs} | sort -k1,1 | awk '{print $1"\t"$2}' > ${dir}/quickTemp 
  mv ${dir}/quickTemp ${dir}/${facs}
fi

cp uploads/qnXwCIkHowiwsN725giF/clust.xml uploads/_t10breast_navin/

# ------------------------------------------------------------------------------
# -- Map Reads & Prepare Samples For Processing
# ------------------------------------------------------------------------------

timeA=`date +%s.%N`

total=`wc -l < ${dir}/${inFile}`

export init="1";

if [ "$init" == "1" ];
then

  if [ -z "$clustThresh" ];
  then
    echo "Cleaning directory"
    # Clean directory
    rm -f ${dir}/*_mapped ${dir}/*.jpeg ${dir}/*.newick ${dir}/*.xml ${dir}/*.cnv ${dir}/Seg* ${dir}/results.txt
    rm -r -f ${dir}/CNVprofiles
  fi

  # Map user bed files to appropriate bins
  cnt=0
  while read file;
  do
    ${home}/scripts/status ${dir}/${statFile} 1 $file $cnt $total

    # Unzip gzip files if necessary
    if [[ "${file}" =~ \.gz$ ]];
    then
      firstLineChr=$(zcat ${dir}/${file} | head -n 1 | cut -f1 | grep "chr")
      if [[ "${firstLineChr}" == "" ]];
      then
        awk '{print "chr"$0}' <(zcat ${dir}/${file}) > ${dir}/${file}_tmp
        mv ${dir}/${file}_tmp ${dir}/${file/.gz/}
        gzip -f ${dir}/${file/.gz/}
      fi

      if [ -f ${dir}/${file}.cells ];
      then
        echo "${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` <(zcat -cd ${dir}/${file}) `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mapped_ `wc -l < ${dir}/${file}.cells` ${dir}/${file}.cells"
        ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` <(zcat -cd ${dir}/${file}) `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mapped_ `wc -l < ${dir}/${file}.cells` ${dir}/${file}.cells
      else
        echo "${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` <(zcat -cd ${dir}/${file}) `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mapped 1"
        ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` <(zcat -cd ${dir}/${file}) `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mapped 1
      fi

      if [ "$improvebounds" == "1" ]; then
        if [ -f ${dir}/${file}.cells ];
        then
          ${home}/scripts/binUnsorted ${genome}/${binMethFine} `wc -l < ${genome}/${binMethFine}` <(zcat -cd ${dir}/${file}) `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mappedfine_ `wc -l < ${dir}/${file}.cells` ${dir}/${file}.cells
        else
          ${home}/scripts/binUnsorted ${genome}/${binMethFine} `wc -l < ${genome}/${binMethFine}` <(zcat -cd ${dir}/${file}) `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mappedfine 1
        fi
      fi

    # 
    else
      firstLineChr=$( head -n 1 ${dir}/${file} | cut -f1 | grep "chr")
      if [[ "${firstLineChr}" == "" ]];
      then
        awk '{print "chr"$0}' ${dir}/${file} > ${dir}/${file}_tmp
        mv ${dir}/${file}_tmp ${dir}/${file}
      fi

      if [ -f ${dir}/${file}.cells ];
      then
        ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` ${dir}/${file} `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mapped_ `wc -l < ${dir}/${file}.cells` ${dir}/${file}.cells
      else
        ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` ${dir}/${file} `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mapped 1
      fi

      if [ "$improvebounds" == "1" ]; then
        if [ -f ${dir}/${file}.cells ];
        then
          ${home}/scripts/binUnsorted ${genome}/${binMethFine} `wc -l < ${genome}/${binMethFine}` ${dir}/${file} `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mappedfine_ `wc -l < ${dir}/${file}.cells` ${dir}/${file}.cells
        else
          ${home}/scripts/binUnsorted ${genome}/${binMethFine} `wc -l < ${genome}/${binMethFine}` ${dir}/${file} `echo ${file} | awk -F ".bed" '{print $1}'` ${dir}/${file}_mappedfine 1
        fi
      fi
      gzip ${dir}/${file}
    fi

    cnt=$(($cnt+1))
  done < ${dir}/${inFile}

  # Concatenate binned reads to central file  
  if [ "$improvebounds" == "1" ]; then
    paste ${dir}/*_mappedfine* > ${dir}/dataFine
    rm -f ${dir}/*_mappedfine* ${dir}/*_binned*
  fi
  paste ${dir}/*_mapped* > ${dir}/data
  #rm -f ${dir}/*_mapped* ${dir}/*_binned
fi

# ------------------------------------------------------------------------------
# -- Map User Provided Reference/Segmentation Sample
# ------------------------------------------------------------------------------

if [ "$segMeth" == "2" ]; then
  if [ -f ${dir}/${ref}.cells ];
  then
    ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` ${dir}/${ref} Reference ${dir}/${ref}_mapped_ `wc -l < ${dir}/${ref}.cells` ${dir}/${ref}.cells
  else
    ${home}/scripts/binUnsorted ${genome}/${binMeth} `wc -l < ${genome}/${binMeth}` ${dir}/${ref} Reference ${dir}/${ref}_mapped 1
  fi
else
    ref=refDummy.bed
    touch ${dir}/${ref}_mapped
fi

timeB=`date +%s.%N`
runtime=$(echo "$timeB-$timeA" | bc)
rm ${dir}/timing.txt
echo -e $runtime >> ${dir}/timing.txt
#echo -e  "Time to map reads: " $runtime >> ${dir}/timing.txt

# ------------------------------------------------------------------------------
# -- Run Mapped Data Through Primary Pipeline
# ------------------------------------------------------------------------------

if [ "$process" == "1" ]; then
  if [ "$improvebounds" == "1" ]; then
    echo "Launching process.R $genome $dir $statFile data $segMeth $binMeth $clustMeth $distMet $color ${ref}_mapped $f $facs $sex $rmbadbins $maxploidy $minbinwidth $improvebounds $binMethFine dataFine $clustThresh"
    ${home}/scripts/process.R $genome $dir $statFile data $segMeth $binMeth $clustMeth $distMet $color ${ref}_mapped $f $facs $sex $rmbadbins $maxploidy $minbinwidth $improvebounds $binMethFine dataFine $clustThresh
  else
    echo "Launching process.R $genome $dir $statFile data $segMeth $binMeth $clustMeth $distMet $color ${ref}_mapped $f $facs $sex $rmbadbins $maxploidy $minbinwidth $improvebounds $clustThresh"
    ${home}/scripts/process.R $genome $dir $statFile data $segMeth $binMeth $clustMeth $distMet $color ${ref}_mapped $f $facs $sex $rmbadbins $maxploidy $minbinwidth $improvebounds $clustThresh
  fi
fi

timeC=`date +%s.%N`
runtime=$(echo "$timeC-$timeB" | bc)
echo -e $runtime >> ${dir}/timing.txt
#echo -e  "Time to segment:   " $runtime >> ${dir}/timing.txt

# ------------------------------------------------------------------------------
# -- Recreate Clusters/Heat Maps (With New Parameters)
# ------------------------------------------------------------------------------

if [ "$fix" == "1" ]; then
  echo "Launching reclust.R $genome $dir $statFile $binMeth $clustMeth $distMet $f $facs $sex"
  ${home}/scripts/reclust.R $genome $dir $statFile $binMeth $clustMeth $distMet $f $facs $sex
fi

timeD=`date +%s.%N`
runtime=$(echo "$timeD-$timeC" | bc)
echo -e $runtime >> ${dir}/timing.txt
#echo -e  "Time to recluster: " $runtime >> ${dir}/timing.txt

# ------------------------------------------------------------------------------
# -- Create CNV profiles
# ------------------------------------------------------------------------------

nbCols=$(awk '{ print NF; exit; }' $dir/SegCopy)
for (( i=1; i<=$nbCols; i++ ));
do
  if [ -z "$clustThresh" ];
  then
    currCell=$(cut -f$i $dir/SegCopy | head -n 1 | tr -d '"')
    if [ "$currCell" == "" ]; then
      continue;
    fi
    cut -f$i $dir/SegCopy | tail -n+2 | awk '{if(NR==1) print "1,"$1; else print NR","prev"\n"NR","$1;prev=$1; }' > $dir/$currCell.cnv
  else
    currCell=$(cut -f$i $dir/SegCopyClust | head -n 1 | tr -d '"')
    if [ "$currCell" == "" ]; then
      continue;
    fi
    cut -f$i $dir/SegCopyClust | tail -n+2 | awk '{if(NR==1) print "1,"$1; else print NR","prev"\n"NR","$1;prev=$1; }' > $dir/$currCell.cnv
  fi
done

timeE=`date +%s.%N`
runtime=$(echo "$timeE-$timeD" | bc)
echo -e $runtime >> ${dir}/timing.txt
#echo -e  "Time to profile:   " $runtime >> ${dir}/timing.txt

# ------------------------------------------------------------------------------
# -- Call CNVs
# ------------------------------------------------------------------------------

echo "Launching ${home}/scripts/CNVcaller ${dir}/SegCopy ${dir}/CNV1 ${dir}/CNV2"
${home}/scripts/CNVcaller ${dir}/SegCopy ${dir}/CNV1 ${dir}/CNV2

mkdir -p ${dir}/CNVprofiles
cp ${dir}/*_CN.jpeg ${dir}/CNVprofiles/
tar -zcvf ${dir}/CNVprofiles.tar.gz ${dir}/CNVprofiles

timeF=`date +%s.%N`
runtime=$(echo "$timeF-$timeE" | bc)
echo -e $runtime >> ${dir}/timing.txt
echo -e "Done" >> ${dir}/timing.txt
#echo -e  "Time to call CNVs: " $runtime >> ${dir}/timing.txt

# ------------------------------------------------------------------------------
# -- Email notification of completion
# ------------------------------------------------------------------------------

if [ "$email" != "" ]; then
	echo -e "Your analysis on Ginkgo is complete! Check out your results at $permalink" | mail -s "Your Analysis Results" $email -- -F "Ginkgo"
fi
