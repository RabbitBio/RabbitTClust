#!/bin/bash 
#set -x 

cd ref/
resName="refList"
if [ -f "$resName" ]; then
	echo "the file exist, remove"
	rm $resName
fi

for dir in archaea bacteria fungi viral plant protozoa human vertebrate_mammalian vertebrate_other
do
	#echo $dir
	ls $dir/*.fna.gz >$dir.gz.list
	cat $dir.gz.list | while read line
	do 
		#echo $line
		gunzip $line
	done
	rm $dir.gz.list

	ls $dir/*.fna > $dir.list
	cat  $dir.list | while read line
	do
		echo `pwd`/$line >>$resName
	done
	rm $dir.list
done

mv $resName ../
cd ../



