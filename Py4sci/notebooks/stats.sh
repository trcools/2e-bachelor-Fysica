#!/usr/bin/env bash

#This is based on a suggestion found here
# https://stackoverflow.com/questions/4589731/git-blame-commit-statistics

rm -rf stats.tmp
for i in *.ipynb;
	do echo $i;
	git blame  --line-porcelain $i |\
	grep -ae "^author "|\
	sed 's/author//'|\
	sort|\
	sed 's/JLeliaert/Jonathan Leliaert/'|\
	sed 's/jleliaer/Jonathan Leliaert/' |\
	sed 's/yxcheng/Yingching Cheng/' |\
	sed 's/wamissia/Warre Missiaen/' |\
	uniq -c|\
	sort -nr | tee -a stats.tmp;
	done

echo "\nSUM TOTAL"
sed -i 's/^[ \t]*//' stats.tmp
sed -i 's/ /\t/1' stats.tmp
awk -F"\t" '{a[($2)]+=($1)}END{for(x in a)print("   ",a[x],x)}' stats.tmp|sort -nr
rm -rf stats.tmp
