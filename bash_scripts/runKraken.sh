#!/bin/bash

~/tools/kraken/kraken --db ~/tools/kraken/minikraken_20141208/ $1 --fastq-input --quick --min-hits=2 > temp 
~/tools/kraken/kraken-translate --db ~/tools/kraken/minikraken_20141208/ temp > "$1-kr-sp"
awk 'FS=";" {print $8	$9}' "$1-kr-sp" | sort | uniq -c | sort -nr > "$1-kr-ov"
rm temp
