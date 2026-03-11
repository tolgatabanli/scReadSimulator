java -jar target/scReadSimulator-1.0-SNAPSHOT.jar \
  -readcounts /path/to/readcounts.tsv \
  -fasta /path/to/genome.fa \
  -fidx /path/to/genome.fa.fai \
  -gtf /path/to/annotations.gtf \
  -od /path/to/output_dir \
  -length 75 \
  -tailLength 200 \
  -mutationrate 1.0