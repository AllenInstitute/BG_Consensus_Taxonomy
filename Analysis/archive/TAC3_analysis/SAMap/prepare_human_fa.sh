
## From gtf extract protein coding transcripts and their CDS
zcat genes/genes.gtf.gz | awk '$3=="CDS"' > genes/genes.CDS.gtf
awk '$1 ~ /^chr([0-9]+|X|Y|M)$/' genes/genes.CDS.gtf > genes/genes.CDS_only.gtf ## Remove non-standard chromosomes/contigs

## Extract genome fasta for gffread
gffread -g fasta/genome.fa -x samap/species.cds.fa genes/genes.CDS_only.gtf

## Map transcript to gene and clean headers
awk '{
    attr="";
    for(i=9;i<=NF;i++) attr=attr $i " ";
    match(attr,/gene_id "([^"]+)"/,g); gid=g[1];
    match(attr,/transcript_id "([^"]+)"/,t); tid=t[1];
    if(gid!="" && tid!="") print tid"\t"gid
}' genes/genes.CDS_only.gtf > samap/tx2gene.tsv

awk 'NR==FNR{m[$1]=$2; next}
/^>/{sub(/^>/,"",$0); tid=$1; gid=m[tid]; print ">"gid"; next}
{print}' samap/tx2gene.tsv samap/species.cds.fa > samap/species.cds.gene.tx.fa

## Keep only the longest transcript per gene
awk '/^>/{if(seq){print length(seq)"\t"hdr"\n"seq; seq=""} hdr=$0; next}{seq=seq$0}END{print length(seq)"\t"hdr"\n"seq}' \
  samap/species.cds.gene.tx.fa \
| awk 'BEGIN{RS=">"; ORS=""} NR>1{
      split($0,h,"\n"); hdr=h[1]; seq="";
      for(i=2;i<=NF;i++) seq=seq h[i];
      len=length(seq);
      split(hdr,ids,"|"); gid=ids[1];
      if(!(gid in keep) || len > maxlen[gid]){
          keep[gid]=">"hdr"\n"seq"\n"; maxlen[gid]=len
      }
} END{for(k in keep) printf "%s", keep[k]}' > samap/species.cds.longest.fa

## Optional
sed '/^>/! s/[^ACGTNacgtn]//g' samap/species.cds.longest.fa \
  | awk '/^>/ {print; next} {print toupper($0)}' \
  | fold -w 80 > samap/species.cds.cleaned.fa