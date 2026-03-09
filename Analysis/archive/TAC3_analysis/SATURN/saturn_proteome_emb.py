## --------------------------------------------
## Human proteome embeddings
## --------------------------------------------
CUDA_VISIBLE_DEVICES=0
python esm/scripts/extract_single_gpu.py \
    esm1b_t33_650M_UR50S \
    Homo_sapiens.GRCh38.pep.all_clean.fa \
    ../proteome/embeddings/Homo_sapiens.GRCh38.pep.all_clean.fa_esm1b \
    --include mean \
    --gpu_id 0 \
    --truncate

python ../SATURN/protein_embeddings/map_gene_symbol_to_protein_ids.py \
    --fasta_path Homo_sapiens.GRCh38.pep.all.fa \
    --save_path Homo_sapiens.GRCh38.gene_symbol_to_protein_ID.json

python ../SATURN/protein_embeddings/convert_protein_embeddings_to_gene_embeddings.py \
    --embedding_dir embeddings/Homo_sapiens.GRCh38.pep.all_clean.fa_esm1b \
    --gene_symbol_to_protein_ids_path Homo_sapiens.GRCh38.gene_symbol_to_protein_ID.json \
    --embedding_model ESM1b \
    --save_path embeddings/Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM1b.pt

## --------------------------------------------
## Macaque proteome embeddings
## --------------------------------------------
python esm/scripts/extract.py \
    esm1b_t33_650M_UR50S \
    Macaca_mulatta.Mmul_10.pep.all_clean.fa \
    ../proteome/embeddings/Macaca_mulatta.Mmul_10.pep.all_clean.fa_esm1b \
    --include mean \
    --truncate \
    --nogpu

python ../SATURN/protein_embeddings/map_gene_symbol_to_protein_ids.py \
    --fasta_path Macaca_mulatta.Mmul_10.pep.all.fa \
    --save_path Macaca_mulatta.Mmul_10.gene_symbol_to_protein_ID.json

python ../SATURN/protein_embeddings/convert_protein_embeddings_to_gene_embeddings.py \
    --embedding_dir embeddings/Macaca_mulatta.Mmul_10.pep.all_clean.fa_esm1b \
    --gene_symbol_to_protein_ids_path Macaca_mulatta.Mmul_10.gene_symbol_to_protein_ID.json \
    --embedding_model ESM1b \
    --save_path embeddings/Macaca_mulatta.Mmul_10.gene_symbol_to_embedding_ESM1b.pt

## --------------------------------------------
## Marmoset proteome embeddings
## --------------------------------------------
CUDA_VISIBLE_DEVICES=1
python esm/scripts/extract_single_gpu.py \
    esm1b_t33_650M_UR50S \
    Callithrix_jacchus.mCalJac1.pat.X.pep.all_clean.fa \
    ../proteome/embeddings/Callithrix_jacchus.mCalJac1.pat.X.pep.all_clean.fa_esm1b \
    --include mean \
    --truncate \
    --gpu_id 1

python ../SATURN/protein_embeddings/map_gene_symbol_to_protein_ids.py \
    --fasta_path Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa \
    --save_path Callithrix_jacchus.mCalJac1.pat.X.gene_symbol_to_protein_ID.json

python ../SATURN/protein_embeddings/convert_protein_embeddings_to_gene_embeddings.py \
    --embedding_dir embeddings/Callithrix_jacchus.mCalJac1.pat.X.pep.all_clean.fa_esm1b \
    --gene_symbol_to_protein_ids_path Callithrix_jacchus.mCalJac1.pat.X.gene_symbol_to_protein_ID.json \
    --embedding_model ESM1b \
    --save_path embeddings/Callithrix_jacchus.mCalJac1.pat.X.gene_symbol_to_embedding_ESM1b.pt

## --------------------------------------------
## Mouse proteome embeddings
## --------------------------------------------
python esm/scripts/extract_single_gpu.py \
    esm1b_t33_650M_UR50S \
    Mus_musculus.GRCm39.pep.all_clean.fa \
    ../proteome/embeddings/Mus_musculus.GRCm39.pep.all_clean.fa_esm1b \
    --include mean \
    --truncate \
    --gpu_id 0

python ../SATURN/protein_embeddings/map_gene_symbol_to_protein_ids.py \
    --fasta_path Mus_musculus.GRCm39.pep.all.fa \
    --save_path Mus_musculus.GRCm39.gene_symbol_to_protein_ID.json

python ../SATURN/protein_embeddings/convert_protein_embeddings_to_gene_embeddings.py \
    --embedding_dir embeddings/Mus_musculus.GRCm39.pep.all_clean.fa_esm1b \
    --gene_symbol_to_protein_ids_path Mus_musculus.GRCm39.gene_symbol_to_protein_ID.json \
    --embedding_model ESM1b \
    --save_path embeddings/Mus_musculus.GRCm39.gene_symbol_to_embedding_ESM1b.pt