
## Activate virtual env
source /Users/hru/Dropbox/lag3/applications/cellphone_db_venv/bin/activate

cd /Users/hru/Dropbox/MelanoMAP/results/scrnaseq/cellphonedb/

## Use cellphonedb

cellphonedb method statistical_analysis --iterations=100 --threads=4 \
  --counts-data hgnc_symbol \
  --project-name li \
  input_files/li_meta.txt input_files/li_counts.txt




############################ newer

## Activate virtual env
source /Users/janihuuh/Dropbox/AML_TIM3/applications/cpdb-venv/bin/activate




## Use cellphonedb on li et al (new)
cd /Users/janihuuh/Dropbox/MelanoMAP/results/scrnaseq/li/new/cellphonedb/

cellphonedb method statistical_analysis --iterations=1000 --threads=25 \
  --counts-data hgnc_symbol \
  --project-name li_seurat_new \
  input_files/li_seurat_new_meta.txt input_files/li_seurat_new_counts.txt

## Use cellphonedb on durante et al (new)
cd /Users/janihuuh/Dropbox/MelanoMAP/results/scrnaseq/durante/new/cellphonedb/

cellphonedb method statistical_analysis --iterations=1000 --threads=25 \
  --counts-data hgnc_symbol \
  --project-name uv_total_new \
  input_files/uv_total_new_meta.txt input_files/uv_total_new_counts.txt





  ## Use cellphonedb on durante et al (new clusters new)
  cd /Users/janihuuh/Dropbox/MelanoMAP/results/scrnaseq/durante/new/cellphonedb/

  cellphonedb method statistical_analysis --iterations=1000 --threads=25 \
    --counts-data hgnc_symbol \
    --project-name uv_total_new \
    input_files/uv_total_new_clusters_new_meta.txt input_files/uv_total_new_clusters_new_counts.txt
