# CellTypist models

CellTypist annotation is **optional** and disabled by default.

## Listing available models

```python
import celltypist
celltypist.models.models_description()
```

## Downloading a model

```python
import celltypist
celltypist.models.download_models(model="Immune_All_Low.pkl")
```

Models default to `~/.celltypist/data/models/`.

## Running the pipeline with annotation

```bash
nextflow run main.nf \
    --samplesheet test/samplesheet.csv \
    --prebuilt_index ref/index.idx --t2g ref/t2g.txt \
    --run_celltypist \
    --celltypist_model_path ~/.celltypist/data/models/Immune_All_Low.pkl
```

The model file must already be downloaded; the pipeline does not fetch models
to keep runs fully offline-reproducible.
