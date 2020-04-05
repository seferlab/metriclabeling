# Data

`mlabel download` fills `data/raw/<source>/`; `mlabel build` caches processed
datasets in `data/processed/<source>/`.  Both folders are git-ignored.

* `original` (default): the 2011 snapshot published with the paper at
  <https://github.com/emresefer/metriclabeling/tree/master/data>
  (BioGRID 3.0.65, GO annotation files of January 2011, GO OBO of May 2010,
  Gasch et al. 2000 expression compendium).
* `current`: latest BioGRID / GO releases (URLs in `configs/species.yaml`).

To use the exact GO term lists of the paper instead of the automatic
selection, put one GO id per line in a file and pass `terms_file: path` in the
dataset entry of an experiment config.
