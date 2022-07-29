# Allele frequency analysis

```
# Extract and format CODIS/disease freqs
./extract.sh codis
./format.py known_afreqs.txt $(ls freqs/codis*.tab) > called_freqs_codis.tab

./extract.sh disease
./format.py disease_afreqs.txt $(ls freqs/disease*.tab) > called_freqs_disease.tab
```

Notes:
* /storage/hziaeija/ensemble/notebooks/codis.ipynb # get called_freqs.txt for the notebook
