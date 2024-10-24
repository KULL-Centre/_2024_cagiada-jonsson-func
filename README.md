# _2024_cagiada-jonsson-func

# Decoding molecular mechanisms for loss of function variants in the human proteome

This repository contains Python code, [Jupyter](http://jupyter.org) Notebooks and data to generate new predictions and reproduce the results of the scientific paper `Decoding molecular mechanisms for loss of function variants in the human proteome` by M. Cagiada, N. Jonsson and K. Lindorff-Larsen, 2024, biorxiv ([link](https://doi.org/10.1101/2024.05.21.595203)).

## Layout
- `FunC-ESMs.ipynb` Jupyter Notebook to generate new predictions with FunC-ESMs, it can be run using the Google Colaboratory system.
- `Download_predictions.ipynb` Jupyter Notebook to download FunC-ESMs predictions of the human proteome, it can be run using the Google Colaboratory system.

- `script_and_analyses` Folder with all the data, notebooks and figures used to generate the current version of the manuscript. The data are collected in different folders:
  - `data` folder with all the data to reproduce the results of the analyses;
  - `notebooks` folder with all the Jupyter notebooks used to process the FunC-ESMs prediction and for the different analyses.
  - `scripts` folder with all the Python scripts used to process the FunC-ESMs predictions.

**N.B.: due to the large size, the FunC-ESM predictions (please ADD how to download the full df here), the thermodynamic stability prediction from [Tsuboyama,2023](https://www.nature.com/articles/s41586-023-06328-6) and MAVE score assays from [Proteingym](https://proteingym.org/) are not included in the Github repository. If you wish to fully reproduce the analyses, please download the data from the original source (KU ERDA or zenodo) or manuscript and move the files into the correct folders.

## FunC-ESMS human proteome predictions

### Unified Dataframe
The input data for the model and the predictions used in this study are available at the Electronic Research Data Archive at University of Copenhagen (KU/UCPH) (ERDA) and available via this [link](https://sid.erda.dk/cgi-sid/ls.py?share_id=DUWFpyjZp0) in a single H5PY file. 

To process the H5PY unified proteome prediction dataframe, due to the large size of the file, we suggest using the Python library [Vaex](https://vaex.io/), which works in a similar way to the Pandas library.

Basic commands for read and manipulate the H5PY dataframe with Vaex:
- to open your file you can use the command: 
```
df_unified=vaex.open(your_h5py_df_location)
```
- to extract column names use the command:
```
df_unified.get_column_names()
```
- to extract certain columns or to filter in a similar way to the Pandas, using the commands

```
column_data=df_unified['column_name']
filter_df=df_unified[df_unified['column_name']>=filter_value]
```
For all other advance command combinations and how to install Vaex, please refer to the software documentation.

### Single protein prediction analysis and download

Single prediction for each protein in the human proteome, alongside with a quick overview of the results can also be accessed via using this [colaboratory notebook](https://colab.research.google.com/github/KULL-Centre/_2024_cagiada-jonsson-func/blob/main/Download_predictions.ipynb)

    
## FunC-ESMs prediction notebook
### Usage
To run the prediction Notebook, we suggest to open it in [Google Colaboratory](https://colab.research.google.com/), it can be open clicking [here](https://colab.research.google.com/github/KULL-Centre/_2024_cagiada-jonsson-func/blob/main/FunC_ESMs.ipynb).
To run the Notebook is necessary to follow different steps:
1 Run the cells with the PRELIMINARY OPERATIONS prefix, this will install all the dependencies and load the functions required to run the predictor.
2 Upload all the structure and additional information in the `DATA UPLOADING` cell. Specifically you are required to:
  - set up a job name;
  - enter the query sequence for the the target protein;
  - download or upload an AlphaFold2 structure (via download from AF2 DB or upload from local storage);
  - specify the chain to be used as the main target in the predidction from the uploaded PDB (chain sequence and input sequence have to match);
  - if you wish to use a complex, tick the is_complex box and specify which chains should be used to use in the prediction (optional).
3 Run the `MODEL RUN` cell to first run the ESM models (ESM-1b and ESM-IF) and then to use FunC-ESMs and generate the variant and residue labels.
NOTE: for a new prediction, the `PRELIMINARY OPERATIONS` does not need to be run again and the new prediction can be run by loading the new target protein information (including a different job name) into the `DATA UPLOADING` cell and running the  `MODEL RUN` cell again.
4 Run the `VISUALISE RESULTS` cells  to display the different predictions from the ESMs models and FunC-ESMs as heatmaps (both at variant and residue level);
5 Run the `DOWNLOAD ALL RESULTS` cell to compress and download all the files generated during the current session.

### Input requirements:

- The input structure must be a predicted structure from AlphaFold2 if possible. This is necessary to ensure maximum prediction quality as ESM-IF has been trained using the backbone from AlphaFold2.

- The maximum size of proteins that can be used in the ESM-1B ESM-IF Colab implementation is 1023 residues, larger proteins cannot currently be handled by the Colab implementation.

### Output:

When a prediction is complete, the files associated with the run are stored in the `/content/FunC-ESMS_outputs` folder. When the `DOWNLOAD RESULTS` cell is executed, all files are downloaded at once.

The output files saved for each run are
- the AF2 pdb used as input
- the plots generated for each protein
- the prediction output files for each protein (ESM-1b,ESM-IF,FunC-ESMs)

The output files for each run will be label using the chosen `jobname`.

### Prediction file Format
The output can be generated as a csv or xlsx file and consists of several columns:
- For variant level prediction, the first column is the target 'MUTATION' which reports: wt amino acid, target substitution and position. The remaining columns describe the different scores produced (the column headers give the description of the score).
- For residue level predictions, the first column is the target 'RESIDUE' which reports: wt amino acid and position. The remaining columns describe the different scores produced (the column headers give the description of the score).


### Extra
#### License:

The source code and parameters of model are licensed under the permissive MIT licence.

#### Bugs:

For any bugs please report the issue on the project [Github](https://github.com/KULL-Centre/_2024_cagiada-jonsson-func) or contact one of the listed authors in the connected [manuscript](https://doi.org/).

#### Citing this work:

If you use our model please cite:

Cagiada, M., Jonsson, N., & Lindorff-Larsen, K. (2024). Decoding molecular mechanisms for loss of function variants in the human proteome. In bioRxiv (p. 2024.05.21.595203). https://doi.org/10.1101/2024.05.21.595203

```
@ARTICLE{Cagiada2024-wp,
  title    = "Decoding molecular mechanisms for loss of function variants in
              the human proteome",
  author   = "Cagiada, Matteo and Jonsson, Nicolas and Lindorff-Larsen, Kresten",
  journal  = "bioRxiv",
  pages    = "2024.05.21.595203",
  month    =  may,
  year     =  2024,
  language = "en"
```
