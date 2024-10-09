# pFBA simulations of Sulfolobus acidocaldarius
Code and data associated with the publication:
Sedlmayr, V. L., Széliová, D., De Kock, V., Gansemans, Y. G., Van Nieuwerburgh, F., Peeters, E., ... & Spadiut, O. Impact of nutrient excess on physiology and metabolism of Sulfolobus acidocaldarius. Frontiers in Microbiology, 15, 1475385.
DOI: https://doi.org/10.3389/fmicb.2024.1475385

The data contains pFBA simulations with Sulfolobus solfataricus (Sso) GSMM constrained with data from Sulfolobus acidocaldarius (Saci) from three experimental conditions: Low Cell Density, High Cell Density and Overfeeding.


## Folder structure
* code
    * data_preprocessing.ipynb - merging of experimental data, unit conversion
    * model_preprocessing.ipynb - GSMM of Sso corrected, metabolite names updated for escher visualization
    * FBA_sampling.ipynb - pFBA simulations
    * draw_fluxes_escher.ipynb - visualisation of predictions from FBA_sampling.ipynb
    * model_class.py - definition of Model() class + other functions

* data
    * experimental - fermentation data from Saci fermentations
        * ferm_data.csv
        * sulfolobus_rates.csv
        * sulfolobus_stdev.csv
    * central_fluxes_pfba.csv - pFBA results
    * escher_map.json - template for Escher plots
    * molar_masses.csv - molar masses of metabolites for unit conversion
    * node_names_escher.csv - metabolite names to be used in Escher
    * processed_rates.csv, processed_rates.pkl - rates form "data/experimental" folder processed with "data_preprocessing.ipynb", units mmol/(gDW h)
    * Sulfolobus_central_metabolism.csv - reaction IDs that we want to plot
    * sulfolobus_metabolite_ids.pkl, sulfolobus_reaction_ids.pkl - metabolite and reaction IDs from the original model before reanming in "model_preprocessing.ipynb"
    * Sulfolobus-solfataricus_corrected.json, Sulfolobus-solfataricus_corrected.xml - corrected GSMMs, outpu of "model_preprocessing.ipynb"

* figures
    * outputs of FBA_sampling.ipynb and draw_fluxes_escher.ipynb

* environment.yml - used python environment
