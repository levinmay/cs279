
# Project README

## Overview

This project utilizes `pymol_code` to generate a script that is executed through PyMOL. The script populates data in the `static_folder`, which is then used by `app.py` to create a Flask-based webpage. The webpage is currently hosted [here](http://mayl.pythonanywhere.com). To apply this to your own protein data, you need to populate the `data` folder with your PDB files. Running the `pymol_code` from the PyMOL desktop application triggers the update process, and the Flask webpage will automatically reflect these updates on your hosting platform.

## Steps to Use

1. **Prepare Your Protein Data:**
   - Place your PDB files in the `data` folder.

2. **Running the PyMOL Script:**
   - Open the PyMOL application.
   - Execute the `pymol_code` script.
   - This script processes your PDB files and fills the `static_folder` with the necessary data.

3. **Generating the Flask Webpage:**
   - The `app.py` file uses the data from `static_folder` to create a Flask-oriented webpage.

## Authors

- May Levin
- Emma Williamson
- Joseph Allen

## References

- Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
  [TM-align](http://zhanglab.ccmb.med.umich.edu/TM-align/)
