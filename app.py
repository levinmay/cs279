import numpy as np
import os
import csv
from flask import Flask, request, render_template


app = Flask(__name__, static_url_path='',
            static_folder='static_folder',
            template_folder='template_folder')


# Set the data name to the name of the folder in the data folder that has your PDB files
dataName = "Hemo"
image_folder = os.path.join(app.static_folder, f'data/{dataName}/images')
df_folder = os.path.join(app.static_folder, f'data/{dataName}/dfs')


app.config["DEBUG"] = True

@app.route('/')
def homepage():
    """
    Renders the homepage of the Flask application.

    Returns:
        str: Rendered HTML template.
    """
    # Get a list of image files in the image folder
    image_files = [f for f in os.listdir(image_folder) if os.path.isfile(os.path.join(image_folder, f)) and f.lower().endswith('.png')]

    # Generate RMSD data for bar plot
    data_rmsd, labels_rmsd = generateBarPlotData("rmsd")
    # Generate SASA data for bar plot
    data_sasa, labels_sasa = generateBarPlotData("sasa")
    # Generate TM data for bar plot
    data_TM, labels_TM = generateBarPlotData("TM")
    
    return render_template(
        'index.html',
        dataName = dataName,
        data_rmsd = data_rmsd,
        labels_rmsd = labels_rmsd,
        data_sasa = data_sasa,
        labels_sasa = labels_sasa,
        data_TM = data_TM,
        images=image_files
    )

def generateBarPlotData(typeName):
    """
    Generates bar plot data from a CSV file.

    Args:
        typeName (str): Type of data to generate (e.g., "rmsd", "sasa", "TM").

    Returns:
        tuple: A tuple containing two lists - data and labels.
    """
    data = []
    labels =  []
    with open(f'{df_folder}/{typeName}_data.csv') as csvfile:
        reader = csv.reader(csvfile)
        next(reader) 
        for row in reader:
            label, score = row
            labels.append(label) 
            data.append(float(score)) 
    return data, labels



if __name__ == '__main__':
    app.run(debug=True)