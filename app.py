import numpy as np
import os
import csv
from flask import Flask, request, render_template
from flask import render_template, redirect, url_for



app = Flask(__name__, static_url_path='',
            static_folder='static_folder',
            template_folder='template_folder')


dataName = "TTR"
image_folder = os.path.join(app.static_folder, f'data/{dataName}/images')
df_folder = os.path.join(app.static_folder, f'data/{dataName}/dfs')


app.config["DEBUG"] = True

@app.route('/')
def homepage():
    

    labels2 = ['April', 'May', 'June']
    data2 = [15, 25, 35]

    

    image_files = [f for f in os.listdir(image_folder) if os.path.isfile(os.path.join(image_folder, f))]

    #generate RMSD data
    data_rmsd, labels_rmsd = generateBarPlotData("rmsd")
    data_sasa, labels_sasa = generateBarPlotData("sasa")
    
    return render_template(
        'index.html',
        dataName = dataName,
        data_rmsd = data_rmsd,
        labels_rmsd = labels_rmsd,
        data_sasa = data_sasa,
        labels_sasa = labels_sasa,
        images=image_files
    )

def generateBarPlotData(typeName):
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