import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import Omics.OmicsData as OD
import os

import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, FeatureAgglomeration
from sklearn.metrics import silhouette_score, calinski_harabaz_score

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.linear_model import SGDClassifier, LogisticRegression, RidgeClassifier, RandomizedLasso
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier, IsolationForest
from sklearn.tree  import DecisionTreeClassifier
import xgboost as xgb
from sklearn.neighbors import KNeighborsClassifier
from sklearn import svm
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score

from sklearn.feature_selection import RFE, f_regression
from minepy import MINE

RS = 20170628

app = dash.Dash()

CancerDataStore = pd.HDFStore('data/CancerDataStore_psn.h5')
cancerTypes = [x.replace('/', '') for x in CancerDataStore.keys()]
CancerDataStore.close()

#available_indicators = df['Indicator Name'].unique()

app.layout = html.Div([
    html.Div([

        html.Div([
            dcc.Dropdown(
                id='cancer-type',
                options=[{'label': i, 'value': i} for i in cancerTypes],
                value='BLCA'
            ),
        ],
        style={'width': '48%', 'display': 'inline-block'}),
    ]),

    dcc.Graph(id='pca-plot'),
])

@app.callback(
    dash.dependencies.Output('pca-plot', 'figure'),
    [dash.dependencies.Input('cancer-type', 'value')]
)
def update_graph(cancer_type):
    CancerDataStore = pd.HDFStore('data/CancerDataStore_psn.h5')
    dfCancerType = CancerDataStore.get(cancer_type)
    CancerDataStore.close()

    ClassVar = 'CancerStatus'

    dfCancerType = OD.dropNaNs(dfCancerType, ClassVar)

    dfAnalysis = dfCancerType.copy()

    VarLevelsToKeep = ['Primary solid Tumor', 'Solid Tissue Normal']
    
    logTransOffset = 1

    # filter and transform the data as specified above
    dfAnalysis_fl = OD.FilterLevels(dfAnalysis, ClassVar, VarLevelsToKeep, printStats='no')
    dfAnalysis_fl = OD.prepareDF(dfAnalysis_fl, ClassVar)


    # no gene filtering
    dfAnalysis_fl_cd = dfAnalysis_fl

    # Perform label encoding for the ClassVar and scale data using log transform
    dfAnalysis_fl_cd, ClassVarEncOrder = OD.mapClassVar(dfAnalysis_fl_cd,ClassVar)
    X, y = OD.fitLogTransform(dfAnalysis_fl_cd,logTransOffset)


    pca = PCA(n_components=2)
    X_r = pca.fit(X).transform(X)

    l = len(VarLevelsToKeep)
    colors = ['darkcyan', 'indigo', 'darkgreen', 'darkgoldenrod', 'darkblue','darkorange', 'mediumvioletred', 'crimson', 'darksalmon', 'darkred', 'cyan', 'orange','green']   
    colors = colors[0:l]


    traces = []
    for i in list(range(0,l)):
        traces.append(go.Scatter(
            x=X_r[y == i, 0],
            y=X_r[y == i, 1],
            mode='markers',
            opacity=0.7,
            marker={'size': 15,
                    'line': {'width': 0.5, 'color': 'white'}
            },
            name=VarLevelsToKeep[i]
        ))

    return {
        'data': traces,
        'layout': go.Layout(
            xaxis={'title': 'PCA 1'},
            yaxis={'title': 'PCA 2'},
            margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest'
        )
    }


if __name__ == '__main__':
    app.run_server(debug=True)