import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_html_components as html
import plotly.graph_objects as go

import pandas as pd
import numpy as np

from models.expression import expression
from helper_functions import *
from tabs.default_data import *

from make_app import app

cluster_dict_levels = {
    '1': [
         {
            'x' : 0,
            'y' : 0,
            'z' : 0
         },
         {
            'x' : 1,
            'y' : 0,
            'z' : 0
         }       
    ],
    '2': [
        {
            'x' : 0,
            'y' : 0,
            'z' : 0
        },
        {
            'x' : 1,
            'y' : 0,
            'z' : 0
        },
        {
            'x' : 0,
            'y' : 1,
            'z' : 0
        },
        {
            'x' : 1,
            'y' : 1,
            'z' : 0
        }
    ],
    '3': [
        {
            'x' : 0,
            'y' : 0,
            'z' : 0
        },
        {
            'x' : 1,
            'y' : 0,
            'z' : 0
        },
        {
            'x' : 0,
            'y' : 1,
            'z' : 0
        },
        {
            'x' : 1,
            'y' : 1,
            'z' : 0
        },
        {
            'x' : 0,
            'y' : 0,
            'z' : 1
        },
        {
            'x' : 1,
            'y' : 0,
            'z' : 1
        },
        {
            'x' : 0,
            'y' : 1,
            'z' : 1
        },
        {
            'x' : 1,
            'y' : 1,
            'z' : 1
        }
    ]
}
cluster_to_count_dict = {
    1:2,
    2:4,
    3:8
}

placeholder_dd = default_exp.gene_names_dd
   
cluster_label = html.Label(
    children = 'Select Gene',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
cluster_dd = dcc.Dropdown(
    id='cluster_selection_dd',
    options=default_exp.gene_names_dd,
    multi=False
)

cluster_center_label = html.Label(
    children = 'Select Group as Center (for coloring)',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
cluster_center_dd = dcc.Dropdown(
    id='cluster_center_selection_dd',
    options=default_exp.group_columns_dd,
    multi=False
)

update_button = html.Button(
    'Update', 
    id='update_button_cluster',
)

cluster_layout = go.Layout(
	title='',
    height=800
)
cluster_fig = go.Figure(
	data = None,
    layout = cluster_layout
)
cluster_graph = dcc.Graph(
    id='cluster',
    figure=cluster_fig,
)

input_combo_label  = html.Label(
    children = 'Input number of combinations in dataset',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
input_combo =  dcc.Input(
    id='combo_number',
    type='number',
    value=3,
    min=1,
    max=3
)

combo_labels_and_dds = []
combo_name_list = [
 	'Control Group', 
	'Condition 1',
	'Condition 2', 
    'Condition 1 and Condition 2', 
	'Condition 3', 
	'Condition 1 and Condition 3', 
	'Condition 2 and Condition 3', 
	'Condition 1 Condition 2, and Condition 3', 
]

for i, combo in enumerate(combo_name_list):
	combo_labels_and_dds.append(
		html.Label(
        	children = combo,
        	style={
            	'textAlign': 'left',
           		'color': 'black',
        	}
    	)
	)
	combo_labels_and_dds.append(
		dcc.Dropdown(
        	id=f'input_combo_{i+1}',
        	options=default_exp.group_columns_dd,
        	multi=False,
    	)
	)

body_2 = html.Div(
	[ 
		dbc.Container(
        	[   
                dbc.Row(
                    [
                    	dbc.Col(
	                        [
	                            cluster_label,
	                            cluster_dd,	      
                                cluster_center_label,
                                cluster_center_dd                
	                        ]
                    	),
                	]
            	),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                cluster_graph
                            ]
                    	)
                	]
            	),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                input_combo_label,
                                html.Br(),
                                input_combo,
                                html.Br(),
                            ] + 
                            combo_labels_and_dds +
                            [
                            	update_button,
                            ]
                        )
                    ]
                )
            ]
        )
    ]
)

tab_2 = dcc.Tab(
	label='Hypercube', 
	children=[
        body_2
    ]
)

@app.callback(
    [
        Output('cluster', 'figure')
    ],
    [
        Input('cluster_selection_dd', 'value'),
        Input('cluster_center_selection_dd', 'value'),
    ],
    [
        State('cluster_data', 'data'),

    ]
)
def update_cluster_graph(cluster_dd, center_group, cluster_data):
    if(cluster_data is None or cluster_dd is None):
        raise PreventUpdate
    print(cluster_data['meta'])
    print(cluster_data['data'].keys())
    mean_vals = []

    figure_params = {
        'x' : [],
        'y' : [],
        'z' : [],
        'text' : [],
        'hovertext' : [],
        'hoverinfo' : 'text',
        'mode' : 'markers+text',
        'textfont_size' : 15,
        'marker' : {
            'size' : 15,
            'line' : {
                'width' : 1
            },
        },
    }
    size_multiplier = 2

    for cluster in cluster_data['data'].keys():
        coord = cluster_data['meta']['coord'][cluster]
        figure_params['x'].append(coord['x'])
        figure_params['y'].append(coord['y'])
        figure_params['z'].append(coord['z'])
        figure_params['text'].append(cluster)

        value_dict = cluster_data['data'][cluster][cluster_dd]
        mean_vals.append(value_dict['mean'])
        hovertext = f"Mean Value: {value_dict['mean']:.2f}<br>Sample Count: {value_dict['count']} \
            <br>Min: {value_dict['min']:.2f}<br>Max: {value_dict['max']:.2f}<br>Range: {value_dict['range']:.2f}"
        figure_params['hovertext'].append(hovertext)
    figure_params['marker']['color'] = mean_vals
    if(center_group == 'Auto'):
        figure_params['marker']['colorscale'] = 'RdBu_r'
    else:
        center_val = cluster_data['data'][center_group][cluster_dd]['mean']
        adj_center_val = abs((center_val - min(mean_vals))/(max(mean_vals) - min(mean_vals)))
        figure_params['marker']['colorscale'] = [
            [0, 'rgb(0, 0, 255)'],   
            [adj_center_val, 'rgb(255, 255, 255)'],  
            [1, 'rgb(255, 0, 0)']]


    figure_params['marker']['colorbar'] = {
        "thickness" : 15
    }

    cluster_layout = go.Layout(
        height=800,
        title={
            'text': cluster_dd,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        }
    )
    cluster_fig = go.Figure(
        data = go.Scatter3d(figure_params), 
        layout = cluster_layout
    )
    cluster_fig.update_layout(
        plot_bgcolor='rgb(230, 235, 235, 0.5)',
        scene=dict(
            xaxis=dict(showticklabels=False, showgrid=False, title=cluster_data['meta']["conditions"][1], showline= False,showspikes = False, showbackground=False),
            yaxis=dict(showticklabels=False, showgrid=False, title=cluster_data['meta']["conditions"][2], showline= False,showspikes = False, showbackground=False),
            zaxis=dict(showticklabels=False, showgrid=False, title=cluster_data['meta']["conditions"][4], showline= False,showspikes = False, showbackground=False),
        ),
        font=dict(
            size=18
        )
    )

    return [cluster_fig]


@app.callback(
    [
        Output('cluster_data', 'data'),
        Output('cluster_center_selection_dd', 'options'),
        Output('cluster_center_selection_dd', 'value'),
    ],
    [
        Input('update_button_cluster', 'n_clicks'),
        Input('current_data', 'data'),
    ],
    [
        State('combo_number', 'value'),
        State('input_combo_1', 'value'),
        State('input_combo_2', 'value'),
        State('input_combo_3', 'value'),
        State('input_combo_4', 'value'),
        State('input_combo_5', 'value'),
        State('input_combo_6', 'value'),
        State('input_combo_7', 'value'),
        State('input_combo_8', 'value'),
    ]
)
def update_cluster_data(n_clicks, current_data, combo_number, ic_1, ic_2, ic_3, ic_4, ic_5, ic_6, ic_7, ic_8):
    if(current_data is None):
        raise PreventUpdate
    if(current_data == 'default'):
        exp = default_exp
    else:
        exp = expression(pd.DataFrame(None))
        exp.organize_from_store_upload(current_data)

    ic_clusters = [ic_1, ic_2, ic_3, ic_4, ic_5, ic_6, ic_7, ic_8]

    if(not any(ic_clusters)):
        ic_clusters = [exp.group_dict[x] for x in exp.group_dict.keys()]

    ic_clusters = ic_clusters[:cluster_to_count_dict[combo_number]]
    ic_cluster_names = [find_group_names(x, exp.group_dict) for x in ic_clusters]
    ic_clusters = [clean_column_names(x) for x in ic_clusters]

    cluster_data = {}
    cluster_meta = {
        'combo_number' : combo_number,
        'coord' : {},
        'conditions': []
    }

    count = 0
    for cluster, name in zip(ic_clusters, ic_cluster_names):
        r_df = exp.df[cluster]
        mean = r_df.mean(axis=1).values
        max_vals = r_df.max(axis=1).values
        min_vals = r_df.min(axis=1).values
        range_vals = max_vals - min_vals
        temp_df = pd.DataFrame(
            {
                'mean' : mean,
                'max' : max_vals,
                'min' : min_vals,
                'range' : range_vals,
                'count' : [len(cluster)] * len(max_vals)
            }
        )
        cluster_data[name] = dict(zip(list(exp.gene_names), temp_df.to_dict(orient = 'records')))
        cluster_meta['coord'][name] = cluster_dict_levels[str(combo_number)][count]
        cluster_meta['conditions'].append(name)
        count += 1

    while(len(cluster_meta['conditions']) < 8):
        cluster_meta['conditions'].append("")
    cluster_mem = {
        'meta':cluster_meta,
        'data':cluster_data
    }
    columns = []

    center_options = [{'label': i, 'value': i} for i in ['Auto'] + cluster_meta['conditions']]
    return [cluster_mem, center_options, 'Auto']

