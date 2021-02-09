# -*- coding: utf-8 -*-
import pandas as pd
import dash_bio as dashbio
import pathlib
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import numpy as np
from scipy import stats
from dash.dependencies import Input, Output, State
from dash_extensions import Download
import os
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash.exceptions import PreventUpdate
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import math
import io
import base64
import plotly.figure_factory as ff
import plotly.express as px
from itertools import groupby
import datetime


NUMBER_OF_GENES_IN_VIOLIN_PLOT = 5

df = pd.read_csv(str(pathlib.Path("Gene_Expression_Data/GSE134896_FSKrepExpTable.csv")))

HEAT_GENES = 50

cluster_dict_three = {
    '3':[
    [[19, 20, 21], [19, 20, 21], [9, 10, 11], [13, 14, 15], [29, 30, 31], [25, 26, 27], [19, 20, 21], [19, 20, 21]],
    [[13, 13, 13], [40, 40, 40], [20, 20, 20],[30, 30, 30] ,[20, 20, 20], [30, 30, 30], [20, 20, 20], [30, 30, 30]] ],
    '2':[
    [[19, 20, 21], [13,14,15], [25, 26, 27], [19, 20, 21]],
    [[20, 20, 20], [30,30,30], [30,30,30],   [30, 30, 30]]
    ],
    '1':[
    [[19, 20, 21], [19, 20, 21]],
    [[20, 20, 20], [30, 30, 30]]
    ],
}
def load_t_test_cache():
    file_dir = pathlib.Path('Preterm Birth Gene Data/t_test_cache/')
    all_files = os.listdir(file_dir)
    cache = {}
    for i in range(len(all_files)):
        hold = pd.read_csv(file_dir / all_files[i])
        ##print(all_files[i])
        cache[all_files[i][:-4]] = hold

    return cache

#cache = load_t_test_cache()

row_labels = list(df['Gene_Name'])
row_labels_arr = np.array(row_labels)
row_dict = { row_labels[i] : i for i in range(0, len(row_labels) ) }

df_with_names = df
df = df.drop('Gene_Name', 1)
df_row_means = df.mean(axis=1)
df_row_stds = df.std(axis=1)

default_groups = {
    'CONTROL':['CTRL_0', 'CTRL_1', 'CTRL_2'],
    'P4':['P4_0', 'P4_1', 'P4_2'],
    'IL1B':['IL1B_0', 'IL1B_1', 'IL1B_2'],
    'FSK':['FSK_0', 'FSK_1', 'FSK_2'],
    'FSK_P4':['FSK_P4_0', 'FSK_P4_1', 'FSK_P4_2'],
    'P4_IL1B':['P4_IL1B_0', 'P4_IL1B_1', 'P4_IL1B_2'],
    'FSK_IL1B': ['FSK_IL1B_0','FSK_IL1B_1', 'FSK_IL1B_2',],
    'FSK_P4_IL1B': ['FSK_P4_IL1B_0', 'FSK_P4_IL1B_1','FSK_P4_IL1B_2']
}

default_options = [{'label': i, 'value': ' '.join(default_groups[i])} for i in default_groups.keys()]
column_options = [{'label': i, 'value': i} for i in df.columns]
gene_drop_down_options = default_options + column_options

single_gene_options = [{'label': i, 'value': i} for i in row_labels]

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.GRID])
app.config['suppress_callback_exceptions'] = True
server = app.server

tran_df = df.T

heatmap_Data = go.Heatmap(
        z=df,
        x=df.columns,
        y=row_labels, 
        colorscale= [
            [0, 'rgb(250, 250, 250)'],        #0
            [1./10000, 'rgb(250, 200, 200)'], #10
            [1./1000, 'rgb(250, 150, 150)'],  #100
            [1./100, 'rgb(250, 100, 100)'],   #1000
            [1./10, 'rgb(250, 50, 50)'],       #10000
            [1., 'rgb(250, 0, 0)'],             #100000
        ],    # y-axis labels
    )
layout_heatmap = go.Layout(
        title='',
        autosize=True,
        height=500,
        xaxis=go.layout.XAxis(
            showgrid=True,
        ),
        yaxis=go.layout.YAxis(
            showgrid=True,
        )
    )
fig_heatmap = go.Figure(data=None, layout=layout_heatmap)

fig_scatter = go.Figure(
        data = None,
            layout =  go.Layout(
            title='Scatter Plot',
            autosize=True,
            height=500,
            xaxis=go.layout.XAxis(
                showgrid=True,
            ),
            yaxis=go.layout.YAxis(
                showgrid=True,       
            ),
        )
    )

colors = {
        'background': 'white',
        'text': 'black'
    }

top_title = html.H1(
        children='Multi Group Gene Explorer (MUGGLE)',
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    )

heatmap_graph = dcc.Graph(
        id='heatmap',
        figure=fig_heatmap,
        style={
            'width':'49%',
            'display': 'inline-block'
        }
    )
large_scatter_graph = dcc.Graph(
        id='scatter_plot',
        style={
            'width':'49%',
            'display': 'inline-block'
        }
    )
x_drop_label = html.Label(
        children = 'Select x-axis grouping',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
x_drop_graph = dcc.Dropdown(
        id='x_axis_selection_dd',
        options=gene_drop_down_options,
        multi=True,
    )
y_drop_label = html.Label(
        children = 'Select y-axis grouping',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
y_drop_graph = dcc.Dropdown(
        id='y_axis_selection_dd',
        options=gene_drop_down_options,
        multi=True,
    )
combo_type_label = html.Label(
        children = 'Select combination type',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
combo_type_graph = dcc.Dropdown(
        id='combination_selection_dd',
        options=[
            {'label': 'mean', 'value': 'mean'},
            {'label': 'median', 'value': 'median'}
        ],
        value= 'mean',
        multi=False,
    )
update_button = html.Button(
        'Update', 
        id='update_button',
    )
violin_graph1 = dcc.Graph(
        id='violin_plot1',
        figure=fig_scatter
    )
violin_graph2 = dcc.Graph(
        id='violin_plot2',
        figure=fig_scatter
    )

violin_graph3 = dcc.Graph(
        id='violin_plot3',
        figure=fig_scatter
    )

violin_graph4 = dcc.Graph(
        id='violin_plot4',
        figure=fig_scatter
    )

x_group_id = html.Label(
        children = [html.Br(), 'Current X group'],
        id = 'x_group_id',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
y_group_id = html.Label(
        children = [html.Br(), 'Current Y group'],
        id = 'y_group_id',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
single_gene_selection_label = html.Label(
        children = 'Gene Focused or Condition Focused Comparison',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
single_gene_drop_graph = dcc.Dropdown(
        id='gene_selection_dd',
        options=single_gene_options,
        value='RNU86',
        multi=True,
    )
indivdual_or_gene_drop = dcc.Dropdown(
        id='invidual_or_gene_dd',
        options=[{'label': i, 'value': i} for i in ['Gene Focused', 'Condition Focused']],
        value='Condition Focused',
        multi=False,
    )
number_of_genes_label = html.Label(
        children = 'Bottom Right Display as Scatter Plot or Histogram',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
number_of_violin_genes = dcc.Dropdown(
        id='violin_gene_number',
        options=[{'label': i, 'value': i} for i in ['Histogram', 'Scatter Plot']],
        value='Scatter Plot',
        multi=False,
    )
genes_look_label = html.Label(
        children = 'Gene Selection',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
genes_look_dd = dcc.Dropdown(
        id='gene_select_dd',
        options=[{'label': i, 'value': i} for i in row_labels],
        value=row_labels[0],
        multi=True,
    )
var_label = html.Label(
        children = 'Select Variance Cut Off',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
fold_label = html.Label(
        children = 'Select Fold Cut Off',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
var_filter = html.Div(
    [
        dcc.Input(
            id='var_filter',
            type='number',
            value=1,
            min=0,
            max=100,
        )
    ]
)
download_heat_button = html.Div([html.Button("Download csv", id="btn"), Download(id="download")])

fold_filter = html.Div(
    [
        dcc.Input(
            id='fold_filter',
            type='number',
            value=1,
            min=0,
            max=1000
        )
    ]
)

body_1 = dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dcc.Store(id='heatmap_memory'),
                            dcc.Store(id='scatter_plot_memory'),
                            dcc.Store(id='name_memory'),
                            dcc.Store(id='main_data'),
                            dcc.Store(id='temp_data'),
                            dcc.Store(id='input_data'),
                            dcc.Store(id='default_data'),

                            top_title
                        ]
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            heatmap_graph,
                            large_scatter_graph
                        ],
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            single_gene_selection_label,
                            indivdual_or_gene_drop,
                            x_drop_label,
                            x_drop_graph, 
                            y_drop_label,
                            y_drop_graph,
                            genes_look_label,
                            genes_look_dd,
                            combo_type_label,
                            combo_type_graph,
                            number_of_genes_label,
                            number_of_violin_genes,
                            fold_label,
                            fold_filter,
                            var_label,
                            var_filter,
                            update_button,
                            x_group_id,
                            y_group_id,
                            download_heat_button,
                        ],
                        width = 4
                    ),
                    dbc.Col(
                        [
                            violin_graph1,
                        ], 
                        width = 8
                    ),
                ]
            )
        ],
        fluid = True
    )

top_title_upload = html.H1(
        children='Upload Data as a CSV',
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    )

update_button2 = html.Button(
        'Update', 
        id='update_button2',
    )
upload = html.Div([
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select Files')
        ]),
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    ),
    html.Div(id='output-data-upload'),
])
number_of_groups = dcc.Dropdown(
        id='number_of_groups_id',
        options=[{'label': i, 'value': i} for i in range(1, 10)],
        value=5,
        multi=False,
    )
number_of_groups = dcc.Dropdown(
        id='number_of_groups_id',
        options=[{'label': i, 'value': i} for i in range(1, 10)],
        value=5,
        multi=False,
    )
number_of_groups_labels = html.Label(
        children = 'Number of Custom Groups',
        style={
            'textAlign': 'center',
            'color': colors['text'],
        }
    )
update_msg = html.Label(
        children = 'Press the update button once all your custom groups and gene name column has been selected',
        style={
            'textAlign': 'center',
            'color': colors['text'],
        }
    )
input_msg = html.Div(id='input_msg')
variable_dropdowns = html.Div(id='variable_dropdowns')
variable_inputs = html.Div(id='variable_inputs')
gene_column_select_label = html.Div(id='gene_column_select_label')
gene_column_select = html.Div(id='gene_column_select')

loop_dict = {
    '1':2,
    '2':4,
    '3':8
    }
cluster_drop_label = html.Label(
        children = 'Select Gene',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
cluster_drop_graph = dcc.Dropdown(
        id='cluster_selection_dd',
        options=[{'label': i, 'value': i} for i in row_labels],
        value=row_labels[0],
        multi=False
    )
input_combo_label  = html.Label(
        children = 'Input number of combinations in dataset',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
control_label  = html.Label(
        children = 'Control Group',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_1_label  = html.Label(
        children = 'Condition 1',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_2_label  = html.Label(
        children = 'Condition 2',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_3_label  = html.Label(
        children = 'Condition 3',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_1_2_label  = html.Label(
        children = 'Condition 1 and Condition 2',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_1_3_label  = html.Label(
        children = 'Condition 1 and Condition 3',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_2_3_label  = html.Label(
        children = 'Condition 2 and Condition 3',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
condition_1_2_3_label  = html.Label(
        children = 'Condition 1, Condition 2, and Condition 3',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
input_combo =  dcc.Input(
            id='combo_number',
            type='number',
            value=3,
            min=1,
            max=1000
        )

input_combo_1 = dcc.Dropdown(
        id='input_combo_1',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )

input_combo_2 = dcc.Dropdown(
        id='input_combo_2',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
input_combo_3 = dcc.Dropdown(
        id='input_combo_3',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
input_combo_4 = dcc.Dropdown(
        id='input_combo_4',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
input_combo_5 = dcc.Dropdown(
        id='input_combo_5',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
input_combo_6 = dcc.Dropdown(
        id='input_combo_6',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
input_combo_7 = dcc.Dropdown(
        id='input_combo_7',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
input_combo_8 = dcc.Dropdown(
        id='input_combo_8',
        options=gene_drop_down_options,
        value="CONTROL",
        multi=False,
        #disabled=True,
    )
update_button_cluster = html.Button(
        'Update', 
        id='update_button_cluster',
    )
clusters = dcc.Graph(
        id='cluster',
        figure=fig_scatter
    )
body_3 = dbc.Container(
        [   
                dbc.Row(
                    [
                    dbc.Col(
                        [
                            cluster_drop_label,
                            cluster_drop_graph,
                            dcc.Store(id='cluster_info'),
                            dcc.Store(id='temp_cluster_info')
                            
                        ]
                    ),
                ]
            ),
                    dbc.Row(
                        [
                            dbc.Col(
                            [
                                clusters
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
                                control_label,
                                input_combo_1,
                                condition_1_label,
                                input_combo_2,
                                condition_2_label,
                                input_combo_3,
                                condition_1_2_label,
                                input_combo_4,
                                condition_3_label,
                                input_combo_5,
                                condition_1_3_label,
                                input_combo_6,
                                condition_2_3_label,
                                input_combo_7,
                                condition_1_2_3_label,
                                input_combo_8,
                                update_button_cluster,
                            ]
                    )
                ]
            ),
                    #dbc.Row(
                    #    [
                    #        dbc.Col(
                    #        [
                    #            html.Div(html.Img(src=app.get_asset_url('study_img1.jpg')))
                    #        ]
                   # )
                #]
            #)
        ]
    )
body_2 = dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            top_title_upload,
                        ]
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            upload,

                        ]
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            number_of_groups_labels,
                            number_of_groups,
                            input_msg,
                            gene_column_select_label,
                            gene_column_select,
                        ]
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            variable_inputs,
                        ]
                    ),
                    dbc.Col(
                        [
                            variable_dropdowns,
                        ]
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            update_msg,
                        ]
                    ),
                ]
            ),  
            dbc.Row(
                [
                    dbc.Col(
                        [
                            update_button2,
                        ]
                    ),
                ]
            ),  

        ])    
first_tab = html.Div([body_1])
second_tab = html.Div([body_2])
third_tab = html.Div([body_3])
app.layout = html.Div([
    dcc.Tabs([
        dcc.Tab(label='Specific Visuals', children=[
            first_tab
        ]),
        dcc.Tab(label='Clusters', children=[
            third_tab
        ]),
        dcc.Tab(label='Custom Data', children=[
            second_tab
        ]),

    ])
])



#Callback functions

@app.callback(
        [
            Output('heatmap_memory', 'data'),
            Output('scatter_plot_memory', 'data'),
            Output(component_id='x_group_id', component_property='children'),
            Output(component_id='y_group_id', component_property='children'),
            Output('name_memory', 'data'),
            Output("scatter_plot", "selectedData"),
        ],
        [
            Input('update_button', 'n_clicks'),

        ],
        state=[
            State('main_data', 'data'),
            State('x_axis_selection_dd', 'value'),
            State('y_axis_selection_dd', 'value'),
            State('gene_select_dd', 'value'),
            State('invidual_or_gene_dd', 'value'),
            State('violin_gene_number', 'value'),
            State('combination_selection_dd', 'value'),
            State('name_memory', 'data'),
            State('scatter_plot_memory', 'data'),
            State('fold_filter', 'value'),
            State('var_filter', 'value'),
            State('default_data', 'data'),

        ]
)
def update_button_main(n_clicks, main_data, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, violin_gene_number, combination_dd_value, name_memory, scatter_mem, fold_filter, var_filter, default_data):
    '''
        This callback takes all the dropdown menus and will calculate the differentially expressed genes and update the memory objects
    '''
    if(1==1):
        if(x_dd_value is None or y_dd_value is None):
            raise PreventUpdate
        x_num = len(x_dd_value)
        y_num = len(y_dd_value)

        if(x_num == 0 or y_num == 0):
            raise PreventUpdate
        if(default_data is None):
            defualt_data = load_default_data(1,1)
        if(main_data is None):
            main_data = defualt_data
        x_dd_value = combine_columns(x_dd_value, main_data)
        y_dd_value = combine_columns(y_dd_value, main_data)
        x_num = len(x_dd_value)
        y_num = len(y_dd_value)



        

        orig_fold = main_data['fold_filter']
        orig_var = main_data['var_filter']

        #print("orig_fold", type(orig_fold), orig_fold)
        #print("orig_var", type(orig_var), orig_var)

        #print("fold_filter", type(fold_filter), fold_filter)
        #print("var_filter", type(var_filter), var_filter)

        main_df = pd.DataFrame(main_data['main_df'])

        if(len(main_df) == 0 or not (set(x_dd_value).issubset(set(main_df.columns))) or not (set(y_dd_value).issubset(set(main_df.columns)))):
            raise PreventUpdate

        main_data = filter_by_var_fold(main_data, x_dd_value, y_dd_value, float(fold_filter), float(var_filter))
        main_df = pd.DataFrame(main_data['main_df'])

        x_dff = main_df.loc[:, x_dd_value]
        y_dff = main_df.loc[:, y_dd_value]
        if(individ_or_gene_dd == 'Condition Focused'):

            if(name_memory is not None and x_dd_value == name_memory['x_dd'] and y_dd_value == name_memory['y_dd']
                and combination_dd_value == name_memory['combo'] and individ_or_gene_dd == name_memory['type']
                and violin_gene_number == name_memory['violin_gene_number'] and abs(orig_fold - fold_filter) < 0.001
                and abs(orig_var - var_filter) < 0.001):

                gene_diff_df = pd.DataFrame(scatter_mem['df1'])

            else:
                gene_diff_df = find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value, main_data, True)

        else:

            gene_diff_df = find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value, main_data)
        start1 = datetime.datetime.now()
        heat_dict = create_heatmap_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data)
        start2 = datetime.datetime.now()
        #print("scatter_memeory time ",str(start2 - start1))
        scatter_dict = create_scatter_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data)

        start3 = datetime.datetime.now()

        #print("scatter_memeory time ",str(start3 - start2))
        x_display, y_display = create_x_y_group_display_names(x_dd_value, y_dd_value)

        name_memory = {
            'x_dd' : x_dd_value,
            'y_dd' : y_dd_value,
            'type' : individ_or_gene_dd,
            'combo' : combination_dd_value,
            'violin_gene_number':violin_gene_number
            }
    else:
        raise PreventUpdate
    return heat_dict, scatter_dict, x_display, y_display, name_memory, None

@app.callback(
        [
            Output('main_data', 'data'),
            Output('update_button', 'n_clicks'),
            Output('x_axis_selection_dd', 'options'),
            Output('y_axis_selection_dd', 'options'),
            Output('gene_select_dd', 'options'),
            Output('input_combo_1', 'options'),
            Output('input_combo_2', 'options'),
            Output('input_combo_3', 'options'),
            Output('input_combo_4', 'options'),
            Output('input_combo_5', 'options'),
            Output('input_combo_6', 'options'),
            Output('input_combo_7', 'options'),
            Output('input_combo_8', 'options'),
            Output('cluster_selection_dd', 'options'),
        ],
        [
            Input('temp_data', 'data'),
        ],
        state=[
            State('update_button', 'n_clicks')
        ]
        
)
def update_main_memory(temp_data, n_clicks):
    #print("IS THIS HAPPENING??????????????????????????")
    if(n_clicks is None):
        n_clicks = 0
    if(temp_data is None):
        main_data = load_default_data(1,1)
    else:
        main_data = temp_data

    return main_data, n_clicks+1, main_data['gene_drop_down_options'], main_data['gene_drop_down_options'], main_data['single_gene_options'], main_data['gene_drop_down_options'],main_data['gene_drop_down_options'],main_data['gene_drop_down_options'],main_data['gene_drop_down_options'],main_data['gene_drop_down_options'],main_data['gene_drop_down_options'],main_data['gene_drop_down_options'],main_data['gene_drop_down_options'], main_data['single_gene_options']

@app.callback(
        
            Output(component_id='heatmap', component_property='figure'),
        [
            Input('heatmap_memory', 'data'),
        ],

)
def update_heatmap(memory):
    '''
        This callback uses the heat map memory store and any selected data in the scatter plot to update the heatmap
    '''
    #print("staritng heat call")
    if(memory is None):
        raise PreventUpdate

    start1 = datetime.datetime.now()

    name_adj_heat_df = pd.DataFrame(memory['original'])
    name_adj_normal_heat_df = pd.DataFrame(memory['normal'])

    current_row_labels = memory['current_row_labels']
    #print("staritng heat call1")
    try:    
        max_val = name_adj_normal_heat_df.values.max()
        min_val = name_adj_normal_heat_df.values.min()
    except:
        max_val = 2
        min_val = -2
    heatmap_Data = go.Heatmap(
            z=name_adj_normal_heat_df,
            x=name_adj_normal_heat_df.columns,
            y=current_row_labels, 
            hoverinfo='text',
            text=get_hover_text(name_adj_heat_df, current_row_labels),
            colorscale='RdBu_r',
            zmax=max(max_val + 0.2, 0.5),
            zmin= min(min_val - 0.2, -0.5),
        )
    #print("staritng heat call2")
    '''
    [
                [0.0, "rgb(49,54,149)"],
                [0.1111111111111111, "rgb(69,117,180)"],
                [0.2222222222222222, "rgb(116,173,209)"],
                [0.3333333333333333, "rgb(171,217,233)"],
                [0.4444444444444444, "rgb(224,243,248)"],
                [0.5555555555555556, "rgb(254,224,144)"],
                [0.6666666666666666, "rgb(253,174,97)"],
                [0.7777777777777778, "rgb(244,109,67)"],
                [0.8888888888888888, "rgb(215,48,39)"],
                [1.0, "rgb(165,0,38)"]
                ]
    '''
    layout_heatmap = go.Layout(
            title={
                'text': "Heatmap",
                'y':0.9,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            autosize=True,
            height=500,
            xaxis=go.layout.XAxis(
                showgrid=True,
            ),
            yaxis=go.layout.YAxis(
                showgrid=True,
            )
        )
    #print("staritng heat call3")
    fig_heatmap = go.Figure(data=[heatmap_Data], layout=layout_heatmap)
    start2 = datetime.datetime.now()
    #print("heatmap callback graph", start2 - start1)
    #print("staritng heat call4")
    return fig_heatmap

one_combo = {
    'x':[0],
    'y':[0],
    'show':0
}
two_combo = {
    'x':[-1, 0],
    'y':[0, 0],
    'show': 1
}
three_combo = {
    'x':[-1, 0, 1],
    'y':[0, 0, 0],
    'show': 1
}
four_combo = {
    'x':[-1, 0, 1, -1],
    'y':[0, 0, 0, -1],
    'show': 1
}
five_combo = {
    'x':[-1, 0, 1, -1, 0],
    'y':[0, 0, 0, -1, -1],
    'show': 1
}
six_combo = {
    'x':[-1, 0, 1, -1, 0, 1],
    'y':[0, 0, 0, -1, -1, -1],
    'show': 1
}
combo_locs = {
    "1":one_combo,
    "2":two_combo,
    "3":three_combo,
    "4":four_combo,
    "5":five_combo,
    "6":six_combo,
}
def arrange_cluster_points(x_location_center, y_location_center, combo_group, start_num):
    combo_length = len(combo_group)
    combo_dict = combo_locs[str(combo_length)]
    x = list(map(lambda x : x + x_location_center, combo_dict['x']))
    y = list(map(lambda x : x + y_location_center, combo_dict['y']))
    show_val = combo_dict['show'] + start_num
    return x, y, show_val
@app.callback(
            Output('cluster_info', 'data'),
        [
            Input('update_button_cluster', 'n_clicks')
        ],
        [
            State('main_data', 'data'),
            State('input_combo_1', 'value'),
            State('input_combo_2', 'value'),
            State('input_combo_3', 'value'),
            State('input_combo_4', 'value'),
            State('input_combo_5', 'value'),
            State('input_combo_6', 'value'),
            State('input_combo_7', 'value'),
            State('input_combo_8', 'value'),
            State('combo_number', 'value'),

        ]
)
def update_cluster_button(n_clicks, main_data, input_combo_1, input_combo_2,input_combo_3, input_combo_4, input_combo_5, input_combo_6, input_combo_7, input_combo_8, input_combo):

    #print("AHHHHHHHHHHHHH")

    if(n_clicks == 0 or main_data is None):
        raise PreventUpdate
    if(not ((input_combo == 1 and input_combo_1 is not None and input_combo_2 is not None) or
       (input_combo == 2 and input_combo_1 is not None and input_combo_2 is not None and input_combo_3 is not None and input_combo_4 is not None) or
       (input_combo == 3 and input_combo_1 is not None and input_combo_2 is not None and input_combo_3 is not None and input_combo_4 is not None and input_combo_5 is not None and input_combo_6 is not None and input_combo_7 is not None and input_combo_8 is not None))):
        raise PreventUpdate

    input_combo_1 = combine_columns(input_combo_1.split(), main_data)
    input_combo_2 = combine_columns(input_combo_2.split(), main_data)
    input_combo_3 = combine_columns(input_combo_3.split(), main_data)
    input_combo_4 = combine_columns(input_combo_4.split(), main_data)
    input_combo_5 = combine_columns(input_combo_5.split(), main_data)
    input_combo_6 = combine_columns(input_combo_6.split(), main_data)
    input_combo_7 = combine_columns(input_combo_7.split(), main_data)
    input_combo_8 = combine_columns(input_combo_8.split(), main_data)
    x_locs = []
    y_locs = []
    single_row_names = []
    x_location_template = cluster_dict_three[str(input_combo)][0]
    y_location_template = cluster_dict_three[str(input_combo)][1]

    loop_num = loop_dict[str(input_combo)]
    input_combos = [input_combo_1, input_combo_2, input_combo_3, input_combo_4, input_combo_5, input_combo_6, input_combo_7, input_combo_8]
    orientation = ['middle left', 'top center', 'middle right'] * (loop_num)
    final_pos = []
    show_text = []
    total_num = 0
    for i in range(loop_num):
        current_combo = input_combos[i]
        x, y, show_val = arrange_cluster_points(x_location_template[i][1], y_location_template[i][1], current_combo, total_num)
        x_locs.extend(x)
        y_locs.extend(y)
        single_row_names.extend(current_combo)
        final_pos.extend(orientation[:len(current_combo)])
        show_text.append(show_val)

        total_num = total_num + len(current_combo)

    #print("finished updating cluster data")
    #print(orientation)

    return {'single_row_names':single_row_names, 'x_locs':x_locs, 
            'y_locs':y_locs, 'orientation':orientation, 'show_text': show_text}




@app.callback(
        
            Output(component_id='cluster', component_property='figure'),

        [
            Input('cluster_selection_dd', 'value'),
            Input('cluster_info', 'data'),

        ]
)
def update_cluster(cluster_dd, cluster_data):
    '''
        This callback uses the scatter plot memory store and any selected data in the heatmap to update the scatter plot
    '''
    #print("we doing something")
    if(cluster_data is None):
        #print("right here?")
        single_row = df_with_names[df_with_names['Gene_Name'] == cluster_dd]
        #print(single_row)
        all_columns = df.columns
        layout = go.Layout(
        title = 'Clusters',
        xaxis = go.layout.XAxis(
            title = '',
             tickmode = 'array',
             tickvals = [''],
             ticktext =[''],
             tickfont = go.Font(
                 color = 'rgb(255, 255, 255)'
             )
        ),
        yaxis = go.layout.YAxis(
            title = '',
             tickmode = 'array',
             tickvals = [''],
             ticktext =[''],
             tickfont = go.Font(
                 color = 'rgb(255, 255, 255)'
             )
        ),
        height=700,
        xaxis_showgrid=False, yaxis_showgrid=False
        )
        x = [9, 10, 11, 9, 10, 11, -1, 0, 1, 19, 20, 21, 15, 16, 17, 3, 4, 5, 9, 10, 11, 9, 10, 11]
        y = [-7, -7, -7, 20, 20, 20, 0,0,0, 0, 0, 0, 10, 10, 10, 10, 10, 10, 0, 0, 0, 10, 10, 10]
        x = list(map(lambda x : x + 10, x))
        y = list(map(lambda x : x + 20, y))
        markers = []
        hover_text = []
        single_row_list = list(single_row)
        modes = []
        all_vals = []
        for column in single_row.columns:
            entries = list(single_row[column])
            if(column != 'Gene_Name' and len(entries) > 0):
                ##print(single_row[column])
                ##print(single_row[column][0])
                
                all_vals.append(entries[0])
        counter = 0
        for column in single_row.columns:
            entries = list(single_row[column])
            if(column != 'Gene_Name' and len(entries) > 0):
                #markers.append(f'rgba(60,179,113,{normalize_gene(entries[0] + 0.01)})')#, min_value = min(all_vals), max_value = max(all_vals), log_transform = False)})')
                markers.append(f'rgba(60,179,113,{max([0.01,min([1, normalize_gene(entries[0] + 0.01, min_value = min(all_vals), max_value = max(all_vals), log_transform = False)])])})')
                hover_text.append(f"{column}<br>{round(entries[0], 1)}")
                if(counter == 1):
                    modes.append(f"{column[:-2]}")
                else:
                    modes.append(f"")
                
                counter = counter + 1
                if(counter == 3):
                    counter = 0
        #print(markers)
        fig = go.Figure(go.Scattergl(x =x, y=y, mode='markers+text', 
                                   text=modes,
                                 marker_size=25,
                                textposition= ['middle left', 'top center', 'middle right', 
                                               'middle left', 'bottom center', 'middle right', 
                                               'middle left', 'top center', 'middle right', 
                                               'middle left', 'top center', 'middle right', 
                                               'middle left', 'top center', 'middle right', 
                                               'middle left', 'top center', 'middle right', 
                                               'middle left', 'top center', 'middle right', 
                                               'bottom left', 'top center', 'bottom right'],
                                hovertext=hover_text,
                                marker_color=markers), layout=layout)
        #print("FINISHED")
    else:
        #print("starting updating figure")
        #print(cluster_dd)
        single_row = df_with_names[df_with_names['Gene_Name'] == cluster_dd]
        #print(single_row)
        #print("bah")
        #print(cluster_data['single_row_names'])
        all_columns = df.columns
        layout = go.Layout(
        title = 'Clusters',
        xaxis = go.layout.XAxis(
            title = '',
             tickmode = 'array',
             tickvals = [''],
             ticktext =[''],
             tickfont = go.Font(
                 color = 'rgb(255, 255, 255)'
             )
        ),
        yaxis = go.layout.YAxis(
            title = '',
             tickmode = 'array',
             tickvals = [''],
             ticktext =[''],
             tickfont = go.Font(
                 color = 'rgb(255, 255, 255)'
             )
        ),
        height=700,
        xaxis_showgrid=False, yaxis_showgrid=False
        )
        markers = []
        hover_text = []
        single_row_list = list(single_row)
        modes = []
        all_vals = []
        #print("COLUMN PRINT")
        for column in cluster_data['single_row_names']:
            #print(column)
            entries = list(single_row[column])
            if(column != 'Gene_Name' and len(entries) > 0):
                ##print(single_row[column])
                ##print(single_row[column][0])
                
                all_vals.append(entries[0])
        counter = 0
        for column in cluster_data['single_row_names']:
            entries = list(single_row[column])
            if(column != 'Gene_Name' and len(entries) > 0):
                #print(counter)
                #print(counter in cluster_data['show_text'])
                #markers.append(f'rgba(60,179,113,{normalize_gene(entries[0] + 0.01)})')#, min_value = min(all_vals), max_value = max(all_vals), log_transform = False)})')
                markers.append(f'rgba(60,179,113,{max([0.01,min([1, normalize_gene(entries[0] + 0.01, min_value = min(all_vals), max_value = max(all_vals), log_transform = False)])])})')
                hover_text.append(f"{column}<br>{round(entries[0], 1)}")
                if(counter in cluster_data['show_text']):
                    modes.append(f"{column[:-2]}")
                else:
                    modes.append(f"")
                
                counter = counter + 1

        #print(markers)
        #print(cluster_data['show_text'])
        #print(cluster_data['x_locs'])
        #print(cluster_data['y_locs'])
        #print(modes)
        fig = go.Figure(go.Scattergl(x =cluster_data['x_locs'], y=cluster_data['y_locs'], mode='markers+text', 
                                   text=modes,
                                 marker_size=25,
                                textposition= cluster_data['orientation'],
                                hovertext=hover_text,
                                marker_color=markers), layout=layout)
        #print("FINISHED")

    return fig


@app.callback(

            Output(component_id='scatter_plot', component_property='figure'),

        [
            Input('scatter_plot_memory', 'data'),
        ],
        [
            State('main_data', 'data'),
        ]
)
def update_scatter(memory, main_data):
    '''
        This callback uses the scatter plot memory store and any selected data in the heatmap to update the scatter plot
    '''
    if(memory is None):
        raise PreventUpdate


    start1 = datetime.datetime.now()
    stored_mem1 = memory['df1']
    stored_mem2 = memory.get('df2')

    if(isinstance(stored_mem2, list)):
        gene_diff_df = pd.DataFrame(stored_mem1)
        gene_diff_df.sort_values(by = ["absolute difference measure"], inplace=True, ascending=False)
        scatter_figure = create_scatter_diff_genes(gene_diff_df, stored_mem2[0], stored_mem2[1])

    else:
        stored_mem3 = memory['names']
        scatter_figure = create_gene_based_scatter(pd.DataFrame(stored_mem1), memory['selected_indices'], stored_mem3[0], stored_mem3[1], main_data)

    start2 = datetime.datetime.now()
    #print("scatter callback update time", start2 - start1)

    return scatter_figure

@app.callback(

        Output(component_id='violin_plot1', component_property='figure'),


        [
            Input("scatter_plot", "selectedData"),
            Input('heatmap_memory', 'data'),
        ],
        state=[
            State('violin_gene_number', 'value'),
            State('name_memory', 'data'),
            State('main_data', 'data'),
        ]
)

def update_violin(scatter_selected_data, memory, violin_gene_number, name_memory, main_data):
    '''
    This callback uses the heat map memory store and any selected data in the scatter plot and heat map to update the violin plot
    '''
    if(memory is None or main_data is None or not set(memory['x_dd']).issubset(set(pd.DataFrame((main_data['main_df'])).columns)) or 
        not set(memory['y_dd']).issubset(set(pd.DataFrame((main_data['main_df'])).columns))):
        raise PreventUpdate

    start1 = datetime.datetime.now()
    main_df = pd.DataFrame((main_data['main_df']))

    current_row_labels = memory['current_row_labels']

    selected_scatter_genes = []
    if(scatter_selected_data is not None):
        for point in scatter_selected_data['points']:
            selected_scatter_genes.append(point['text'])


        selected_scatter_set = frozenset(selected_scatter_genes)
        all_genes = [x for x in current_row_labels if x in selected_scatter_set]
        all_genes_set = set(all_genes)
        non_heat_genes = [x for x in selected_scatter_genes if x not in all_genes_set]
        current_row_labels = all_genes + non_heat_genes


    current_row_labels = current_row_labels[:10]

    violin_figure = update_violin_plot(current_row_labels, name_memory['x_dd'], name_memory['y_dd'], violin_gene_number, main_df, main_data['row_dict'])

    start2 = datetime.datetime.now()
    #print("violin callback update time", start2 - start1)

    return violin_figure

@app.callback(  
            [
                Output('input_data', 'data'),
                Output('input_msg', 'children')
            ],
            [
                Input('upload-data', 'contents')
            ],
            [
                State('upload-data', 'filename'),
                State('upload-data', 'last_modified')
            ])
def create_upload_display(contents, filename, last_modified):
    if(contents is None):
        raise PreventUpdate

    if(1==1):
        content_type, content_string = contents[0].split(',')
        decoded = base64.b64decode(content_string)
        input_df = pd.read_csv(
            io.StringIO(decoded.decode('utf-8')))
        return_msg = html.Div([
            f'{filename[0]} successfully uploaded'
        ])
        return_mem = create_new_data(input_df)
    else:# Exception as e:
        return_msg = html.Div([
            'There was an error processing this file.'
        ])
        return_mem = None

    return return_mem, return_msg

'''
@app.callback(
            [
                Output('x_axis_selection_dd', 'options'),
                Output('y_axis_selection_dd', 'options'),
                Output('gene_select_dd', 'options'),
            ],
            [
                Input('input_data', 'data')
            ],
            [ 
                State('x_axis_selection_dd', 'options'),
                State('y_axis_selection_dd', 'options'),
            ]
        )
def change_drop_down_options(input_data, x_dd_value, y_dd_value):
    if(input_data is None):
        raise PreventUpdate

    #print(input_data['gene_drop_down_options'])
    return input_data['gene_drop_down_options'], input_data['gene_drop_down_options'], input_data['single_gene_options']
'''
@app.callback(
    [
        Output('variable_dropdowns', 'children'),
        Output('variable_inputs', 'children'),

       
    ],
    [
        Input('number_of_groups_id', 'value'),
        Input('input_data', 'data')
    ]
)
def update_div(number_groups, input_data):
    if(input_data is None):
        raise PreventUpdate
    dropdowns = [dcc.Dropdown(
        id=f'default_group #{i}',
        options=input_data['gene_drop_down_options'],
        value=5,
        multi=True,
        style={
            'height': 50,
            'width':'100%',
            'display': 'inline-block'
        }
    ) for i in range (number_groups)]
    labels = [dcc.Input(
            id=f"input_{i}",
            type='text',
            placeholder="Input Group Name",
        style={
            'height': 54,
            'width':'100%',
            'display': 'inline-block'
        }
        ) for i in range(number_groups)]


    return dropdowns, labels

@app.callback(
    [
        Output('temp_data', 'data'),
    ],
    [
        Input('update_button2', 'n_clicks'),
    ],
    [
        State('input_data', 'data'),
        State('variable_dropdowns', 'children'),
        State('variable_inputs', 'children'),
    ]
)
def update_input_memory(n_clicks, input_data, dropdowns, input_names):
    if(input_data is None):
        return [load_default_data(1,1)]
        raise PreventUpdate

    gene_drop_down_options_hold = []
    default_groups_hold = {}
    for name, dd in zip(input_names, dropdowns):

        if(name['props'].get('value') is not None and dd['props'].get('value') is not None):
            gene_drop_down_options_hold.append({'label': name['props']['value'], 'value' : ' '.join(list(dd['props']['value']))})
            default_groups_hold[name['props']['value']] = [dd['props']['value']]
    #gene_drop_down_options_hold = [{'label': name['value'], 'value' : ' '.join(list(dd['value']))} for name, dd in zip(input_names, dropdowns)]
    #default_groups_hold = {name['value'] : [dd['value']] for name, dd in zip(input_names, dropdowns)}
    input_data['gene_drop_down_options'] = gene_drop_down_options_hold + input_data['gene_drop_down_options']
    input_data['default_groups'] = default_groups_hold
    
    return [input_data]

@app.callback(
        Output("download", "data"),
        [
            Input("btn", "n_clicks"),
        ],
        [
            State('heatmap_memory', 'data'),
        ]
)
def generate_csv(n_nlicks, data):
    if(data is None):
        raise PreventUpdate
    # Convert data to a string.
    s = io.StringIO()
    d_df = pd.DataFrame(data['original'])
    d_df['Gene_Name'] = data['current_row_labels']
    rename_dict = {}
    for val in d_df.columns:
        if(val.endswith('_X') or val.endswith('_Y')):
            rename_dict[val] = val[:len(val)-2]
    #print(rename_dict)
    d_df = d_df.rename(rename_dict, axis=1)
    column_names_download = list(d_df.columns)

    d_df = d_df[['Gene_Name'] + list(column_names_download[:len(column_names_download) - 1])]
    d_df.to_csv(s, index=False)
    content = s.getvalue()
    # The output must follow this form for the download to work.
    return dict(filename="Heatmap_data.csv", content=content, type="text/csv")    



   

#Helper Functions

def filter_by_var_fold(main_data, x_dd_value, y_dd_value, fold_filter, var_filter):
    x_y_cols = list(set(x_dd_value + y_dd_value))
    main_df = pd.DataFrame(main_data['main_df'])
    new_row_labels = np.array(main_data['row_labels'])


    x_dff = main_df.loc[:, x_dd_value]
    y_dff = main_df.loc[:, y_dd_value]


    var_df = pd.concat([x_dff.reset_index(drop=True), y_dff.reset_index(drop=True)], axis=1)
    var_filter = var_df.var(axis=1) > var_filter 

    main_df = main_df[var_filter]
    x_dff = main_df.loc[:, x_dd_value]
    y_dff = main_df.loc[:, y_dd_value]
    
    new_row_labels = new_row_labels[var_filter]

    x_mean = np.array(x_dff.mean(axis=1))
    y_mean = np.array(y_dff.mean(axis=1))
    x_mean_log = np.log2(x_mean)
    y_mean_log = np.log2(y_mean)

    fold_change = y_mean_log - x_mean_log

    final_fold_filter = np.greater(fold_change, fold_filter) | np.less(fold_change, -fold_filter)

    main_df = main_df[final_fold_filter]
    new_row_labels = new_row_labels[final_fold_filter]
    main_df = main_df.reset_index()

    return load_default_data(fold_filter, var_filter, data = main_df,  row_labels_hold = new_row_labels)

def add_default_groups():
        default_groups_hold = {
        'CONTROL':['CTRL_0', 'CTRL_1', 'CTRL_2'],
        'P4':['P4_0', 'P4_1', 'P4_2'],
        'IL1B':['IL1B_0', 'IL1B_1', 'IL1B_2'],
        'FSK':['FSK_0', 'FSK_1', 'FSK_2'],
        'FSK_P4':['FSK_P4_0', 'FSK_P4_1', 'FSK_P4_2'],
        'P4_IL1B':['P4_IL1B_0', 'P4_IL1B_1', 'P4_IL1B_2'],
        'FSK_IL1B': ['FSK_IL1B_0','FSK_IL1B_1', 'FSK_IL1B_2',],
        'FSK_P4_IL1B': ['FSK_P4_IL1B_0', 'FSK_P4_IL1B_1','FSK_P4_IL1B_2']
    }

def create_new_data(input_df):
    '''
    Create new data object
    '''
    df_hold = input_df
    gene_name_col = list(df_hold.columns)[0]
    
    row_labels_hold = df_hold.loc[:, gene_name_col]
    row_labels_arr_hold = np.array(row_labels_hold)
    row_dict_hold = { row_labels_hold[i] : i for i in range(0, len(row_labels_hold) ) }
    if(gene_name_col in df_hold.columns):
        df_hold = df_hold.drop(gene_name_col, 1)
    df_row_means_hold = df_hold.mean(axis=1)
    df_row_stds_hold = df_hold.std(axis=1)
    df_row_vars_hold = df_hold.var(axis=1)

    default_options_hold = []
    column_options_hold = [{'label': i, 'value': i} for i in df_hold.columns]
    gene_drop_down_options_hold = default_options_hold + column_options_hold

    single_gene_options_hold = [{'label': i, 'value': i} for i in row_labels_hold]
    return_mem = {
    'main_df':df_hold.to_dict('list'),
    'df_row_means':df_row_means_hold,
    'df_row_stds':df_row_stds_hold,
    'df_row_vars': df_row_vars_hold,
    'row_labels':row_labels_hold,
    'row_dict':row_dict_hold,
    'default_groups':{},
    'default_options':default_options_hold,
    'column_options':column_options_hold,
    'gene_drop_down_options':gene_drop_down_options_hold,
    'single_gene_options':single_gene_options_hold,
    'fold_filter':fold_filter,
    'var_filter':var_filter,
    }
    return return_mem

def load_default_data(fold_filter, var_filter, data = None, row_labels_hold = None):
    '''
    Load defualt data objects
    '''
    if(data is None):
        df_hold = pd.read_csv(str(pathlib.Path("Gene_Expression_Data/GSE134896_FSKrepExpTable.csv")))
        row_labels_hold = list(df_hold['Gene_Name'])
    else:
        df_hold = data
        row_labels_hold = row_labels_hold
    
    row_labels_arr_hold = np.array(row_labels_hold)
    row_dict_hold = { row_labels_hold[i] : i for i in range(0, len(row_labels_hold) ) }

    if('Gene_Name' in df_hold.columns):
        df_hold = df_hold.drop('Gene_Name', 1)
    df_row_means_hold = df_hold.mean(axis=1)
    df_row_stds_hold = df_hold.std(axis=1)
    df_row_vars_hold = df_hold.var(axis=1)

    default_groups_hold = {
        'CONTROL':['CTRL_0', 'CTRL_1', 'CTRL_2'],
        'P4':['P4_0', 'P4_1', 'P4_2'],
        'IL1B':['IL1B_0', 'IL1B_1', 'IL1B_2'],
        'FSK':['FSK_0', 'FSK_1', 'FSK_2'],
        'FSK_P4':['FSK_P4_0', 'FSK_P4_1', 'FSK_P4_2'],
        'P4_IL1B':['P4_IL1B_0', 'P4_IL1B_1', 'P4_IL1B_2'],
        'FSK_IL1B': ['FSK_IL1B_0','FSK_IL1B_1', 'FSK_IL1B_2',],
        'FSK_P4_IL1B': ['FSK_P4_IL1B_0', 'FSK_P4_IL1B_1','FSK_P4_IL1B_2']
    }

    default_options_hold = [{'label': i, 'value': ' '.join(default_groups_hold[i])} for i in default_groups_hold.keys()]
    column_options_hold = [{'label': i, 'value': i} for i in df_hold.columns]
    gene_drop_down_options_hold = default_options_hold + column_options_hold

    single_gene_options_hold = [{'label': i, 'value': i} for i in row_labels_hold]
    return_mem = {
    'main_df':df_hold.to_dict('list'),
    'df_row_means':df_row_means_hold,
    'df_row_stds':df_row_stds_hold,
    'df_row_vars': df_row_vars_hold,
    'row_labels':row_labels_hold,
    'row_dict': row_dict_hold,
    'default_groups':default_groups_hold,
    'default_options':default_options_hold,
    'column_options':column_options_hold,
    'gene_drop_down_options':gene_drop_down_options_hold,
    'single_gene_options':single_gene_options_hold,
    'fold_filter':fold_filter,
    'var_filter':var_filter,
    'original_data':True,
    }
    return return_mem

def combine_columns(column_selection, main_data):
    '''
    Ensures that input x and y groups are lists of strings. Grabs the default column groups underlying columns and ensures no repeats within a selection.
    Sorts the columns alphabetically
    '''
    default_groups = main_data['default_groups']
    final_cols = []
    if(isinstance(column_selection, str)):
        if(default_groups.get(column_selection) != None):
            final_cols = default_groups.get(column_selection)
        else:
            final_cols = [column_selection]
    elif(isinstance(column_selection, list)):
        for i in range(len(column_selection)):
            current_column = column_selection[i]
            if(isinstance(current_column, str)):
                if(default_groups.get(current_column) != None):
                    final_cols = final_cols + default_groups.get(current_column)
                else:
                    current_column = current_column.split()
                    final_cols = final_cols + current_column
            else:
                raise PreventUpdate
    final_cols = list(set(final_cols))
    final_cols.sort()
    return final_cols

def normalize_gene(value, min_value = -2.0, max_value = 4.81852346, log_transform = True):
    if(log_transform):
        value = math.log10(value)

    return (value - min_value)/(max_value - min_value)

def find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value, main_data, quick_method = False):
    main_df = pd.DataFrame(main_data['main_df'])
    if(combination_dd_value == 'mean'):
        x_combo = np.array(x_dff.mean(axis=1))
        y_combo = np.array(y_dff.mean(axis=1))
    else:
        x_combo = np.array(x_dff.median(axis=1))
        y_combo = np.array(y_dff.median(axis=1))
    if(len(x_dff.columns) == 1 or len(y_dff.columns) == 1 or quick_method):
        difference_measure = np.subtract(x_combo, y_combo)
        absolute_difference_measure = np.absolute(difference_measure)
    else:
        difference_measure = np.array(main_df.apply(row_t_test, args = [x_dff.columns, y_dff.columns], axis = 1))
        absolute_difference_measure = np.absolute(difference_measure)


    gene_diff_df = pd.DataFrame({
        'original_index': range(0, len(difference_measure)),
        'difference measure' : difference_measure,
        'absolute difference measure' : absolute_difference_measure,
        'x_vals':x_combo,
        'y_vals' : y_combo,
        'gene name' : main_data['row_labels'],
        'means' : main_data['df_row_means'],
        'stds' : main_data['df_row_stds'],
        })

    gene_diff_df.sort_values(by = ["absolute difference measure"], inplace=True, ascending=False)

    gene_diff_df = gene_diff_df.reset_index()

    return gene_diff_df

def row_t_test(row, x_column_list, y_column_list):
    '''
    Performs a t-test for a row of data of x group columns vs y group columns
    '''
    x_distribution = list(row[x_column_list])
    y_distribution = list(row[y_column_list])

    t_statistic, p_value = stats.ttest_ind(x_distribution, y_distribution)

    return abs(t_statistic)

def create_heatmap_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data):
    '''
    Uses the gene_diff_df to create a heat df and normalized df that will be used to create heatmap
    '''
    #print("we updating heatmap mem")
    main_df = pd.DataFrame(main_data['main_df'])
    new_row_labels = main_data['row_labels']
    new_row_dict = main_data['row_dict']
    if(individ_or_gene_dd == 'Gene Focused'):
        if(type(gene_select_dd) == str):
            gene_select_dd = [gene_select_dd]
        gene_select_indices = []
        for val in gene_select_dd:
            try:
                gene_select_indices.append(new_row_dict[val])
            except:
                pass
        heat_df = main_df.iloc[main_df.index.isin(gene_select_indices)]

    else:
        gene_select_indices = list(gene_diff_df['original_index'])[:HEAT_GENES]
        heat_df = main_df.iloc[list(gene_diff_df['original_index'])[:HEAT_GENES]]
    
    current_row_labels = [new_row_labels[index] for index in gene_select_indices]

    name_adj_all_columns = [s + '_X' for s in x_dd_value] + [s + '_Y' for s in y_dd_value]

    #heat_means = heat_df.mean(axis=1).values
    #heat_stds = heat_df.std(axis=1).values

    narrowed_df = heat_df.loc[:, list(set(x_dd_value + y_dd_value))]
    heat_means = narrowed_df.mean(axis=1).values
    heat_stds = narrowed_df.std(axis=1).values

    x_df, x_df_normal = normalize_and_select_heat_df(heat_df, x_dd_value, heat_means, heat_stds)
    y_df, y_df_normal = normalize_and_select_heat_df(heat_df, y_dd_value, heat_means, heat_stds)

    name_adj_heat_df = combine_and_rename_heat_dfs(x_df, y_df)
    name_adj_normal_heat_df = combine_and_rename_heat_dfs(x_df_normal, y_df_normal)

    name_adj_heat_df = name_adj_heat_df[name_adj_all_columns]
    name_adj_normal_heat_df = name_adj_normal_heat_df[name_adj_all_columns]

    heat_dict = {
        'original' : name_adj_heat_df.to_dict('list'),
        'normal' : name_adj_normal_heat_df.to_dict('list'),
        'current_row_labels': current_row_labels,
        'x_dd' : x_dd_value,
        'y_dd' : y_dd_value,

    }
    #print("finsihed update heatmap mem")

    return heat_dict

def create_scatter_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data):
    '''
    Creates a dictionary that will be used to create the scatter plot memory object
    '''
    main_df = pd.DataFrame(main_data['main_df'])
    new_row_dict = main_data['row_dict']
    if(individ_or_gene_dd == 'Condition Focused'):
        scatter_dict = {
        'df1' : gene_diff_df.to_dict('list'),
        'df2' : [x_dd_value, y_dd_value],
        }
    else:
        if(isinstance(gene_select_dd, str)):
            gene_select_dd = [gene_select_dd]


        selected_indices_hold = []
        for val in gene_select_dd:
            try:
                selected_indices_hold.append(new_row_dict[val])
            except:
                pass


        scatter_dict = {
            'df1' : gene_diff_df.to_dict('list'),
            'names' : [x_dd_value, y_dd_value],
            'selected_indices':selected_indices_hold
        }

    
    return scatter_dict

def adjust_std_mean(heat_df, heat_means, heat_stds):
    '''
    Produces the z-score for every value based on its row
    '''
    new_df = {}
    for col in heat_df.columns:
        new_df[col] = np.divide(np.subtract(heat_df[col].values, heat_means), heat_stds)
    return pd.DataFrame(new_df)

def normalize_and_select_heat_df(heat_df, z_dd_value, heat_means, heat_stds):
    '''
    Selects the subset of data we are interested in and normalizes it
    '''
    z_df = heat_df[z_dd_value].reset_index(drop=True)

    z_df_normal = adjust_std_mean(z_df, heat_means, heat_stds).reset_index(drop=True)
    return z_df, z_df_normal

def combine_and_rename_heat_dfs(x_df, y_df):
    '''
    creates the final dataframe based off of the x and y heat dfs
    '''
    x_df = x_df.add_suffix('_X')
    y_df = y_df.add_suffix('_Y')

    heat_df = pd.concat([x_df, y_df], axis = 1)
    return heat_df

def get_hover_text(z_df, current_row_labels):
    '''
    Creates the displayed hover data for the heatmap
    '''
    hover_text = []
    z_df_columns = list(z_df.columns)
    for index, row in z_df.iterrows():
        hover_text.append(list())
        for col in z_df_columns:
            hover_text[-1].append('Individual: ' + col + '<br>' + 'Gene: ' + str(current_row_labels[index]) + '<br>' + 'Value: ' + str(row[col]))

    return hover_text

def create_scatter_diff_genes(gene_diff_df, x_dd_value, y_dd_value):
    '''
    Creates the scatter plot based upon the differentially expressed genes
    '''
    gene_diff_df.sort_values(by = ["absolute difference measure"], inplace=True, ascending=False)
    top_hundred_indices = list(gene_diff_df['original_index'])[:HEAT_GENES]
    bottom_gene = gene_diff_df[~gene_diff_df['original_index'].isin(top_hundred_indices)]
    bottom_gene_indices = gene_diff_df['original_index']


    final_x_values = list(gene_diff_df['x_vals'])
    final_y_values = list(gene_diff_df['y_vals'])
    final_names = list(gene_diff_df['gene name'])

    top_x_values= final_x_values[:HEAT_GENES]
    top_y_values= final_y_values[:HEAT_GENES]
    top_row_labels = final_names[:HEAT_GENES]


    bottom_x_values= final_x_values[HEAT_GENES:]
    bottom_y_values= final_y_values[HEAT_GENES:]
    bottom_row_labels = final_names[HEAT_GENES:]

    updated_scatter = go.Figure()

    updated_scatter.add_trace(go.Scattergl(
        x = bottom_x_values,
        y = bottom_y_values,
        text = bottom_row_labels,
        mode='markers',
        marker_color='blue',
        marker={
            'size':3
        },
        legendgroup='Differentially Expressed Genes',
        name='Non Differential Genes',

    ))
    
    updated_scatter.add_trace(go.Scattergl(
        x = top_x_values,
        y = top_y_values,
        text = top_row_labels,
        mode='markers',
        marker_color='red',
        marker={
            'size':7
        },
        legendgroup='Non Differentially Expressed Genes',
        name='Differential Genes',

    ))

    updated_scatter.update_layout(
        title={
            'text': "Scatter Plot",
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        autosize=True,
        xaxis=go.layout.XAxis(
            title=str(', '.join(x_dd_value)),
            showgrid=True,
            ),
        yaxis=go.layout.YAxis(
            title=str(', '.join(y_dd_value)),
            showgrid=True,  
            ),
        xaxis_type="log", 
        yaxis_type="log",
        showlegend=True

        )
    return updated_scatter

def create_gene_based_scatter(gene_diff_df, selected_indices, x_dd_value, y_dd_value, main_data):
    '''
    Creates the scatter plot based upon the differentially expressed genes
    '''
    new_row_dict = main_data['row_dict']



    final_x_values = list(gene_diff_df['x_vals'])
    final_y_values = list(gene_diff_df['y_vals'])
    final_names = list(gene_diff_df['gene name'])
    top_x_values = []
    top_y_values = []
    top_row_labels = []
    bottom_x_values= []
    bottom_y_values=  []
    bottom_row_labels = []

    for i, name in enumerate(final_names):
        if(new_row_dict[name] in selected_indices):
            top_x_values.append(final_x_values[i])
            top_y_values.append(final_y_values[i])
            top_row_labels.append(final_names[i])
        else:
            bottom_x_values.append(final_x_values[i])
            bottom_y_values.append(final_y_values[i])
            bottom_row_labels.append(final_names[i])


    updated_scatter = go.Figure()

    updated_scatter.add_trace(go.Scattergl(
        x = bottom_x_values,
        y = bottom_y_values,
        text = bottom_row_labels,
        mode='markers',
        marker_color='blue',
        marker={
            'size':3
        },
        legendgroup='Differentially Expressed Genes',
        name='Non Differential Genes',

    ))
    
    updated_scatter.add_trace(go.Scattergl(
        x = top_x_values,
        y = top_y_values,
        text = top_row_labels,
        mode='markers',
        marker_color='red',
        marker={
            'size':7
        },
        legendgroup='Non Differentially Expressed Genes',
        name='Differential Genes',

    ))

    updated_scatter.update_layout(
        title={
            'text': "Scatter Plot",
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        autosize=True,
        xaxis=go.layout.XAxis(
            title=str(', '.join(x_dd_value)),
            showgrid=True,
            ),
        yaxis=go.layout.YAxis(
            title=str(', '.join(y_dd_value)),
            showgrid=True,  
            ),
        xaxis_type="log", 
        yaxis_type="log",
        showlegend=True

        )
    return updated_scatter

def create_x_y_group_display_names(x_dd_value, y_dd_value):
    '''
    Creates x and y group displays
    '''
    x_names = " ".join(x_dd_value)
    y_names = " ".join(y_dd_value)
    final_x_display = [html.Br(), 'Current X group: ']
    final_x_display.append(x_names)
    final_y_display = [html.Br(), 'Current Y group: ']
    final_y_display.append(y_names)
    final_x_display.append(html.Br())
    final_y_display.append(html.Br())
    return final_x_display, final_y_display

def update_violin_plot(current_row_labels, x_dd_value, y_dd_value, violin_gene_number, main_df, row_dict_hold):
    '''
    Graphs the selected indices in a violin plot
    '''
    color_option_1 = ['green', 'orange', 'yellow', 'white', 'grey']
    color_option_2 = ['purple', 'black', 'pink', 'red', 'blue']
    violin_type = violin_gene_number
    updated_violin = go.Figure()
    ##print(violin_gene_number)
    reduced_x_df = main_df[x_dd_value]
    reduced_y_df = main_df[y_dd_value]
    x_names = " ".join(x_dd_value)
    y_names = " ".join(y_dd_value)

    violin_gene_number = 10
    all_figs = []
    show_legend = [True] + [False] * (violin_gene_number - 1)
    colors = ['green','purple']
    fig = go.Figure()
    if(violin_type == 'Scatter Plot'):
        for index, value  in enumerate(current_row_labels[:5]):
            if(index <= violin_gene_number):
                index_in_df = row_dict_hold[value]
                for choice, r_colors, selection_type in zip([reduced_x_df, reduced_y_df], [color_option_1, color_option_2], ['X Selection', 'Y Selection']):
                    sub_types = [name[:-2] for name in choice.columns]
                    reduced_sub_types = list(set(sub_types))
                    count = 0
                    for sub in reduced_sub_types:
                        length_of_sub = sub_types.count(sub)
                        cols = [x for i, x in enumerate(choice.columns) if x[:-2] == sub]
                        rr_df = choice[cols]
                        ##print(cols)
                        rrr_df = list(rr_df.iloc[index_in_df])

                        ##print("ahhhhh")
                        
                        fig.add_trace(go.Scattergl(
                            x=rrr_df,
                            y=[value] * len(rr_df),
                            name=cols[0][:-2],
                            showlegend = show_legend[index],
                            marker=dict(
                                color=r_colors[count],
                                line_color=r_colors[count],
                            ),
                            hovertemplate =
                                '<b>Expression</b>: %{x}'+
                                '<br><b>Gene</b>: %{y}<br>'+
                                '<b>Condition</b>  %{text}',
                            text = cols,
                        ))
                        count += 1

                fig.update_traces(mode='markers', marker=dict(line_width=3, symbol=42, size=16))
                fig.update_layout(xaxis_type="log")
    else:
        all_dfs = []
        counter = 0
        for index, value  in enumerate(current_row_labels[:3]):
            ##print(value)
            if(index < 3):
                counter +=1
                index_in_df = row_dict_hold[value]
                rr_x_df = list(reduced_x_df.iloc[index_in_df])
                rr_y_df = list(reduced_y_df.iloc[index_in_df])

                figure={
                    'data': [
                        {
                            'x': list(reduced_x_df.iloc[index_in_df]),
                            'name': 'X Group',
                            'type': 'histogram'
                        },
                        {
                            'x': list(reduced_y_df.iloc[index_in_df]),
                            'name': 'Y Group',
                            'type': 'histogram'
                        }
                    ],
                    'layout': {
                        'barmode':'stack', 
                        'title_text':value
                    }
                }

                grouping = []
                values = []
                hover_info = []
                colors = []
                for count, col_name in enumerate(reduced_x_df.columns):
                    colors.append('green')
                    hover_info.append(col_name)
                    grouping.append('X Group')
                    values.append(rr_x_df[count])
                for count, col_name in enumerate(reduced_y_df.columns):
                    colors.append('purple')
                    hover_info.append(col_name)
                    grouping.append('Y Group')
                    values.append(rr_y_df[count])

                new_df = pd.DataFrame({
                    'value':values,
                    'Group':grouping,
                    'hover_info':hover_info,
                    'color':colors,
                    'Gene Name':[value] * len(colors)

                    })
                all_dfs.append(new_df)

                new_colors = {
                    'X Group':'green',
                    'Y Group':'purple'
                }
        final_df = pd.concat(all_dfs)
        #print(len(final_df))
        count_list = list(final_df['value'])
        start_val = min(count_list)
        end_val = max(count_list)
        bins = np.logspace(start_val - 1, end_val + 1, 15)
        fig = px.histogram(final_df, x="value", color="Group", marginal="rug", # can be `box`, `violin`
                                barmode="overlay", opacity=0.75, hover_name="hover_info", color_discrete_map = new_colors, facet_col = 'Gene Name', range_x=[start_val, end_val])
        #fig.update_layout(xaxis_type="log")
        #for count in range(counter):
        #    fig['layout'][f'xaxis{str(counter + 1)}'].update(tickvals = bins)
        #fig.update_xaxes(matches=None)
        #for count in range(counter):
        #    fig.update_xaxes(type="log", row = 1, col = count + 1)







    

    return fig

def removeSuffixs(sentence):

    suffixList = ["_X", "_Y"] #add more as nessecary

    for item in suffixList:
        if item in sentence:

            sentence = sentence.replace(item, "")
            repeatLetters = next((True for char, group in groupby(sentence)
                                  if sum(1 for _ in group) >= 2), False)

            if repeatLetters:

                sentence = sentence[:-1]

    return sentence 

'''
        fig = px.histogram(new_df, x="value", color="Group", marginal="rug", # can be `box`, `violin`
                        barmode="overlay", nbins = 15, opacity=0.75, title = value, hover_name="hover_info", color_discrete_map = new_colors)
        if(index == 1):
            fig = px.histogram(new_df, x="value", color="Group", marginal="box", # can be `box`, `violin`
                        barmode="overlay", nbins = 15, opacity=0.75, title = value, hover_name="hover_info", color_discrete_map = new_colors)
    
    for i, name in enumerate(final_names):
        if(new_row_dict[name] in selected_indices):
            top_x_values.append(final_x_values[i])
            top_y_values.append(final_y_values[i])
            top_row_labels.append(final_names[i])
        else:
            bottom_x_values.append(final_x_values[i])
            bottom_y_values.append(final_y_values[i])
            bottom_row_labels.append(final_names[i])


    updated_scatter = go.Figure()

    updated_scatter.add_trace(go.Scatter(
        x = bottom_x_values,
        y = bottom_y_values,
        text = bottom_row_labels,
        mode='markers',
        marker_color='blue',
        marker={
            'size':3
        },
        legendgroup='Differentially Expressed Genes',
        name='Non Differential Genes',

    ))
    
    updated_scatter.add_trace(go.Scatter(
        x = top_x_values,
        y = top_y_values,
        text = top_row_labels,
        mode='markers',
        marker_color='red',
        marker={
            'size':7
        },
        legendgroup='Non Differentially Expressed Genes',
        name='Differential Genes',

    ))

    updated_scatter.update_layout(
        title={
            'text': "Scatter Plot",
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        autosize=True,
        xaxis=go.layout.XAxis(
            title=str(', '.join(x_dd_value)),
            showgrid=True,
            ),
        yaxis=go.layout.YAxis(
            title=str(', '.join(y_dd_value)),
            showgrid=True,  
            ),
        xaxis_type="log", 
        yaxis_type="log",
        showlegend=True

        )
    return updated_scatter
@app.callback(
        [
            Output(component_id='scatter_plot', component_property='figure'),
            Output(component_id='violin_differential', component_property='figure'),
            Output(component_id='heatmap', component_property='figure'),
            Output(component_id='x_group_id', component_property='children'),
            Output(component_id='y_group_id', component_property='children')
        ],
        [
            Input('update_button', 'n_clicks')
        ],
        state=[
            State('x_axis_selection_dd', 'value'),
            State('y_axis_selection_dd', 'value'),
            State('combination_selection_dd', 'value'),
            State('invidual_or_gene_dd', 'value'),
            State('violin_gene_number', 'value'),
            State('gene_select_dd', 'value'),
            State('scatter_plot', 'figure'),
            State('violin_differential', 'figure'),
        ]
)
def update_button(n_clicks, x_dd_value, y_dd_value, combination_dd_value, individ_or_gene_dd, violin_gene_number, gene_select_dd, scatter, violin):
    x_dd_value = combine_columns(x_dd_value)
    y_dd_value = combine_columns(y_dd_value)
    #print(x_dd_value)
    #print(y_dd_value)
    x_names = " ".join(x_dd_value)
    y_names = " ".join(y_dd_value)

    valid_selection = determine_combination_type(x_dd_value, y_dd_value)
    updated_scatter = scatter
    updated_violin = violin
    updated_heatmap = heatmap_graph

    x_names = " ".join(x_dd_value)
    y_names = " ".join(y_dd_value)
    final_x_display = [html.Br(), 'Current X group: ']
    final_x_display.append(x_names)
    final_y_display = [html.Br(), 'Current Y group: ']
    final_y_display.append(y_names)
    final_x_display.append(html.Br())
    final_y_display.append(html.Br())

    if(individ_or_gene_dd == 'Gene Focused'):
        updated_heatmap = gene_based_analysis(gene_select_dd, x_dd_value, y_dd_value)
        if(len(x_dd_value) > 1 and len(y_dd_value) > 1):
            updated_scatter = pca_by_gene(x_dd_value, y_dd_value)
    else:
        if(valid_selection == 'One Missing'):
            raise PreventUpdate
        else:
            reduced_x_df = df.loc[:, x_dd_value]
            reduced_y_df = df.loc[:, y_dd_value]
            final_x_values, final_y_values = consolidate_data(reduced_x_df, reduced_y_df, combination_dd_value)
            if(valid_selection == 'Singles'):
                df_combined_vals = pd.DataFrame({
                    'final_x_values': final_x_values,
                    'final_y_values': final_y_values,
                })
                difference_measure = np.array(df_combined_vals.apply(differentially_expressed_genes_difference, axis = 1))
                top_hundred_genes = np.argsort(-difference_measure)[:100]
                top_graphed_genes = top_hundred_genes[:violin_gene_number]
            elif(valid_selection == 'Groups'):
                file_path, cache_exists, top_hundred_genes, difference_measure = check_difference_measure_cache(cache, x_dd_value, y_dd_value)
                if(cache_exists):
                    #print('Used Cache')
                    top_graphed_genes = top_hundred_genes[:violin_gene_number]
                else:
                    difference_measure = np.array(df.apply(differentially_expressed_genes_t_test, args = [x_dd_value, y_dd_value], axis = 1))
                    cache_df = pd.DataFrame({
                        'indices': range(0, len(difference_measure)),
                        'difference_measure':difference_measure
                    })
                    cache_df = cache_df.sort_values(by=['difference_measure'], ascending=False)

                    top_hundred_genes = list(cache_df['indices'])[:100]
                    top_graphed_genes = top_hundred_genes[:violin_gene_number]

                    cache_df.to_csv(file_path)


            updated_violin = update_violin_plot(top_graphed_genes, reduced_x_df, reduced_y_df, x_names, y_names, violin_gene_number)

            updated_scatter = update_scatter_plot(list(final_x_values), list(final_y_values), top_hundred_genes, x_dd_value, y_dd_value)

            updated_heatmap = update_heatmap(top_hundred_genes, x_dd_value, y_dd_value)



    #print('Finished Updating')
    return updated_scatter, updated_violin, updated_heatmap, final_x_display, final_y_display 
@app.callback(
        [
            Output(component_id='scatter_plot', component_property='figure'),
        ],
        [
            Input("scatter_plot", "selectedData")
        ],
        state=[
            State('x_axis_selection_dd', 'value'),
            State('y_axis_selection_dd', 'value'),
            State('invidual_or_gene_dd', 'value'),
            State('scatter_plot', 'figure'),
            

        ]
)
def make_individual_figure(main_scatter_hover, x_dd_value, y_dd_value, individ_or_gene, scatter_plot):
    #print(main_scatter_hover)
    return scatter_plot
def gene_based_analysis(gene_select_dd, x_dd_value, y_dd_value):

    if(type(gene_select_dd) == str):
        gene_select_dd = [gene_select_dd]
    gene_select_indices = []
    for val in gene_select_dd:
        gene_select_indices.append(row_dict[val])


    heat_df = df.iloc[df.index.isin(gene_select_indices)]

    all_columns = x_dd_value + y_dd_value
    name_adj_all_columns = [s + '_X' for s in x_dd_value] + [s + '_Y' for s in y_dd_value]

    heat_means = df_row_means[gene_select_indices].values
    heat_stds = df_row_stds[gene_select_indices].values

    x_df, x_df_normal = heatmap_naming_adjustment(heat_df, x_dd_value, heat_means, heat_stds)
    y_df, y_df_normal = heatmap_naming_adjustment(heat_df, y_dd_value, heat_means, heat_stds)

    name_adj_heat_df = combine_heat_dfs(x_df, y_df)
    name_adj_normal_heat_df = combine_heat_dfs(x_df_normal, y_df_normal)


    name_adj_heat_df = name_adj_heat_df[name_adj_all_columns]
    name_adj_normal_heat_df = name_adj_normal_heat_df[name_adj_all_columns]


    heatmap_Data = go.Heatmap(
            z=name_adj_normal_heat_df,
            x=all_columns,
            y=gene_select_dd, 
            hoverinfo='text',
            text=get_hover_text(name_adj_heat_df, gene_select_dd),
            colorscale='RdBu'
        )
    layout_heatmap = go.Layout(
            title='Heatmap',
            autosize=True,
            height=500,
            xaxis=go.layout.XAxis(
                showgrid=True,
            ),
            yaxis=go.layout.YAxis(
                showgrid=True,
            )
        )
    fig_heatmap = go.Figure(data=[heatmap_Data], layout=layout_heatmap)
    return fig_heatmap
def update_heatmap(top_hundred_genes, x_dd_value, y_dd_value):
    heat_df = df.iloc[top_hundred_genes]

    all_columns = x_dd_value + y_dd_value
    name_adj_all_columns = [s + '_X' for s in x_dd_value] + [s + '_Y' for s in y_dd_value]

    heat_means = df_row_means[top_hundred_genes].values
    heat_stds = df_row_stds[top_hundred_genes].values

    x_df, x_df_normal = heatmap_naming_adjustment(heat_df, x_dd_value, heat_means, heat_stds)
    y_df, y_df_normal = heatmap_naming_adjustment(heat_df, y_dd_value, heat_means, heat_stds)

    name_adj_heat_df = combine_heat_dfs(x_df, y_df)
    name_adj_normal_heat_df = combine_heat_dfs(x_df_normal, y_df_normal)


    name_adj_heat_df = name_adj_heat_df[name_adj_all_columns]
    name_adj_normal_heat_df = name_adj_normal_heat_df[name_adj_all_columns]




    current_row_labels = [row_labels[index] for index in top_hundred_genes]
    heatmap_Data = go.Heatmap(
            z=name_adj_normal_heat_df,
            x=name_adj_all_columns,
            y=current_row_labels, 
            hoverinfo='text',
            text=get_hover_text(name_adj_heat_df, current_row_labels),
            colorscale='RdBu'
        )
    layout_heatmap = go.Layout(
            title='Heatmap',
            autosize=True,
            height=500,
            xaxis=go.layout.XAxis(
                showgrid=True,
            ),
            yaxis=go.layout.YAxis(
                showgrid=True,
            )
        )
    fig_heatmap = go.Figure(data=[heatmap_Data], layout=layout_heatmap)
    return fig_heatmap
def update_violin_plot(top_graphed_genes, reduced_x_df, reduced_y_df, x_names, y_names, violin_gene_number):
    updated_violin = go.Figure()

    show_legend = [True] + [False] * (violin_gene_number - 1)
    counter = 0

    for index_value in top_graphed_genes:

        updated_violin.add_trace(go.Violin(
                    x=[row_labels[index_value] for i in range(len(reduced_x_df.columns))],
                    y=list(reduced_x_df.iloc[index_value]),
                    side='negative',
                    legendgroup=x_names,
                    name=x_names,
                    line_color = 'green',
                    showlegend = show_legend[counter]
                )
            )
        updated_violin.add_trace(go.Violin(
                    x=[row_labels[index_value] for i in range(len(reduced_y_df.columns))],
                    y=list(reduced_y_df.iloc[index_value]),
                    side='positive',
                    legendgroup = y_names,
                    name =y_names,
                    line_color = 'purple',
                    showlegend = show_legend[counter]


                )
            )
        counter = counter + 1

    updated_violin.update_layout(

            showlegend=True,
            autosize=True,
            violingap=0, 
            violinmode='overlay',
            yaxis_type="log",

        )
    updated_violin.update_traces(points = 'all')
    return updated_violin
def determine_combination_type(x_dd_value, y_dd_value):
    valid_selection = 'One Missing'
    if(len(x_dd_value) == 1 or len(y_dd_value) == 1):
        valid_selection = 'Singles'
    elif(len(x_dd_value) > 1 or len(y_dd_value) > 1):
        valid_selection = 'Groups'
    return valid_selection
def differentially_expressed_genes_difference(row):

    difference = abs(row['final_x_values'] - row['final_y_values']) 

    return difference
def check_difference_measure_cache(cache, x_dd_value, y_dd_value):
    x_dd_value.sort()
    y_dd_value.sort()
    x_dd_value = ''.join(x_dd_value)
    y_dd_value = ''.join(y_dd_value)
    cols = [x_dd_value, y_dd_value]
    cols.sort()
    current_key = 'sep'.join(cols)
    file_path = str(pathlib.Path('Preterm Birth Gene Data/t_test_cache/' + current_key + '.csv'))
    if(1==2 and cache.get(current_key) is not None):
        #print(current_key)
        cache_df = cache[current_key]
        indices = np.array(cache_df['indices'])
        difference_measure = np.array(cache_df['difference_measure'])
        path_exists = True

    else:
        path_exists = False
        indices = None
        difference_measure = None

    return file_path, path_exists, indices, difference_measure
'''

app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

if __name__ == '__main__':

    app.run_server(debug=False)
