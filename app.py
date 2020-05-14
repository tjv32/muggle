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
import os
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import math
import io
import base64

NUMBER_OF_GENES_IN_VIOLIN_PLOT = 5

df = pd.read_csv(str(pathlib.Path("Preterm Birth Gene Data/GSE134896_FSKrepExpTable.csv")))


def load_t_test_cache():
    file_dir = pathlib.Path('Preterm Birth Gene Data/t_test_cache/')
    all_files = os.listdir(file_dir)
    cache = {}
    for i in range(len(all_files)):
        hold = pd.read_csv(file_dir / all_files[i])
        #print(all_files[i])
        cache[all_files[i][:-4]] = hold

    return cache

#cache = load_t_test_cache()

row_labels = list(df['Gene_Name'])
row_labels_arr = np.array(row_labels)
row_dict = { row_labels[i] : i for i in range(0, len(row_labels) ) }


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

fig_scatter = go.Figure(
        data = [
            go.Scatter(
                    x = list(df.loc[:, df.columns[0]]),
                    y = list(df.loc[:, df.columns[1]]),
                    mode='markers',
                    marker={
                        'size':5
                    },
                )
            ],
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
violin_graph = dcc.Graph(
        id='violin_plot',
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
        children = 'Gene or Individual Comparison',
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
        options=[{'label': i, 'value': i} for i in ['Gene', 'Individual']],
        value='Individual',
        multi=False,
    )
number_of_genes_label = html.Label(
        children = 'Number of Genes Violin Plot Displays',
        style={
            'textAlign': 'left',
            'color': colors['text'],
        }
    )
number_of_violin_genes = dcc.Dropdown(
        id='violin_gene_number',
        options=[{'label': i, 'value': i} for i in range(1, 11)],
        value=5,
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
                        ],
                        width = 4
                    ),
                    dbc.Col(
                        [
                            violin_graph,
                        ], 
                        width = 8
                    )
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
app.layout = html.Div([
    dcc.Tabs([
        dcc.Tab(label='Visuals', children=[
            first_tab
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
            State('var_filter', 'value')

        ]
)
def update_button(n_clicks, main_data, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, violin_gene_number, combination_dd_value, name_memory, scatter_mem, fold_filter, var_filter):
    '''
        This callback takes all the dropdown menus and will calculate the differentially expressed genes and update the memory objects
    '''
    if(1==1):
        if(x_dd_value is None or y_dd_value is None):
            raise PreventUpdate
        x_num = len(x_dd_value)
        y_num = len(y_dd_value)

        if(x_num == 0 or y_num == 0 or main_data is None):
            raise PreventUpdate

        x_dd_value = combine_columns(x_dd_value, main_data)
        y_dd_value = combine_columns(y_dd_value, main_data)
        x_num = len(x_dd_value)
        y_num = len(y_dd_value)

        print('hi')


        

        orig_fold = main_data['fold_filter']
        orig_var = main_data['var_filter']



        main_df = pd.DataFrame(main_data['main_df'])

        if(len(main_df) == 0 or not (set(x_dd_value).issubset(set(main_df.columns))) or not (set(y_dd_value).issubset(set(main_df.columns)))):
            raise PreventUpdate

        main_data = filter_by_var_fold(main_data, x_dd_value, y_dd_value, float(fold_filter), float(var_filter))
        main_df = pd.DataFrame(main_data['main_df'])

        x_dff = main_df.loc[:, x_dd_value]
        y_dff = main_df.loc[:, y_dd_value]
        print('preprocessing done')
        if(individ_or_gene_dd == 'Individual'):

            if(name_memory is not None and x_dd_value == name_memory['x_dd'] and y_dd_value == name_memory['y_dd']
                and combination_dd_value == name_memory['combo'] and individ_or_gene_dd == name_memory['type']
                and violin_gene_number == name_memory['violin_gene_number'] and abs(orig_fold - fold_filter) < 0.001
                and abs(orig_var - var_filter) < 0.001):

                gene_diff_df = pd.DataFrame(scatter_mem['df1'])

            else:
                print('this should happen')
                gene_diff_df = find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value, main_data)
                print('this happened')

        else:

            gene_diff_df = None
        
        heat_dict = create_heatmap_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data)

        scatter_dict = create_scatter_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data)

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
        ],
        [
            Input('temp_data', 'data'),
        ],
        state=[
            State('update_button', 'n_clicks')
        ]
        
)
def update_main_memory(temp_data, n_clicks):
    if(n_clicks is None):
        n_clicks = 0
    if(temp_data is None):
        main_data = load_default_data(1,1)
    else:
        main_data = temp_data

    return main_data, n_clicks+1, main_data['gene_drop_down_options'], main_data['gene_drop_down_options'], main_data['single_gene_options']

@app.callback(
        
            Output(component_id='heatmap', component_property='figure'),
        [
            Input('heatmap_memory', 'data'),
            Input("scatter_plot", "selectedData"),
        ],

)
def update_heatmap(memory, selected_data):
    '''
        This callback uses the heat map memory store and any selected data in the scatter plot to update the heatmap
    '''
    if(memory is None):
        raise PreventUpdate
    name_adj_heat_df = pd.DataFrame(memory['original'])
    name_adj_normal_heat_df = pd.DataFrame(memory['normal'])

    current_row_labels = memory['current_row_labels']

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
            colorscale='RdBu',
            zmax=max(max_val + 0.2, 0.5),
            zmin= min(min_val - 0.2, -0.5),
        )
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
    fig_heatmap = go.Figure(data=[heatmap_Data], layout=layout_heatmap)

    return fig_heatmap


@app.callback(

            Output(component_id='scatter_plot', component_property='figure'),

        [
            Input('scatter_plot_memory', 'data'),
        ],
)
def update_scatter(memory):
    '''
        This callback uses the scatter plot memory store and any selected data in the heatmap to update the scatter plot
    '''
    if(memory is None):
        raise PreventUpdate

    stored_mem1 = memory['df1']
    stored_mem2 = memory.get('df2')

    if(isinstance(stored_mem2, list)):
        gene_diff_df = pd.DataFrame(stored_mem1)
        gene_diff_df.sort_values(by = ["absolute difference measure"], inplace=True, ascending=False)
        scatter_figure = create_scatter_diff_genes(gene_diff_df, stored_mem2[0], stored_mem2[1])

    else:
        stored_mem3 = memory['names']
        scatter_figure = create_pca_scatter(stored_mem1, memory['selected_indices'], stored_mem3[0], stored_mem3[1])

    return scatter_figure

@app.callback(

            Output(component_id='violin_plot', component_property='figure'),

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

    main_df = pd.DataFrame((main_data['main_df']))
    print(pd.DataFrame((main_data['main_df'])).columns)
    print(memory['x_dd'])

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


    current_row_labels = current_row_labels[:violin_gene_number]

    violin_figure = update_violin_plot(current_row_labels, name_memory['x_dd'], name_memory['y_dd'], violin_gene_number, main_df, main_data['row_dict'])

    

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
    print('AHHHHHHHHH')
    #print(contents)
    print(contents is None)
    if(contents is None):
        raise PreventUpdate

    if(1==1):
        print('made it here')
        content_type, content_string = contents[0].split(',')
        decoded = base64.b64decode(content_string)
        input_df = pd.read_csv(
            io.StringIO(decoded.decode('utf-8')))
        return_msg = html.Div([
            f'{filename[0]} successfully uploaded'
        ])
        print('up to this functions')
        return_mem = create_new_data(input_df)
        print('and past it')
    else:# Exception as e:
        return_msg = html.Div([
            'There was an error processing this file.'
        ])
        return_mem = None
    print('made it through this')

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

    print(input_data['gene_drop_down_options'])
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
        raise PreventUpdate

    gene_drop_down_options_hold = []
    default_groups_hold = {}
    for name, dd in zip(input_names, dropdowns):

        if(name['props'].get('value') is not None and dd['props'].get('value') is not None):
            print('we adding')
            gene_drop_down_options_hold.append({'label': name['props']['value'], 'value' : ' '.join(list(dd['props']['value']))})
            default_groups_hold[name['props']['value']] = [dd['props']['value']]
    #gene_drop_down_options_hold = [{'label': name['value'], 'value' : ' '.join(list(dd['value']))} for name, dd in zip(input_names, dropdowns)]
    #default_groups_hold = {name['value'] : [dd['value']] for name, dd in zip(input_names, dropdowns)}
    print(gene_drop_down_options_hold)
    input_data['gene_drop_down_options'] = gene_drop_down_options_hold + input_data['gene_drop_down_options']
    input_data['default_groups'] = default_groups_hold
    
    return [input_data]

        

#Helper Functions

def filter_by_var_fold(main_data, x_dd_value, y_dd_value, fold_filter, var_filter):

    print(fold_filter)
    print(var_filter)

    x_y_cols = list(set(x_dd_value + y_dd_value))
    main_df = pd.DataFrame(main_data['main_df'])
    new_row_labels = np.array(main_data['row_labels'])

    print('starting')
    x_dff = main_df.loc[:, x_dd_value]
    y_dff = main_df.loc[:, y_dd_value]
    print('quarter')

    var_df = pd.concat([x_dff.reset_index(drop=True), y_dff.reset_index(drop=True)], axis=1)
    var_filter = var_df.var(axis=1) > var_filter 

    main_df = main_df[var_filter]
    x_dff = main_df.loc[:, x_dd_value]
    y_dff = main_df.loc[:, y_dd_value]
    
    new_row_labels = new_row_labels[var_filter]
    print('half way')
    x_mean = np.array(x_dff.mean(axis=1))
    y_mean = np.array(y_dff.mean(axis=1))
    x_mean_log = np.log2(x_mean)
    y_mean_log = np.log2(y_mean)

    fold_change = y_mean_log - x_mean_log
    print(fold_change)
    final_fold_filter = np.greater(fold_change, fold_filter) | np.less(fold_change, -fold_filter)

    main_df = main_df[final_fold_filter]
    new_row_labels = new_row_labels[final_fold_filter]
    main_df = main_df.reset_index()
    print('we did it!')
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
    print('just beginning function')
    df_hold = input_df
    gene_name_col = list(df_hold.columns)[0]
    
    row_labels_hold = df_hold.loc[:, gene_name_col]
    row_labels_arr_hold = np.array(row_labels_hold)
    row_dict_hold = { row_labels_hold[i] : i for i in range(0, len(row_labels_hold) ) }
    print('quarter way thru')
    if(gene_name_col in df_hold.columns):
        df_hold = df_hold.drop(gene_name_col, 1)
    df_row_means_hold = df_hold.mean(axis=1)
    df_row_stds_hold = df_hold.std(axis=1)
    df_row_vars_hold = df_hold.var(axis=1)

    print('halfway thru function')
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
        df_hold = pd.read_csv(str(pathlib.Path("Preterm Birth Gene Data/GSE134896_FSKrepExpTable.csv")))
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

def find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value, main_data):
    main_df = pd.DataFrame(main_data['main_df'])
    if(combination_dd_value == 'mean'):
        x_combo = np.array(x_dff.mean(axis=1))
        y_combo = np.array(y_dff.mean(axis=1))
    else:
        x_combo = np.array(x_dff.median(axis=1))
        y_combo = np.array(y_dff.median(axis=1))
    if(len(x_dff.columns) == 1 or len(y_dff.columns) == 1):
        difference_measure = np.subtract(x_combo, y_combo)
        absolute_difference_measure = np.absolute(difference_measure)
    else:
        difference_measure = np.array(main_df.apply(row_t_test, args = [x_dff.columns, y_dff.columns], axis = 1))
        absolute_difference_measure = np.absolute(difference_measure)

    print(len(difference_measure))
    print(len(absolute_difference_measure))
    print(len(main_data['df_row_means']))
    print(len(x_combo))
    print(len(main_data['df_row_stds']))
    print(len(main_data['row_labels']))
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
    main_df = pd.DataFrame(main_data['main_df'])
    new_row_labels = main_data['row_labels']
    if(gene_diff_df is None):
        if(type(gene_select_dd) == str):
            gene_select_dd = [gene_select_dd]
        gene_select_indices = []
        for val in gene_select_dd:
            gene_select_indices.append(row_dict[val])
        heat_df = main_df.iloc[main_df.index.isin(gene_select_indices)]

    else:
        gene_select_indices = list(gene_diff_df['original_index'])[:100]
        heat_df = main_df.iloc[list(gene_diff_df['original_index'])[:100]]
    
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

    return heat_dict

def create_scatter_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, main_data):
    '''
    Creates a dictionary that will be used to create the scatter plot memory object
    '''
    if(individ_or_gene_dd == 'Individual'):
        scatter_dict = {
        'df1' : gene_diff_df.to_dict('list'),
        'df2' : [x_dd_value, y_dd_value],
        }
    else:
        if(isinstance(gene_select_dd, str)):
            gene_select_dd = [gene_select_dd]


        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(df[x_dd_value + y_dd_value])
        principal_df = pd.DataFrame(data = principal_components,
        columns = ['principal component 1', 'principal component 2'])

        
        '''

        if(len(y_dd_value) == 1):
            principal_df_y = df[y_dd_value]
            principal_df_y.columns = ['principal component 1']
        else:
            pca_y = PCA(n_components=2)
            principal_components_y = pca_y.fit_transform(df[y_dd_value])
            principal_df_y = pd.DataFrame(data = principal_components_y
                     , columns = ['principal component 1', 'principal component 2'])
        '''


        scatter_dict = {
            'df1' : principal_df.to_dict('list'),
            'names' : [x_dd_value, y_dd_value],
            'selected_indices':[row_dict[i] for i in gene_select_dd]
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
    top_hundred_indices = list(gene_diff_df['original_index'])[:100]
    bottom_gene = gene_diff_df[~gene_diff_df['original_index'].isin(top_hundred_indices)]
    bottom_gene_indices = gene_diff_df['original_index']


    final_x_values = list(gene_diff_df['x_vals'])
    final_y_values = list(gene_diff_df['y_vals'])
    final_names = list(gene_diff_df['gene name'])

    top_x_values= final_x_values[:100]
    top_y_values= final_y_values[:100]
    top_row_labels = final_names[:100]


    bottom_x_values= final_x_values[100:]
    bottom_y_values= final_y_values[100:]
    bottom_row_labels = final_names[100:]

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

def create_pca_scatter(principal_df, selected_indices, x_dd_value, y_dd_value):
    '''
    Creates the scatter plot based on the PCA data
    '''

    updated_scatter = go.Figure()

    selected_indices = np.array(selected_indices)

    x_values = np.array(principal_df['principal component 1'])
    y_values = np.array(principal_df['principal component 2'])

    x_adj = abs(min(min(x_values), 0))
    y_adj = abs(min(min(y_values), 0))

    x_values = x_values + x_adj
    y_values = y_values + y_adj

    top_x_values = x_values[selected_indices]
    top_y_values = y_values[selected_indices]
    top_row_labels = row_labels_arr[selected_indices.astype(int)]

    unselected_indices = np.delete(range(0,len(x_values)), selected_indices)
    bottom_x_values = x_values[unselected_indices]
    bottom_y_values = y_values[unselected_indices]
    bottom_row_labels = row_labels_arr[unselected_indices]


    updated_scatter.add_trace(go.Scatter(
        x = bottom_x_values,
        y = bottom_y_values,
        text = bottom_row_labels,
        mode='markers',
        marker_color='blue',
        marker={
            'size':2
        },


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


    ))

    updated_scatter.update_layout(
        title={
            'text': "Plot Title",
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        autosize=True,
        xaxis=go.layout.XAxis(
            title=str(', '.join(x_dd_value) + ' First Principal Component'),
            showgrid=True,
            ),
        yaxis=go.layout.YAxis(
            title=str(', '.join(y_dd_value) + ' First Principal Component'),
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
    updated_violin = go.Figure()

    reduced_x_df = main_df[x_dd_value]
    reduced_y_df = main_df[y_dd_value]
    x_names = " ".join(x_dd_value)
    y_names = " ".join(y_dd_value)



    show_legend = [True] + [False] * (violin_gene_number - 1)

    for index, value  in enumerate(current_row_labels):
        index_in_df = row_dict_hold[value]
        updated_violin.add_trace(go.Violin(
                    x=[value for i in range(len(reduced_x_df.columns))],
                    y=list(reduced_x_df.iloc[index_in_df]),
                    side='negative',
                    legendgroup=x_names,
                    name=x_names,
                    line_color = 'green',
                    showlegend = show_legend[index]
                )
            )
        updated_violin.add_trace(go.Violin(
                    x=[value for i in range(len(reduced_y_df.columns))],
                    y=list(reduced_y_df.iloc[index_in_df]),
                    side='positive',
                    legendgroup = y_names,
                    name =y_names,
                    line_color = 'purple',
                    showlegend = show_legend[index]


                )
            )

    updated_violin.update_layout(
            title={
                'text': "Violin Plot",
                'y':0.9,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            showlegend=True,
            autosize=True,
            violingap=0, 
            violinmode='overlay',
            yaxis_type="log",

        )
    updated_violin.update_traces(points = 'all')
    return updated_violin


'''
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
    print(x_dd_value)
    print(y_dd_value)
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

    if(individ_or_gene_dd == 'Gene'):
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
                    print('Used Cache')
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



    print('Finished Updating')
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
    print(main_scatter_hover)
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
        print(current_key)
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


if __name__ == '__main__':

    app.run_server()