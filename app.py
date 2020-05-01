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

cache = load_t_test_cache()

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
        children='CWRU',
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
        value='CTRL_0',
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
        value='CTRL_1',
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
body = dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dcc.Store(id='heatmap_memory'),
                            dcc.Store(id='scatter_plot_memory'),
                            dcc.Store(id='name_memory'),
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

app.layout = html.Div([body])


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
            Input('update_button', 'n_clicks')
        ],
        state=[
            State('x_axis_selection_dd', 'value'),
            State('y_axis_selection_dd', 'value'),
            State('gene_select_dd', 'value'),
            State('invidual_or_gene_dd', 'value'),
            State('violin_gene_number', 'value'),
            State('combination_selection_dd', 'value'),
            State('name_memory', 'data'),
            State('scatter_plot_memory', 'data'),

        ]
)
def update_button(n_clicks, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd, violin_gene_number, combination_dd_value, name_memory, scatter_mem):
    '''
        This callback takes all the dropdown menus and will calculate the differentially expressed genes and update the memory objects
    '''
    x_dd_value = combine_columns(x_dd_value)
    y_dd_value = combine_columns(y_dd_value)
    x_num = len(x_dd_value)
    y_num = len(y_dd_value)



    if(x_num == 0 or y_num == 0):

        raise PreventUpdate

    x_dff = df.loc[:, x_dd_value]
    y_dff = df.loc[:, y_dd_value]

    if(individ_or_gene_dd == 'Individual'):

        if(name_memory is not None and x_dd_value == name_memory['x_dd'] and y_dd_value == name_memory['y_dd']
            and combination_dd_value == name_memory['combo'] and individ_or_gene_dd == name_memory['type']
            and violin_gene_number == name_memory['violin_gene_number']):

            gene_diff_df = pd.DataFrame(scatter_mem['df1'])

        else:

            gene_diff_df = find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value)

    else:

        gene_diff_df = None
    
    heat_dict = create_heatmap_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd)

    scatter_dict = create_scatter_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd)

    x_display, y_display = create_x_y_group_display_names(x_dd_value, y_dd_value)

    name_memory = {
        'x_dd' : x_dd_value,
        'y_dd' : y_dd_value,
        'type' : individ_or_gene_dd,
        'combo' : combination_dd_value,
        'violin_gene_number':violin_gene_number
        }

    return heat_dict, scatter_dict, x_display, y_display, name_memory, None


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

    heatmap_Data = go.Heatmap(
            z=name_adj_normal_heat_df,
            x=name_adj_normal_heat_df.columns,
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
        ]
)

def update_violin(scatter_selected_data, memory, violin_gene_number, name_memory):
    '''
    This callback uses the heat map memory store and any selected data in the scatter plot and heat map to update the violin plot
    '''
    print(scatter_selected_data)
    if(memory is None):
        raise PreventUpdate

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

    violin_figure = update_violin_plot(current_row_labels, name_memory['x_dd'], name_memory['y_dd'], violin_gene_number)

    

    return violin_figure


#Helper Functions

def combine_columns(column_selection):
    '''
    Ensures that input x and y groups are lists of strings. Grabs the default column groups underlying columns and ensures no repeats within a selection.
    Sorts the columns alphabetically
    '''
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

def find_differentially_expressed_genes(x_dff, y_dff, combination_dd_value):
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
        difference_measure = np.array(df.apply(row_t_test, args = [x_dff.columns, y_dff.columns], axis = 1))
        absolute_difference_measure = np.absolute(difference_measure)

    gene_diff_df = pd.DataFrame({
        'original_index': range(0, len(difference_measure)),
        'difference measure' : difference_measure,
        'absolute difference measure' : absolute_difference_measure,
        'x_vals':x_combo,
        'y_vals' : y_combo,
        'gene name' : row_labels,
        'means' : df_row_means,
        'stds' : df_row_stds,
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

def create_heatmap_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd):
    '''
    Uses the gene_diff_df to create a heat df and normalized df that will be used to create heatmap
    '''
    if(gene_diff_df is None):
        if(type(gene_select_dd) == str):
            gene_select_dd = [gene_select_dd]
        gene_select_indices = []
        for val in gene_select_dd:
            gene_select_indices.append(row_dict[val])
        heat_df = df.iloc[df.index.isin(gene_select_indices)]

    else:
        gene_select_indices = list(gene_diff_df['original_index'])[:100]
        heat_df = df.iloc[list(gene_diff_df['original_index'])[:100]]
    
    current_row_labels = [row_labels[index] for index in gene_select_indices]

    name_adj_all_columns = [s + '_X' for s in x_dd_value] + [s + '_Y' for s in y_dd_value]

    heat_means = df_row_means[gene_select_indices].values
    heat_stds = df_row_stds[gene_select_indices].values

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

def create_scatter_memory(gene_diff_df, x_dd_value, y_dd_value, gene_select_dd, individ_or_gene_dd):
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
        title='Scatter Plot',
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
        title='Scatter Plot',
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

def update_violin_plot(current_row_labels, x_dd_value, y_dd_value, violin_gene_number):
    '''
    Graphs the selected indices in a violin plot
    '''
    updated_violin = go.Figure()

    reduced_x_df = df[x_dd_value]
    reduced_y_df = df[y_dd_value]
    x_names = " ".join(x_dd_value)
    y_names = " ".join(y_dd_value)



    show_legend = [True] + [False] * (violin_gene_number - 1)

    for index, value  in enumerate(current_row_labels):
        index_in_df = row_dict[value]
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
    app.run_server(debug=True)