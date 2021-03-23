import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_html_components as html
from dash_extensions import Download
import plotly.graph_objects as go

import pandas as pd
import numpy as np
import io

from models.expression import expression
from models.expression_view import expression_view
from helper_functions import *
from tabs.default_data import *

from make_app import app


   

tab_title = html.H1(
    children='Multi Group Gene Explorer (MUGGLE)',
    style={
        'textAlign': 'center',
        'color': 'black'
        }
    )

heatmap_layout = go.Layout(
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
heatmap_figure = go.Figure(
	data=None, 
	layout=heatmap_layout
)
heatmap_graph = dcc.Graph(
    id='heatmap',
    figure=heatmap_figure,
    style={
        'width':'50%',
        'display': 'inline-block'
    }
)

scatter_layout = go.Layout(
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
scatter_figure = go.Figure(
	data=None, 
	layout=scatter_layout
)
scatter_graph = dcc.Graph(
    id='scatter',
    style={
        'width':'50%',
        'display': 'inline-block'
    }
)

rug_layout = go.Layout(
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
rug_figure = go.Figure(
	data=None,
	layout=rug_layout
)
rug_graph = dcc.Graph(
    id='rug',
    style={
        'width':'100%',
        'display': 'inline-block'
    }
)

comparison_selection_label = html.Label(
    children = 'Gene Focused or Condition Focused Comparison',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
comparison_selection_dd = dcc.Dropdown(
    id='comparison_selection_dd',
   	options=[{'label': i, 'value': i} for i in ['Gene Focused', 'Condition Focused']],
    value='Gene Focused',
    multi=False,
)

group_one_label = html.Label(
    children = 'Select Group One',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
group_one_dd = dcc.Dropdown(
    id='group_one_selection_dd',
    options=default_exp.all_columns_dd,
    multi=True
)
group_one_display = html.Label(
    children = [html.Br(), 'Current Group One'],
    id = 'group_one_display',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)

group_two_label = html.Label(
    children = 'Select Group Two',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
group_two_dd = dcc.Dropdown(
    id='group_two_selection_dd',
    options=default_exp.all_columns_dd,
    multi=True
)
group_two_display = html.Label(
    children = [html.Br(), 'Current Group Two'],
    id = 'group_two_display',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)

gene_selection_label = html.Label(
    children = 'Gene Selection',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
gene_selection_dd = dcc.Dropdown(
    id='gene_select_dd',
    options=default_exp.gene_names_dd,
    multi=True
)

combo_type_label = html.Label(
    children = 'Select combination type',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
combo_type_dd = dcc.Dropdown(
    id='combination_selection_dd',
    options=[
        {'label': 'mean', 'value': 'mean'},
        {'label': 'median', 'value': 'median'}
    ],
    value= 'mean',
    multi=False,
)

bottom_right_visualization_label = html.Label(
    children = 'Bottom Right Visualization',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
bottom_right_visualization_dd = dcc.Dropdown(
    id='bottom_right_visualization',
    options=[{'label': i, 'value': i} for i in ['Histogram', 'Rug']],
    value='Rug',
    multi=False,
)

fold_label = html.Label(
    children = 'Select Fold Cut Off',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
fold_filter = html.Div(
	[
		dcc.Input(
		    id='fold_filter',
		    type='number',
		    value=0,
		    min=0,
		    max=100
		)
	]
)

variance_label = html.Label(
    children = 'Select Variance Cut Off',
    style={
        'textAlign': 'left',
        'color': 'black',
    }
)
variance_filter = html.Div(
	[
		dcc.Input(
		    id='variance_filter',
		    type='number',
		    value=0,
		    min=0,
		    max=10
		)
	]
)

update_button = html.Button(
    'Update', 
    id='update_button_main',
)

download_heat_button = html.Div(
	[
		html.Button(
			"Download csv", 
			id="btn"
		), 
		Download(id="download")
	]
)

body_1 = html.Div(
	[
		dbc.Container(
 			[
            	dbc.Row(
                	[
                    	dbc.Col(
                        	[
                            	tab_title
                        	]
                    	)
                	]
            	),
            	dbc.Row(
                	[
                    	dbc.Col(
                        	[
                            	heatmap_graph,
                            	scatter_graph
                        	]
                    	)
                	]
            	),
            	dbc.Row(
                	[
                    	dbc.Col(
                        	[
	                            comparison_selection_label,
	                            comparison_selection_dd,
	                            group_one_label,
	                            group_one_dd, 
	                            group_two_label,
	                            group_two_dd,
	                            gene_selection_label,
	                            gene_selection_dd,
                                bottom_right_visualization_label,
                                bottom_right_visualization_dd,
	                            combo_type_label,
	                            combo_type_dd,
	                            fold_label,
	                            fold_filter,
	                            variance_label,
	                            variance_filter,
	                            update_button,
	                            group_one_display,
	                            group_two_display,
	                            download_heat_button,
                        	],
                        	width = 4
                    	),
                    	dbc.Col(
                        	[
                            	rug_graph,
                        	], 
                        	width = 8
                    	),
                	]
            	)
        	],
        	fluid = True
    	)
    ]
)

tab_1 = dcc.Tab(
	label='Comparisons', 
	children=[
        body_1
    ]
)


@app.callback(
    [
        Output('current_selection_view', 'data')
    ],
    [
        Input('update_button_main', 'n_clicks')
    ],
    [
        State('current_data', 'data'),
        State('comparison_selection_dd', 'value'),
        State('group_one_selection_dd', 'value'),
        State('group_two_selection_dd', 'value'),
        State('gene_select_dd', 'value'),
        State('bottom_right_visualization', 'value'),
        State('variance_filter', 'value'),
        State('fold_filter', 'value'),
        State('combination_selection_dd', 'value')
    ]
)
def update_main_view_blah(n_clicks, current_data, comparison_selection_dd, group_one_selection_dd, group_two_selection_dd, gene_select_dd, bottom_right_visualization, variance_filter, fold_filter, combination_selection_dd):

    if(not all([current_data, comparison_selection_dd, group_one_selection_dd, group_two_selection_dd, bottom_right_visualization])):
        raise PreventUpdate

    if(gene_select_dd is None and comparison_selection_dd == 'Gene Focused'):
        raise PreventUpdate

    if(current_data == 'default'):
        exp = default_exp
    else:
        exp = expression(pd.DataFrame(None))
        exp.organize_from_store_upload(current_data)
    
    input_dict = {
        'comparison_selection_dd' : comparison_selection_dd,
        'group_one_dd' :  clean_column_names(group_one_selection_dd),
        'group_two_dd' :  clean_column_names(group_two_selection_dd),
        'gene_dd' :  clean_gene_list(gene_select_dd),
        'bottom_right_visualization' : bottom_right_visualization,
        'variance_filter' : variance_filter,
        'fold_filter' : fold_filter,
        'combo_type' : combination_selection_dd,
        'exp' : exp
    }

    exp_view = expression_view()
    exp_view.init_from_initial(input_dict)

    return [exp_view.create_store_obj()]


@app.callback(
    [
        Output('heatmap', 'figure')
    ],
    [
        Input('current_selection_view', 'data')
    ],
)
def update_heatmap(current_selection_view):
    
    if(current_selection_view is None):
        raise PreventUpdate

    exp_view = expression_view()
    exp_view.init_from_store(current_selection_view)

    return [exp_view.create_heatmap_graph()]

@app.callback(
    [
        Output('scatter', 'figure')
    ],
    [
        Input('current_selection_view', 'data')
    ],
)
def update_scatter(current_selection_view):
    
    if(current_selection_view is None):
        raise PreventUpdate

    exp_view = expression_view()
    exp_view.init_from_store(current_selection_view)

    return [exp_view.create_scatter_graph()]

@app.callback(
    [
        Output('rug', 'figure')
    ],
    [
        Input('current_selection_view', 'data')
    ],
)
def update_scatter(current_selection_view):
    
    if(current_selection_view is None):
        raise PreventUpdate

    exp_view = expression_view()
    exp_view.init_from_store(current_selection_view)

    return_figure = None
    if(exp_view.bottom_right_visualization == "Rug"):
        return_figure = exp_view.create_rug_graph()
    else:
        return_figure = exp_view.create_histogram_graph()
    return [return_figure]

@app.callback(
        Output("download", "data"),
        [
            Input("btn", "n_clicks"),
        ],
        [
            State('current_selection_view', 'data'),
        ]
)
def generate_csv(n_nlicks, current_selection_view):
    if(current_selection_view is None):
        raise PreventUpdate
    # Convert data to a string.
    s = io.StringIO()

    exp_view = expression_view()
    exp_view.init_from_store(current_selection_view)

    g_df = exp_view.generate_csv()
    g_df.to_csv(s)
    content = s.getvalue()
    # The output must follow this form for the download to work.
    return dict(filename="Heatmap_data.csv", content=content, type="text/csv")  