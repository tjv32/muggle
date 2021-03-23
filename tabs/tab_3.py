import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_html_components as html
import plotly.graph_objects as go

import pandas as pd
import io
import base64
import pathlib

from models.expression import expression
from tabs.default_data import *

from make_app import app

tab_title = html.H1(
    children='Upload Data as a CSV',
    style={
        'textAlign': 'center',
        'color': 'black'
    }
)

upload = html.Div(
	[
	    dcc.Upload(
	        id='file_upload',
	        children=html.Div(
	        	[
	            	'Drag and Drop or ',
	            	html.A('Select Files')
	        	]
			),
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
	        multiple=False
    	)
	]
)

number_of_groups_labels = html.Label(
    children = 'Number of Custom Groups',
    style={
        'textAlign': 'center',
        'color': 'black',
    }
)
number_of_groups_dd = dcc.Dropdown(
    id='number_of_groups_dd',
    options=[{'label': i, 'value': i} for i in range(1, 10)],
    value=5,
    multi=False,
)

input_msg = html.Div(id='upload_update_msg')

variable_upload_inputs = html.Div(id='variable_upload_inputs')
variable_upload_dropdowns = html.Div(id='variable_upload_dds')

update_button_msg = html.Label(
    children = 'Select Custom Groups After Selecting and Uploading the Data Above',
    style={
        'textAlign': 'center',
        'color': 'black',
    }
)
update_button_groups = html.Button(
    'Update Custom Groups', 
    id='update_button_groups',
)
updated_groups_msg = html.Div(id='updated_groups_msg')

body_3 = dbc.Container(
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
	                    upload,

	                ]
	            )
	        ]
	    ),
    	dbc.Row(
        	[
            	dbc.Col(
                	[
                        update_button_msg,
                        html.Br(),
                        html.Br(),
                    	number_of_groups_labels,
                    	number_of_groups_dd,
                    	input_msg
                	]
            	)
        	]
    	),
    	dbc.Row(
        	[
            	dbc.Col(
                	[
						variable_upload_inputs
                	]
            	),
            	dbc.Col(
                	[
                    	variable_upload_dropdowns
                	]
            	)
        	]
    	),
    	dbc.Row(
        	[
            	dbc.Col(
                	[      
                        update_button_groups,
                        updated_groups_msg
                	]
            	),
        	]
    	),  
	]
)

tab_3 = dcc.Tab(
	label='Custom Data', 
	children=[
        body_3
    ]
)

@app.callback(
    [
        Output('variable_upload_dds', 'children'),
        Output('variable_upload_inputs', 'children')
    ],
    [
        Input('number_of_groups_dd', 'value'),
        Input('temp_upload', 'data')
    ]
)
def update_group_selection(number_groups, current_data):
    if(current_data is None):
        raise PreventUpdate
    exp = expression(pd.DataFrame(None))
    exp.organize_from_store_upload(current_data)

    dropdowns = []
    labels = []
    for num in range(number_groups):
        dropdowns.append(
            dcc.Dropdown(
                id=f'default_group_{num}',
                options=exp.single_columns_dd,
                value=5,
                multi=True,
                style={
                    'height': 50,
                    'width':'100%',
                    'display': 'inline-block'
                }
            )
        ) 
        labels.append(
            dcc.Input(
                id=f"input_{num}",
                type='text',
                placeholder="Input Group Name",
                style={
                    'height': 54,
                    'width':'100%',
                    'display': 'inline-block'
                }
            )
        )

    return dropdowns, labels

@app.callback(  
            [
                Output('temp_upload', 'data'),
                Output('upload_update_msg', 'children')
            ],
            [
                Input('file_upload', 'contents')
            ],
            [
                State('file_upload', 'filename'),
                State('file_upload', 'last_modified')
            ])
def create_upload_display(contents, filename, last_modified):
    if(contents is None):
        raise PreventUpdate
    try:
        content_type, content_string = contents[0].split(',')
        decoded = base64.b64decode(content_string)
            
        input_df = pd.read_csv(
            io.StringIO(decoded.decode('utf-8')))
            
        exp = expression(input_df)
        exp.organize_from_initial_upload(default_groups)
        mem_exp = exp.create_store_obj()

        update_msg = html.Div([
            f'{filename[0]} successfully uploaded'
        ])
    except:
        mem_exp = None
        update_msg = html.Div([
            'There was an error processing this file.'
        ])
    return mem_exp, update_msg


@app.callback(
    [
        Output('gene_select_dd', 'options'),
        Output('cluster_selection_dd', 'options'),
        Output('group_one_selection_dd', 'options'),
        Output('group_two_selection_dd', 'options'),
        Output('input_combo_1', 'options'),
        Output('input_combo_2', 'options'),
        Output('input_combo_3', 'options'),
        Output('input_combo_4', 'options'),
        Output('input_combo_5', 'options'),
        Output('input_combo_6', 'options'),
        Output('input_combo_7', 'options'),
        Output('input_combo_8', 'options'),
    ],
    [
        Input('current_data', 'data'),
    ]
)
def update_dropdowns(current_data):
    if(current_data is None or current_data == 'default'):
        raise PreventUpdate

    exp = expression(pd.DataFrame(None))
    exp.organize_from_store_upload(current_data)

    dd_list = [exp.gene_names_dd] * 2 + [exp.all_columns_dd] * 10

    return dd_list

@app.callback(
    [
        Output('current_data', 'data'),
        Output('updated_groups_msg', 'children'),
    ],
    [
        Input('update_button_groups', 'n_clicks')
    ],
    [
        State('temp_upload', 'data'),
        State('variable_upload_inputs', 'children'),
        State('variable_upload_dds', 'children'),
    ]
)
def update_groups(nclicks, current_data, input_names, input_dd):
    if(current_data is None):
        mem_exp = 'default'
        group_update_msg = html.Div(['Using Default Data'])
        return [mem_exp, group_update_msg]
        
    new_groups = {}
    for name, dd in zip(input_names, input_dd):
        new_dd = dd['props'].get('value')
        new_name = name['props'].get('value')
        if(new_dd is not None and new_name is not None):
            new_members = list(new_dd)
            new_groups[new_name] = new_members

    exp = expression(pd.DataFrame(None))
    exp.organize_from_store_upload(current_data)
    exp.update_groups(new_groups)
    mem_exp = exp.create_store_obj()

    update_items = ["Current Custom Groups: ", html.Br()]
    for key in new_groups.keys():
        update_items.append(
            key + " : " + ", ".join(new_groups[key])
        )

    group_update_msg = html.Div(update_items)


    return [mem_exp, group_update_msg]


