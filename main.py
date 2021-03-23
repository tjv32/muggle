import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
import plotly.graph_objects as go

from make_app import app

from tabs.tab_1 import tab_1
from tabs.tab_2 import tab_2
from tabs.tab_3 import tab_3


server = app.server

app.layout = html.Div(
	children = [
		dcc.Store(id='temp_upload'),
		dcc.Store(id='current_data'),
	    dcc.Store(id='current_selection_view'),
	    dcc.Store(id='cluster_data'),
		dcc.Tabs(
			[
				tab_1,
				tab_2,
				tab_3
			]
		)
	]
)

app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

if __name__ == '__main__':

    app.run_server(debug=True)
