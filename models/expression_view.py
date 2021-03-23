import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from helper_functions import *

dataframe_keys = ['combined_df']
list_dict_keys = ['group_one_dd', 'group_two_dd', 'group_one_dd_suffix', 'group_two_dd_suffix', 
				  'variance_filter', 'fold_filter', 'combo_type', 'bottom_right_visualization', 
				  'gene_dd', 'default_groups']

class expression_view:

	def __init__(self):
		pass

	def init_from_initial(self, input_dict):
		self.group_one_dd = input_dict['group_one_dd']
		self.group_two_dd = input_dict['group_two_dd']
		self.bottom_right_visualization = input_dict['bottom_right_visualization']
		self.comparison_selection_dd = input_dict['comparison_selection_dd']
		self.default_groups = input_dict['exp'].group_dict

		group_one_df = input_dict['exp'].df[input_dict['group_one_dd']].add_suffix("_1")
		group_two_df = input_dict['exp'].df[input_dict['group_two_dd']].add_suffix("_2")

		combined_df = pd.concat([group_one_df, group_two_df], axis=1)
		combined_df['Gene'] = list(input_dict['exp'].gene_names)
		combined_df['mean_vals'] = list(input_dict['exp'].mean_vals)
		combined_df['std_dev_vals'] = list(input_dict['exp'].std_dev_vals)


		combined_df['t_vals'] = [0] * len(combined_df)

		if(input_dict['comparison_selection_dd'] == 'Condition Focused'):
			combined_df['t_vals'] = list(combined_df.apply(row_t_test, args = [group_one_df.columns, group_two_df.columns], axis = 1))
			combined_df = combined_df.sort_values('t_vals', ascending = False).reset_index(drop=True)
			self.gene_dd = None
		else:
			self.gene_dd = input_dict['gene_dd']


		self.combined_df = combined_df

		self.group_one_dd_suffix = group_one_df.columns
		self.group_two_dd_suffix = group_two_df.columns

		self.variance_filter = input_dict['variance_filter']
		self.fold_filter = input_dict['fold_filter']
		self.combo_type = input_dict['combo_type']

	def init_from_store(self, store_obj):
		for key in dataframe_keys:
			setattr(self, key, pd.DataFrame(store_obj[key]))

		for key in list_dict_keys:
			setattr(self, key, store_obj[key])

	def create_store_obj(self):
		store_obj = {}

		for key in dataframe_keys:
			store_obj[key] = getattr(self, key).to_dict()

		for key in list_dict_keys:
			store_obj[key] = getattr(self, key)

		return store_obj

	def create_heatmap_graph(self):
		if(self.gene_dd is not None):
			heatmap_df = self.combined_df[self.combined_df.Gene.isin(self.gene_dd)]
		else:
			heatmap_df = self.combined_df[:50]

		x_vals = list(self.group_one_dd_suffix  + self.group_two_dd_suffix)
		y_vals = list(heatmap_df['Gene'])
		z_vals = heatmap_df[self.group_one_dd_suffix  + self.group_two_dd_suffix].sub(heatmap_df['mean_vals'], axis=0).div(heatmap_df['std_dev_vals'], axis=0).values

		heatmap_data = go.Heatmap(
			x=x_vals,
			y=y_vals, 
			z=z_vals,
			colorscale='RdBu_r',
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
			data=heatmap_data, 
			layout=heatmap_layout
		)

		return heatmap_figure

	def create_scatter_graph(self):

		scatter_df = self.combined_df[list(self.group_one_dd_suffix) + list(self.group_two_dd_suffix)]
		variance_filter = scatter_df.var(axis = 1) > self.variance_filter
		scatter_df['Gene'] = list(self.combined_df['Gene'])
		scatter_df = scatter_df[variance_filter]

		group_one_scatter = scatter_df[list(self.group_one_dd_suffix)]
		group_two_scatter = scatter_df[list(self.group_two_dd_suffix)]

		if(self.combo_type == 'median'):
			one_metric = np.array(group_one_scatter.median(axis=1))
			two_metric = np.array(group_two_scatter.median(axis=1))
		else:
			one_metric = group_one_scatter.mean(axis=1)
			two_metric = group_two_scatter.mean(axis=1)

		fold_change = np.array(np.log2(two_metric)) - np.array(np.log2(one_metric))
		final_fold_filter = np.greater(fold_change, self.fold_filter) | np.less(fold_change, -self.fold_filter)
		
		scatter_df = pd.concat([
			scatter_df[['Gene']], 
			pd.DataFrame(
				{
					'one_metric':one_metric, 
					'two_metric':two_metric
				}
			)
		], axis=1)

		scatter_df = scatter_df[final_fold_filter]

		if(self.gene_dd is not None):
			gene_list = self.gene_dd
		else:
			gene_list = list(self.combined_df[:50].Gene)

		mask = scatter_df.Gene.isin(gene_list)
		
		diff_df = scatter_df[mask]
		non_diff_df = scatter_df[~mask]

		non_diff_df = non_diff_df.sample(frac = 0.35)

		scatter_fig = go.Figure()

		scatter_fig.add_trace(
			go.Scatter(
				x = non_diff_df.one_metric,
				y = non_diff_df.two_metric,
				text = non_diff_df.Gene,
				mode='markers',
				marker_color='blue',
				marker={
					'size':3
				},
				legendgroup='Non-Differentially Expressed Genes',
				name='Non-Differential Genes'
			)
		)

		scatter_fig.add_trace(
			go.Scatter(
				x = diff_df.one_metric,
				y = diff_df.two_metric,
				text = diff_df.Gene,
				mode='markers',
				marker_color='red',
				marker={
					'size':7
				},
				legendgroup='Differentially Expressed Genes',
				name='Differential Genes'
			)
		)
		scatter_fig.update_layout(
			title={
				'x':0.5,
				'xanchor': 'center',
				'yanchor': 'top'
			},
			autosize=True,
			xaxis=go.layout.XAxis(
				title="Group One",
				showgrid=True,
			),
			yaxis=go.layout.YAxis(
				title="Group Two",
				showgrid=True,  
			),
			xaxis_type="log", 
			yaxis_type="log",
			showlegend=True
		)
		return scatter_fig


	def create_rug_graph(self):
		color_options = ['green', 'purple', 'orange', 'black', 'yellow', 'white', 'grey']

		rug_gene_count = 5

		if(self.gene_dd is not None):
			rug_df = self.combined_df[self.combined_df.Gene.isin(self.gene_dd)]
		else:
			rug_df = self.combined_df[:rug_gene_count]
		rug_df = rug_df.reset_index()

		temp_df = rug_df[list(self.group_one_dd_suffix + self.group_two_dd_suffix)].copy()
		temp_df.rename(columns = lambda x : str(x)[:-2], inplace=True)
		temp_df['Gene'] = list(rug_df['Gene'])
		rug_df = temp_df

		group_list = []
		for key in self.default_groups:
			group_found = True
			for val in self.default_groups[key]:
				if(val not in list(rug_df.columns)):
					group_found = False
			if(group_found):
				group_list.append(key)


		show_legend = [True] + [False] * (rug_gene_count - 1)
		rug_figure = go.Figure()
		for i, row in rug_df.iterrows():
			current_gene = row['Gene']
			for j, group_key in enumerate(group_list):
				col_names = self.default_groups[group_key]

				rug_figure.add_trace(go.Scattergl(
					x=row[col_names],
					y=[current_gene] * len(col_names),
					name=group_key,
					showlegend = show_legend[i],
					marker=dict(
						color=color_options[j],
						line_color=color_options[j],
					),
					hovertemplate =
						'<b>Expression</b>: %{x}'+
						'<br><b>Gene</b>: %{y}<br>'+
						'<b>Condition</b>  %{text}',
					text = col_names
					 )
				)


		rug_figure.update_traces(mode='markers', marker=dict(line_width=3, symbol=42, size=16))
		rug_figure.update_layout(xaxis_type="log")

		return rug_figure

	def create_histogram_graph(self):
		hist_count = 4
		
		if(self.gene_dd is not None):
			hist_df = self.combined_df[self.combined_df.Gene.isin(self.gene_dd)]
		else:
			hist_df = self.combined_df[:hist_count]
		hist_df = hist_df.reset_index()
		
		hist_figure = make_subplots(rows=1, cols=len(hist_df), subplot_titles=list(hist_df['Gene']))
		show_legend = [True] + [False] * len(hist_df)
		for i, row in hist_df.iterrows():
			hist_figure.add_trace(
				go.Histogram(
					x = list(row[self.group_one_dd_suffix]),
					showlegend = show_legend[i],
					nbinsx=10,
					marker_color = "Green",
					name = "Group One",
				),
				row=1,
				col = i + 1,
			)
			hist_figure.add_trace(
				go.Histogram(
					x = list(row[self.group_two_dd_suffix]),
					showlegend = show_legend[i],
					nbinsx=10,
					marker_color = "Purple",
					name = "Group Two",
				),
				row=1,
				col = i + 1,
			)
		hist_figure.update_layout(barmode='overlay')
		hist_figure.update_traces(opacity=0.65)
		return hist_figure
	def generate_csv(self):
		g_df = self.combined_df[list(self.group_one_dd_suffix + self.group_two_dd_suffix)].copy()
		g_df['Gene'] = list(self.combined_df['Gene'])

		if(self.gene_dd is not None):
			gene_list = self.gene_dd
		else:
			gene_list = list(self.combined_df[:50].Gene)

		mask = g_df.Gene.isin(gene_list)
		g_df = g_df[mask]

		g_df = g_df.set_index('Gene')

		return g_df