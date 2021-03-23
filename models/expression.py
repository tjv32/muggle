import numpy as np
import pandas as pd

dataframe_keys = ['df_with_gene_names', 'df']
list_dict_keys = ['gene_names', 'gene_names_dd',
					'single_columns', 'single_columns_dd',
					'group_columns_dd', 'all_columns_dd',
					'mean_vals', 'std_dev_vals',
					'group_dict']

class expression:
	def __init__(self, df):
		self.df_with_gene_names = df.reset_index().drop("index", axis=1)
		
	def organize_from_initial_upload(self, group_dict):

		self.group_dict = group_dict
		df = self.df_with_gene_names.drop('Gene_Name', 1)
		self.df = df

		gene_names = self.df_with_gene_names['Gene_Name']
		self.gene_names = gene_names
		self.gene_names_dd = [{'label': i, 'value': i} for i in gene_names]

		single_columns = list(df.columns)
		self.single_columns = single_columns

		single_columns_dd = [{'label': i, 'value': i} for i in single_columns]
		self.single_columns_dd = single_columns_dd

		group_columns_dd = [{'label': i, 'value': ' '.join(group_dict[i])} for i in group_dict.keys()] 
		self.group_columns_dd = group_columns_dd

		self.all_columns_dd = group_columns_dd + single_columns_dd

		self.mean_vals = df.mean(axis=1)
		self.std_dev_vals = df.std(axis=1)

	def organize_from_store_upload(self, info_dict):
		for key in dataframe_keys:
			setattr(self, key, pd.DataFrame(info_dict[key]))

		for key in list_dict_keys:
			setattr(self, key, info_dict[key])

	def create_store_obj(self):
		store_obj = {}

		for key in dataframe_keys:
			store_obj[key] = getattr(self, key).to_dict()

		for key in list_dict_keys:
			store_obj[key] = getattr(self, key)

		return store_obj

	def update_groups(self, new_groups):
		self.group_dict = new_groups
		group_columns_dd = [{'label': i, 'value': ' '.join(new_groups[i])} for i in new_groups.keys()] 
		self.group_columns_dd = group_columns_dd

		self.all_columns_dd = group_columns_dd + self.single_columns_dd
