import sys
import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go

class rockout_visualizer:
	def __init__(self):
		pass
		
	#Proportionally sample reads to show positive/negative patterns without bloated file size.
	def sample_to_reasonable(self, df, reasonable_size = 12000):
		max_pos = int(reasonable_size/3)
		positives = df.loc[df['classifier'] == "Positive"]
		if len(positives.index) > max_pos:
			positives = positives.sample(n = max_pos, replace = False)
		positives = positives.reset_index(drop = True)
		
		df = df.loc[df['classifier'] != "Positive"]
		
		df = df.reset_index(drop=True)
		a_reasonable_number = reasonable_size - len(positives.index)
		
		total_rows = len(df.index)
		frac = (a_reasonable_number/total_rows)
		
		#Proportionally sample categories down to 12k reads minus # positives
		df = df.groupby('classifier').apply(pd.DataFrame.sample, frac=frac, replace=False).reset_index(level='classifier', drop=True)
		
		df = pd.concat([positives, df])
		
		df = df.reset_index(drop=True)
		
		return df
		
	def visualize_reads_and_roc_curves(self, df, x, y, x_curve, y_curve, xname = "Position in Multiple Alignment", yname = "Bitscore"):
		x = np.array(x)
		y = np.array(y)

		reads_plot = px.scatter(df, x=df[x], y=df[y], 
						color = df['classifier'], hover_data = [df['read_id'], df['target']],
						color_discrete_sequence = ["blue", "aqua", "darkorange", "green"], 
						category_orders={"classifier": ["Positive", "Homology_Target", "Non_Target", "Negative"]})
		cutoff_line = px.line(x = x_curve, y = y_curve)
		
		final_figure = go.Figure(data=reads_plot.data + cutoff_line.data)
		final_figure.update_xaxes(title=dict(text=xname))
		final_figure.update_yaxes(title=dict(text=yname))
		
		return final_figure