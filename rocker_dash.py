import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
from modules.rocker_4_refiner import plot_data
import os

class rocker_gui:
	def __init__(self):
		external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
		self.app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
		self.server = self.app.server
		
		self.plot_data = plot_data()
		self.selected_rl = None
		
		
	def set_layout(self):
		self.app.layout = html.Div([
			dcc.Tabs(id="tabs", value='loader', children=[
				dcc.Tab(label='ROCkOut Management', value='loader'),
				dcc.Tab(label='Protein Manager', value='prot_sel'),
				dcc.Tab(label='2D Scatter Plots', value='twoD'),
				dcc.Tab(label='3D Plot', value='threeD'),
				#dcc.Tab(label='Multiple Alignment', value='MA'),
			]),
			html.Div(id='control_tabs'),
			
			#To avoid throwing errors in app callbacks, we need to pre-declare all of the ids we'll be using.
			#However, we don't want these to ever display - the style={display:none} flags these as hidden, 
			#effectively permanently, since we assign new values to these same IDs and don't update these initial declarations.
			html.Div(id='loader_init', children=[
			dcc.Upload(id = 'project_selector', filename=None),
			dcc.Store(id = 'project_index', data={"index_file":None}),
			
			dcc.Checklist(id='active_proteins', options = ['Load a project!'], value = ['Load a project!']),
			dcc.Store(id = 'active_store', data={"active_prots":None}),
			
			
			dcc.Dropdown([""], value = "", id='existing_rocker_models'),
			html.Button('Get an existing model', id='download_model', n_clicks=0),
			], style={'display': 'none'}),
			
			html.Div(id='2d_init', children=[
			dcc.Graph(id="2d_plot",figure={}),
			dcc.Dropdown(id="2d_rl_dd")
			], style={'display': 'none'})
		])
			
			
		@self.app.callback(Output('control_tabs', 'children'),
					  Input('tabs', 'value'))
					  
		def render_content(tab):
			if tab == 'loader':
				left_half = ""
				return html.Div([
				html.H3("Welcome to ROCkOut!"),
				html.P('''
				This module of ROCkOut is designed to allow you to inspect the results of your model building efforts so far, to refine the 
				protein set used by your ROCkOut model, and to produce a final filter based on your selections. You can also view the data
				that was used to create other ROCkOut models or to edit their final protein sets to your needs.
				'''),
				html.Br(),
				html.P('''ROCkOut projects contain an index file in their top-level. To load a ROCkOut project in this GUI, select or drag-and-drop 
				the index file using the buttons below. When a project is loaded, the other tabs will update to show data from the project.
				'''),
				html.Br(),
				html.P('''You can also download an existing ROCkOut project through the menu below. Since existing projects already contain a filter,
				you do not need to use this GUI to use an existing model as a filter for your reads.
				'''),
				html.Br(),
				
				dcc.Upload(id = 'project_selector', filename=None,
				children=html.Div(["Drag and Drop or ", html.A("Select a file")]),
				 style={
                "width": "90%",
                "height": "60px",
                "lineHeight": "60px",
                "borderWidth": "1px",
                "borderStyle": "dashed",
                "borderRadius": "5px",
                "textAlign": "center",
                "margin": "10px",
            }),
			dcc.Store(id = 'project_index', data={"index_file":None}),
			dcc.Dropdown(["some", "options"], value = "some", id='existing_rocker_models'),
			html.Button('Get an existing model', id='download_model', n_clicks=0),
			
			html.Button('Save the current model', id='save_model', n_clicks = 0),
			
				])
			
			elif tab == "prot_sel":
				if self.plot_data.readlens is not None:
					return html.Div([
						html.P('''Using the interactive plots, you can see the effect of removing particular proteins from the dataset.
						This page controls the final selection of proteins in your model. When you want to remove a protein to improve model performance,
						that must happen on this page, not on the plots.
						'''),
						html.Br(),
						dcc.Store(id = 'active_store', data={"active_prots":None}),
						dcc.Checklist(id='active_proteins', options = list(self.plot_data.active_prot_dict.keys()), value = self.plot_data.active_proteins,inline=True)
					])
					
				else:
					return html.Div([
						dcc.Store(id = 'active_store', data={"active_prots":None}),
						dcc.Checklist(id='active_proteins', options = ['Load a project to see active proteins.'], value = ['Load a project to see active proteins.']),
					])
				
				
			elif tab == 'twoD':
				if self.plot_data.readlens is not None:
					if self.selected_rl is None:
						if self.plot_data is not None:
							self.selected_rl = min(self.plot_data.readlens)
						#Both of these elses are fail conditions for unloaded data.
						else:
							return html.Div([
								html.H3('You need to load a ROCkOut model on the ROCkOut Management tab first!')
							])
					else:
						dropdown_opts = ["Show Read Length "+ str(rl) for rl in self.plot_data.readlens]
						
						return html.Div([
							dcc.Graph(id='2d_plot', style={'width': '100vw', 'height': '80vh'}),
							dcc.Dropdown(dropdown_opts, value=self.selected_rl, id="2d_rl_dd")
						])
						
					dropdown_opts = ["Show Read Length "+ str(rl) for rl in self.plot_data.readlens]
					
					return html.Div([
						dcc.Graph(id='2d_plot', style={'width': '100vw', 'height': '80vh'}),
						dcc.Dropdown(dropdown_opts, value=self.selected_rl, id="2d_rl_dd")
					])
					
				else:
					return html.Div([
						html.H3('You need to load a ROCkOut model on the ROCkOut Management tab first!')
					])
					
			elif tab == 'threeD':
				#Probably need an oops case.
				if "scatter_3d" in self.plot_data.figs:
					return html.Div([
						dcc.Graph(
							id='3d_graph',
							figure=self.plot_data.figs["scatter_3d"],
							#Set size to 100% of screen.
							style={'width': '93vw', 'height': '93vh'}
						)
					])
				else:
					return html.Div([
						html.H3('You need to load a ROCkOut model on the ROCkOut Management tab first!')
					])
					
					
			elif tab == 'MA':
				if "MA" in self.plot_data.figs:
					return html.Div([
						html.H3('Tab content 4'),
						dcc.Graph(
							id='graph-4-tabs-dcc',
							figure={
								'data': [{
									'x': [1, 2, 3],
									'y': [5, 10, 6],
									'type': 'bar'
								}]
							}
						)
					])
				else:
					return html.Div([
						html.H3('You need to load a ROCkOut model on the ROCkOut Management tab first!')
					])
		
		@self.app.callback(Output('active_store', 'data'),
			Input('active_proteins', 'value'))
		def update_actives(value):
			if value[0] == 'Load a project!' or value is None:
				pass
			else:
				self.plot_data.active_proteins = value
			return {"active_prots":value}
		
		@self.app.callback(Output('2d_plot', 'figure'),
			Input('2d_rl_dd', 'value'))
		def update_2d(read_len):
			try:
				#Remove leading dropdown text
				stripped = read_len.split("Show Read Length ")[1]
				self.selected_rl = stripped
				cur_plot = "scatter_2d_read_len_" + self.selected_rl
				return self.plot_data.figs[cur_plot]
			except:
				return {}
				
		@self.app.callback(Output('existing_rocker_models', 'value'),
			Input('download_model', 'n_clicks'),
			State('existing_rocker_models', 'value')
		)
		def get_old_model(n_clicks, value):
			#We can download something here.
			return value
			
		@self.app.callback(Output('project_index', 'data'),
			Input("project_selector", 'filename'))
		def load_project(filename):
			if filename is not None:
				wor_dir = os.listdir(".")
				proj_dir = None
				for d in wor_dir:
					if os.path.isdir(d):
						if os.path.exists(os.path.normpath(d+"/"+filename)):
							proj_dir = d
							break
				
				try:
					#rocker.plot_data.set_project(proj_dir)
					self.plot_data.run_plotter(proj_dir)
				except:
					pass
					
				return {"index_file":filename}
			else:
				pass
		
		
def main():
	rocker = rocker_gui()
	rocker.set_layout()
	#rocker.plot_data.run_plotter("rockout_out/")
	rocker.app.run_server(debug=True,port=8056)
		
if __name__ == '__main__':
    main()