# -*- coding: utf-8 -*-

import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
import numpy as np
import six.moves.urllib.request as urlreq
from six import PY3
import dash_daq as daq
import dash
from dash.exceptions import PreventUpdate
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
# from dash_bio_utils import xyz_reader
import dash_bootstrap_components as dbc
import base64
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.express as px
import sys
import ase
import ase.io as aio
import scipy as sp
from ase import neighborlist
import networkx as nx


def xyz_to_json(mol):
    cutOff = neighborlist.natural_cutoffs(mol)
    cutOff = [cc - 0.2 for cc in cutOff]
    neighborList = neighborlist.NeighborList(
        cutOff, self_interaction=False, bothways=True)
    neighborList.update(mol)

    neighborList.get_neighbors(0)

    bmatrix = neighborList.get_connectivity_matrix(sparse=False)

    dmatrix = sp.spatial.distance_matrix(mol.positions, mol.positions)
    dmatrix[dmatrix > 2] == 0
    g = nx.Graph(dmatrix)
    g2 = nx.Graph(bmatrix)
    list(g.edges(data=True))

    symbols = mol.get_chemical_symbols()

    nodes_data = [{'id': node, 'atom': symbols[node]} for node in g.nodes]

    edge_data = [{'id': num, 'source': at[0], 'target': at[1],
                  'strength': 2 if (at[0], at[1]) in g2.edges() else 0.1,
                  'bond': 1 if (at[0], at[1]) in g2.edges() else 0,
                  'distance':at[2]['weight'] * 20} for num, at in enumerate(list(g.edges(data=True)))]

    data = {'nodes': nodes_data, 'links': edge_data}
    return data


def create_app(data_file='blaskovits2021_example/data.csv',
               structures_dir='blaskovits2021_example/structures'):

    DATA = pd.read_csv(
        './data.csv')


    DATA = DATA.round(5)
    DATA['id'] = list(DATA.index)
    info_dict = [
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**Compound index**: COMP_names'''],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**S1<sub>exc</sub>**: S1_exc &emsp;&emsp;''',
                         '''&emsp; **T1<sub>exc</sub>**: T1_exc eV&emsp;&emsp;''',
                         '''**&Delta;E<sub>ST</sub>**: S1_T1_split'''],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**&Omega;<sup>S<sub>1</sub></sup><sub>A&#8594;A</sub>**: S1AtoA&emsp;''',
                         '''**&Omega;<sup>S<sub>1</sub></sup><sub>A&#8594;D</sub>**: S1AtoD&emsp;''',
                         '''**&Omega;<sup>S<sub>1</sub></sup><sub>D&#8594;A</sub>**: S1DtoA&emsp;''',
                         '''**&Omega;<sup>S<sub>1</sub></sup><sub>D&#8594;D</sub>**: S1DtoD''',
                     ], style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**&Omega;<sup>T<sub>1</sub></sup><sub>A&#8594;A</sub>**: T1AtoA&emsp;''',
                         '''**&Omega;<sup>T<sub>1</sub></sup><sub>A&#8594;D</sub>**: T1AtoD&emsp;''',
                         '''**&Omega;<sup>T<sub>1</sub></sup><sub>D&#8594;A</sub>**: T1DtoA&emsp;''',
                         '''**&Omega;<sup>T<sub>1</sub></sup><sub>D&#8594;D</sub>**: T1DtoD''',
                     ], style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**D<sub>HOMO</sub> **: D_HOMO &emsp;&emsp;''',
                         '''**D<sub>LUMO</sub> **: D_LUMO &emsp;&emsp;''',
                         '''**D<sub>GAP</sub> **: D_GAP''',

                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**A<sub>HOMO</sub> **: A_HOMO &emsp;&emsp;''',
                         '''**A<sub>LUMO</sub> **: A_LUMO &emsp;&emsp;''',
                         '''**A<sub>GAP</sub> **: A_gap ''',

                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**D<sub>S1</sub> **: D_S1 &emsp;&emsp;''',
                         '''**D<sub>T1</sub> **: D_T1 &emsp;&emsp;''',
                         '''**D(&Delta;E<sub>ST</sub>) **: D_split''',

                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**A<sub>S1</sub> **: A_S1 &emsp;&emsp;''',
                         '''**A<sub>T1</sub> **: A_T1 &emsp;&emsp;''',
                         '''**A(&Delta;E<sub>ST</sub>) **: A_split''',

                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''**&phi;<sub>D-A</sub>**: dft_dihedral_norm''',
                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"}),
        html.P([html.Strong("SMILES: "), "SMILES"],
               style={"margin-bottom": "10pt"}),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''*All energy values are reported in eV*''',
                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%", "text-align": "center"},),

    ]

    info_text = [dcc.Markdown('''
        # App instructions
        - Click once to select a point, click twice to deselect it.
        - Change datapoint size with the sliding scale below the plot.
        - **Details**-labels dictionary:
        '''),
                 dbc.CardBody(info_dict,
                              id="information_card",
                              className="mb-3",
                              style={
                                  "max-height": "230pt",
                                  "width": "800px",
                                  "border": "1pt solid #d5d5d5",
                                  "border-radius": "10pt",
                                  "margin": "auto",
                                  "overflow-y": "auto",
                                  "overflow-x": "auto",
                              }
                              ),
                 dcc.Markdown('''
        # Data generation details:

        - 3D structures and molecular properties are obtained from DFT-optimized geometries    (wB97XD/6-31G*)
        - Singlet-triplet splitting &Delta;E<sub>ST</sub> is calculated as follows: &Delta;E<sub>ST</sub> = E(S<sub>1</sub>) - 2E(T<sub>1</sub>)
        - The character of the excited states are evaluated using the charge transfer numbers (&Omega;<sup>E</sup><sub>i&#8594;j</sub>) gathered
        from the transition density matrices of a given excited state E, which express the
        accumulation of hole and electron density on molecular fragments *i* and *j*, respectively.
        Here,   the dimer is partitioned into the donor (D) and acceptor (A) fragments.
        - The dihedral formed between the donor and acceptor cores in the dimer (&phi;<sub>D-A</sub>) is to be between 0° and 90&deg;.
        - Computational data for each dimer can be found in the directory labelled with that compound’s index in the Materials Cloud archive.

        Please see the paper for further details:
        ''', dangerously_allow_html=True, style={"margin-bottom": "-12pt", "fontsize": "10%"}), dcc.Link(
        href='https://chemrxiv.org/articles/preprint/Identifying_the_Trade-off_between_Intramolecular_Singlet_Fission_Requirements_in_Donor-Acceptor_Copolymers/13333475/1',
        target='https://chemrxiv.org/articles/preprint/Identifying_the_Trade-off_between_Intramolecular_Singlet_Fission_Requirements_in_Donor-Acceptor_Copolymers/13333475/1',
        children=['Identifying the Trade-off between Intramolecular Singlet Fission Requirements in Donor-Acceptor Copolymers', ],
        refresh=False,
        style={"margin-bottom": "25pt"},),
        dcc.Markdown(dangerously_allow_html=True,
                     children=[
                         '''The app was built using :''',
                     ],
                     style={"margin-bottom": "-12pt", "fontsize": "10%"},),
        dcc.Link(
        href='https://github.com/lcmd-epfl/molecular_data_explorer',
        target='https://github.com/lcmd-epfl/molecular_data_explorer',
        children=['Molecular Data Explorer', ],
        refresh=False,
    )
    ]

    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    app = dash.Dash(
        __name__,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
        suppress_callback_exceptions=True)

    select_x_options = [
        {"label": column, "value": column}
        for column in DATA.columns
    ]
    select_z_options = [
        {"label": column, "value": column}
        for column in ['None'] + list(DATA.columns)
    ]

    select_colormap_options = [{'label': val, 'value': val} for val in ['Viridis', 'Hot', 'Inferno', 'Jet',
                                                                        'Rainbow', 'RdBu', 'Discrete']]
    x_selection = dbc.FormGroup(
        [
            dbc.Col(dbc.Label("Select X"), style={
                    "text-align": "right"}, width=5),
            dbc.Col(
                dbc.Select(
                    id="selectX",
                    options=select_x_options,
                    value='S1DtoA'
                ), width=7
            )
        ], row=True
    )
    y_selection = dbc.FormGroup(
        [
            dbc.Col(dbc.Label("Select Y"), style={
                    "text-align": "right"}, width=5),
            dbc.Col(
                dbc.Select(
                    id="selectY",
                    options=select_x_options,
                    value='T1AtoA'
                ), width=7
            )
        ], row=True
    )
    z_selection = dbc.FormGroup(
        [
            dbc.Col(dbc.Label("Select Z"), style={
                    "text-align": "right"}, width=5),
            dbc.Col(
                dbc.Select(
                    id="selectZ",
                    options=select_z_options,
                    value='None'
                ), width=7
            )
        ], row=True
    )

    select_size_options = [
        {"label": column, "value": column}
        for column in ['Fix value'] + list(DATA.columns)
    ]
    size_selection = dbc.FormGroup(
        [
            dbc.Col(dbc.Label("Size"), style={
                    "text-align": "right"}, width=5),
            dbc.Col(
                dbc.Select(
                    id="selectSize",
                    options=select_size_options,
                    value='Fix value'
                ), width=7
            )
        ], row=True
    )
    color_selection = dbc.FormGroup(
        [
            dbc.Col(dbc.Label("Select color"), style={
                    "text-align": "right"}, width=5),
            dbc.Col(
                dbc.Select(
                    id="selectColor",
                    options=select_x_options,
                    value='S1_T1_split'
                ), width=7
            )
        ], row=True
    )
    colormap_selection = dbc.FormGroup(
        [
            dbc.Col(dbc.Label("Color scale"), style={
                    "text-align": "right"}, width=5),
            dbc.Col(
                dbc.Select(
                    id="selectColorMap",
                    options=select_colormap_options,
                    value='Discrete'
                ), width=7
            )
        ], row=True
    )

    app.layout = html.Div(children=[
        # html.Div(id='JmolDiv'),
        dcc.Store(id='memory'),
        html.H1(
            children='Donor-Acceptor Systems for Intramolecular Singlet Fission',
            style={
                "margin-left": "10pt",
                "margin-top": "10pt",
                "text-align": "center",
                "margin-bottom": "0pt"
            }),
        html.Div([
            dbc.Row([
                dbc.Col([
                    dbc.Button(
                        "Show data table", id="showTable", color="primary",
                        className="mr-1", style={"margin-bottom": "10pt"}),
                    html.Div(
                        [
                            dbc.Button("Information", id="open",
                                       style={"margin-bottom": "10pt"}),
                            dbc.Modal(
                                [
                                    dbc.ModalHeader("Information"),
                                    dbc.ModalBody(info_text, style={'padding-left': '50pt',
                                                                    'padding-right': '50pt'}),
                                    dbc.ModalFooter(
                                        dbc.Button("Close", id="close",
                                                   className="ml-auto")
                                    ),
                                ],
                                id="modal", size="xl",
                            ),
                        ], style={'display': 'inline'}
                    ),
                    html.Div(
                        id="tableContainer",
                        style={
                            "overflow-x": "auto",
                            "border": "1pt solid #d5d5d5",
                            "border-radius": "10pt",
                        }
                    ),
                ], width=12)
            ], style={"padding": "40pt", "padding-top": "0pt", "padding-bottom": "0pt"}),
            dbc.Row([
                dbc.Col([
                    dbc.Row(id="scatterContainer", style={
                            "padding": "40pt", "padding-top": "0pt"}),
                    dbc.Row([
                        dbc.Col(x_selection, width=4),
                        dbc.Col(y_selection, width=4),
                        dbc.Col(z_selection, width=4)
                    ], style={"padding-left": "40pt", "padding-right": "40pt"}),
                    dbc.Row([
                        dbc.Col(size_selection, width=4),
                        dbc.Col(color_selection, width=4),
                        dbc.Col(colormap_selection, width=4)
                    ], style={"padding-left": "40pt", "padding-right": "40pt"}),
                    dbc.Row([daq.Slider(
                        min=0.01,
                        max=2,
                        value=0.5,
                        handleLabel={"showCurrentValue": False, "label": "Scale"},
                        step=0.02,
                        id="size_scaler"
                    )
                    ], style={"padding-left": "80pt", "padding-right": "40pt", "padding-top": "10pt"})
                ], width=6, style={"padding": "10pt"}),
                dbc.Col([
                    dbc.Row([
                        dbc.Card([
                            html.H4(
                                "Details",
                                className="card-title",
                                style={"margin": "auto", "margin-bottom": "5pt", "margin-top": "5pt"}),
                            dbc.CardBody(
                                [
                                    html.P(
                                        "Select a datapoint to display the details")
                                ],
                                id="detailsContainer",
                                className="mb-3",
                                style={
                                    "max-height": "200pt",
                                    "width": "800px",
                                    "border": "1pt solid #d5d5d5",
                                    "border-radius": "10pt",
                                    "margin": "auto",
                                    "overflow-y": "auto",
                                    "overflow-x": "auto",
                                }
                            ),
                            html.H4(
                                "Molecular view",
                                className="card-title",
                                style={"margin": "auto", "margin-bottom": "5pt"}),
                            dbc.CardBody(children=[
                                html.P(
                                    "Select a datapoint to view a molecule")
                            ],
                                id="moleculeContainer",
                                className="mb-3",
                                style={
                                    # "max-width": "100%",
                                    "overflow-x": "auto",
                                    "width": "494px",
                                    "border": "1pt solid #d5d5d5",
                                    "border-radius": "10pt",
                                    "margin": "auto"
                            }
                            ),
                        ], style={"width": "100%", "margin-right": "60pt", "border": "0pt"}),
                    ]),
                ], width=6)
            ])
        ], id="container")
    ])


    @app.callback(
        [
            Output("scatterContainer", "children")
        ],
        [
            Input("selectX", "value"),
            Input("selectY", "value"),
            Input("selectZ", "value"),
            Input("selectSize", "value"),
            Input("selectColor", "value"),
            Input("selectColorMap", "value"),
            Input("size_scaler", "value"),

        ],
        [
            State("selectX", "value"),
            State("selectY", "value"),
            State("selectZ", "value"),
            State("selectSize", "value"),
            State("selectColor", "value"),
            State("selectColorMap", "value"),
            State("size_scaler", "value"),
        ]

    )
    def update_scatter_plot(input_x_column, input_y_column, input_z_column,
                            input_size_column, input_color_column, input_colormap,
                            input_size_scale,
                            x_column, y_column, z_column,
                            size_column, color_column, colormap_column,
                            size_scale):
        ctx = dash.callback_context
        if not ctx.triggered:
            button_id = 'No clicks yet'
        else:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]

        if button_id == "selectX":
            x_column = input_x_column
        elif button_id == "selectY":
            y_column = input_y_column
        elif button_id == "selectZ":
            z_column = input_z_column
        elif button_id == "selectSize":
            size_column = input_size_column
        elif button_id == "selectColor":
            color_column = input_color_column
        elif button_id == "selectColorMap":
            colormap_column = input_colormap
        elif button_id == "size_scale":
            size_scale = input_size_scale
        if size_column != 'Fix value':
            size_data = DATA[size_column].values.copy()
            size_data -= np.min(size_data) - 0.4
            size_data = (size_data / size_data.max())
            size_data *= size_scale
        else:
            size_data = np.ones(len(DATA))
        if colormap_column == 'Discrete':
            c1 = np.round(
                np.array([0.75361062, 0.83023285, 0.96087116, 1.]) * 255, 0).astype(int)
            c2 = np.round(
                np.array([0.61931795, 0.74412073, 0.99893092, 1.]) * 255, 0).astype(int)
            c3 = np.round(
                np.array([0.34832334, 0.46571115, 0.88834616, 1.]) * 255, 0).astype(int)
            c4 = np.round(
                np.array([0.70567316, 0.01555616, 0.15023281, 1.]) * 255, 0).astype(int)
            c5 = np.round(
                np.array([0.83936494, 0.32185622, 0.26492398, 1.]) * 255, 0).astype(int)
            c6 = np.round(
                np.array([0.96849975, 0.67397738, 0.55664926, 1.]) * 255, 0).astype(int)

            colors = [c1, c2, c3, c4, c5, c6]
            colors = sum([[ccc] * 2 for ccc in colors], [])
            n_bins = 7
            tiks = sum([[ccc] * 2 for ccc in np.linspace(0, 1, n_bins)], [])
            color_values = [[ttt, "rgb({},{},{})".format(
                ccc[0], ccc[1], ccc[2])] for ttt, ccc in zip(tiks[1:-1], colors)]

        if x_column and y_column and color_column:
            if z_column == 'None':
                fig = px.scatter(
                    DATA, x=x_column, y=y_column, color=color_column,
                    color_continuous_scale=color_values if colormap_column == 'Discrete' else colormap_column,
                    range_color=[DATA[color_column].min(),
                                 DATA[color_column].max()] if color_column != 'S1_T1_split' else [-2.5, 0.5],
                    opacity=0.8,
                    hover_data={x_column: True, y_column: True, color_column: True}
                )
                fig.update_layout(clickmode='event+select')

                if size_column == 'Fix value':
                    fig.update_traces(
                        marker_size=30 * size_scale,
                        marker_line_width=0, marker_line_color='DarkSlateGrey'
                    )
                else:
                    fig.update_traces(
                        marker_size=50 * size_data * size_scale,
                        marker_line={'width': 0},
                    )
                scatter_children = dcc.Graph(
                    id='mainScatter',
                    figure=fig,
                    style={"height": "500pt", "width": "100%"},
                    config={"scrollZoom": True},
                ),
            else:
                fig = px.scatter_3d(
                    DATA, x=x_column, y=y_column, z=z_column, color=color_column,
                    opacity=1,
                    color_continuous_scale=color_values if colormap_column == 'Discrete' else colormap_column,
                    range_color=[DATA[color_column].min(),
                                 DATA[color_column].max()] if color_column != 'S1_T1_split' else [-2.5, 0.5],
                    hover_data={x_column, y_column, color_column}
                )
                fig.update_layout(clickmode='event+select')
                if size_column == 'Fix value':
                    fig.update_traces(
                        marker_size=10 * size_scale,
                        marker_line={'width': 0},
                    )
                else:
                    fig.update_traces(
                        marker_size=50 * size_data * size_scale,
                        marker_line={'width': 0},
                    )

                scatter_children = dcc.Graph(
                    id='mainScatter',
                    figure=fig,
                    style={"height": "500pt", "width": "100%"},
                    config={"scrollZoom": True},
                )
        else:
            scatter_children = html.P(
                "Select X, Y and color columns to display the plot")
        return [
            scatter_children
        ]


    @app.callback(
        [
            Output('mainScatter', 'clickData'),
        ],
        [
            # Input('table', 'selected_rows'),derived_virtual_selected_row_ids
            # Input('table', 'derived_viewport_selected_rows'),
            # Input('table', 'derived_virtual_selected_rows'),
            Input('table', 'derived_viewport_selected_row_ids'),
        ],
        [
            State('selectX', 'value'),
            State('selectY', 'value'),
            # State('selectColor', 'value'),
            # State('mainScatter', 'figure'),
        ]
    )
    def select_point(
        selected_rows, xcol, ycol
    ):
        if xcol is not None and ycol is not None:
            if not selected_rows:
                selected_rows = [0]

            return [{
                    'points': [{
                        'pointNumber': selected_rows[0],
                        'pointIndex': selected_rows[0],
                        # 'x': pdata[xcol],
                        # 'y': pdata[ycol],
                        # 'marker.color': pdata[colorcol],
                    }
                    ]
                    }]


    @app.callback(
        [
            Output("detailsContainer", "children"),
            Output("moleculeContainer", "children"),
            Output("memory", "data")
        ],
        [
            Input("mainScatter", "clickData"),
        ]
    )
    def update_point_details(click_data):
        memory = {}
        details = html.P("Click on a point to see the details")
        vis_tabs = html.P("Click on a point to see the 3D molecular structure")
        if click_data:
            cdata = click_data["points"][0]
            if "pointIndex" in cdata.keys():
                memory["point"] = click_data["points"][0]["pointIndex"]
                pdata = DATA.iloc[click_data["points"][0]["pointIndex"]].to_dict()
            else:
                memory["point"] = click_data["points"][0]["pointNumber"]
                pdata = DATA.iloc[click_data["points"][0]["pointNumber"]].to_dict()

            details = [
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**Compound index**: {}'''.format(pdata['COMP_names'])],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**S1<sub>exc</sub>**: {0:.3f} eV &emsp;&emsp;'''.format(pdata[
                                     'S1_exc']),
                                 '''&emsp; **T1<sub>exc</sub>**: {0:.3f} eV&emsp;&emsp;'''.format(pdata[
                                     'T1_exc']),
                                 '''**&Delta;E<sub>ST</sub>**: {0:.3f} eV'''.format(pdata['S1_T1_split'])],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**&Omega;<sup>S<sub>1</sub></sup><sub>A&#8594;A</sub>**: {0:.3f}&emsp;'''.format(pdata[
                                     'S1AtoA']),
                                 '''**&Omega;<sup>S<sub>1</sub></sup><sub>A&#8594;D</sub>**: {0:.3f}&emsp;'''.format(pdata[
                                     'S1AtoD']),
                                 '''**&Omega;<sup>S<sub>1</sub></sup><sub>D&#8594;A</sub>**: {0:.3f}&emsp;'''.format(pdata[
                                     'S1DtoA']),
                                 '''**&Omega;<sup>S<sub>1</sub></sup><sub>D&#8594;D</sub>**: {0:.3f}'''.format(pdata[
                                     'S1DtoD']),
                             ], style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**&Omega;<sup>T<sub>1</sub></sup><sub>A&#8594;A</sub>**: {0:.3f}&emsp;'''.format(pdata[
                                     'T1AtoA']),
                                 '''**&Omega;<sup>T<sub>1</sub></sup><sub>A&#8594;D</sub>**: {0:.3f}&emsp;'''.format(pdata[
                                     'T1AtoD']),
                                 '''**&Omega;<sup>T<sub>1</sub></sup><sub>D&#8594;A</sub>**: {0:.3f}&emsp;'''.format(pdata[
                                     'T1DtoA']),
                                 '''**&Omega;<sup>T<sub>1</sub></sup><sub>D&#8594;D</sub>**: {0:.3f}'''.format(pdata[
                                     'S1DtoD']),
                             ], style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**D<sub>HOMO</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(
                                     pdata['D_HOMO']),
                                 '''**D<sub>LUMO</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(
                                     pdata['D_LUMO']),
                                 '''**D<sub>GAP</sub> **: {0:.3f} eV'''.format(
                                     pdata['D_gap']),

                             ],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**A<sub>HOMO</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(
                                     pdata['A_HOMO']),
                                 '''**A<sub>LUMO</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(
                                     pdata['A_LUMO']),
                                 '''**A<sub>GAP</sub> **: {0:.3f} eV'''.format(
                                     pdata['A_gap']),

                             ],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**D<sub>S1</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(pdata[
                                     'D_S1']),
                                 '''**D<sub>T1</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(pdata[
                                     'D_T1']),
                                 '''**D(&Delta;E<sub>ST</sub>)**: {0:.3f} eV'''.format(
                                     pdata['D_split']),

                             ],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**A<sub>S1</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(pdata[
                                     'A_S1']),
                                 '''**A<sub>T1</sub> **: {0:.3f} eV &emsp;&emsp;'''.format(pdata[
                                     'A_T1']),
                                 '''**A(&Delta;E<sub>ST</sub>)**: {0:.3f} eV'''.format(
                                     pdata['A_split']),

                             ],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                dcc.Markdown(dangerously_allow_html=True,
                             children=[
                                 '''**&phi;<sub>D-A</sub>**: {0:.3f}&deg; '''.format(
                                     pdata['dft_dihedral_norm']),
                             ],
                             style={"margin-bottom": "-12pt", "fontsize": "10%"}),
                html.P([html.Strong("SMILES: "), "{}".format(
                    pdata['SMILES'])], style={"margin-bottom": "0pt"}),
            ]

            vis_tabs = dbc.Card(
                [
                    dbc.CardHeader(
                        dbc.Tabs(
                            [
                                dbc.Tab(label="2D", tab_id="tab-jmol"),
                                dbc.Tab(label="3D", tab_id="tab-speck"),
                            ],
                            id="vis-tabs",
                            card=True,
                            active_tab="tab-speck",
                        )
                    ),
                    dbc.CardBody(id="vis-content", style={"padding": 0}),
                ]
            )

        return [

            details, vis_tabs, memory
        ]


    @app.callback(
        Output("vis-content", "children"),
        [
            Input("vis-tabs", "active_tab")
        ],
        State("memory", "data")
    )
    def tab_content(active_tab, memory):
        if "point" not in memory:
            raise PreventUpdate

        mol_name = DATA.iloc[memory["point"]].COMP_names
        path = './structures/' + mol_name + '.xyz'

        if active_tab == "tab-speck":

            mol = aio.read(path)
            mol_data = [{'symbol': a, 'x': xyz[0],
                         'y': xyz[1], 'z': xyz[2]} for a, xyz in zip(
                mol.get_chemical_symbols(), mol.positions)]

            mol_plot = dashbio.Speck(
                id='my-speck', data=mol_data,
                view={
                    'zoom': 0.06,
                    'resolution': 450,
                    'ao': 0.0001,
                    # 'outline': 1,
                    'atomScale': 0.15,
                    'relativeAtomScale': 0.33,
                    'bonds': True
                },
                presetView='stickball',
            )
        else:
            svg_file = './images/' + mol_name + '.svg'
            encoded = base64.b64encode(open(svg_file, 'rb').read())
            svg = 'data:image/svg+xml;base64,{}'.format(encoded.decode())
            mol_plot = html.Img(src=svg, style={'width': '100%'})
        return mol_plot


    @app.callback(
        Output("modal", "is_open"),
        [Input("open", "n_clicks"), Input("close", "n_clicks")],
        [State("modal", "is_open")],
    )
    def toggle_modal(n1, n2, is_open):
        if n1 or n2:
            return not is_open
        return is_open


    @app.callback(
        [
            Output("tableContainer", "children"),
            Output("showTable", "children")
        ],
        [
            Input("showTable", "n_clicks"),
        ],
        [
            State("showTable", "children")
        ]
    )
    def show_table(bt, text):
        if "show" in text.lower() and bt:
            children = [
                dash_table.DataTable(
                    id='table',
                    columns=[{'name': i, 'id': i}
                             for i in DATA.columns if i != 'id'],
                    data=DATA.head(10).to_dict("rows"),
                    style_cell={
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'maxWidth': '300px',
                        'minWidth': '100px'
                    },
                    sort_action="custom",
                    sort_mode="multi",
                    row_selectable='single',
                    page_action="custom",
                    selected_rows=[],
                    page_current=0,
                    page_size=10,
                    sort_by=[],
                    style_data_conditional=[{
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(248, 248, 248)'
                    }],
                    style_header={
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'fontWeight': 'bold'
                    },
                    style_table={
                        'padding-left': '11pt',
                        'padding-right': '20pt',
                        'padding-top': '5pt',
                    },
                    css=[
                        {
                            'selector': 'dash-fixed-content',
                            'rule': 'height: 100%;'
                        },
                    ],
                )
            ]
            button_text = "Hide data table"
        else:
            children = []
            button_text = "Show data table"

        return [children, button_text]


    @app.callback(
        [
            Output('table', 'data'),
            Output('table', 'columns'),
            Output('table', 'editable'),
            Output('table', 'row_deletable'),
            Output('table', 'page_count')
        ],
        [
            Input('table', 'page_size'),
            Input('table', 'page_current'),
            Input('table', 'data_timestamp'),
            Input('table', 'sort_by'),
        ],
        [
            State("table", "data")
        ]
    )
    def update_table(page_size, page_current, data_timestamp, sort_by, data):
        if sort_by and len(sort_by):
            result_sorted = DATA.sort_values(
                [col['column_id'] for col in sort_by],
                ascending=[
                    col['direction'] == 'asc'
                    for col in sort_by
                ],
                inplace=False
            )
        else:
            result_sorted = DATA

        result_paginated = result_sorted.iloc[
            page_current * page_size:(page_current + 1) * page_size
        ]

        page_count = len(result_sorted) // page_size

        columns = [
            {
                "name": i,
                "id": i,
                "clearable": True,
                "selectable": True,
                "renamable": False,
                "hideable": True,
                "deletable": False
            }
            for i in DATA.columns
        ]

        return [
            result_paginated.to_dict('records'),
            columns, True, False, page_count
        ]

    return app


def create_server():
    """Shorthand function for running through gunicorn."""
    return create_app().server
