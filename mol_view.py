# -*- coding: utf-8 -*-

# Run this app with `python mol_viewer.py` and
# visit http://127.0.0.1:8060/ in your web browser.

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


data_file = sys.argv[1]
structures_dir = sys.argv[2]

DATA = pd.read_csv(
    data_file)

DATA['id'] = list(DATA.index)
# DATA = DATA.round(6)

# dtypes = DATA.dtypes

# # round floats to N significant decimals
# N = 4
# for i, col in enumerate(DATA.columns):
#     if dtypes[i] == float:
#         DATA[col] = DATA[col].apply(lambda x: np.round(x, N - int(np.floor(np.log10(abs(x))))))

# pd.options.display.float_format = '${:.2f}'.format

info_text = [dcc.Markdown('''
    # App instructions
    - Click once in the scatter plot to select a point and see the corresponding molecular structure, click twice to deselect it.
    - Select a row in the data table to see the corresponding molecular structure.
    - Change datapoint size with the sliding scale below the plot.
    - Select atoms in the 2D view to compute distances, angles and dihedrals of the corresponding 3D geometry.
    ''')]

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True)

select_x_options = [
    {'label': column, 'value': column}
    for column in DATA.columns
]
select_z_options = [
    {'label': column, 'value': column}
    for column in ['None'] + list(DATA.columns)
]

select_colormap_options = [{'label': val, 'value': val} for val in ['Viridis', 'Hot', 'Inferno', 'Jet',
                                                                    'Rainbow', 'RdBu']]
x_selection = dbc.FormGroup(
    [
        dbc.Col(dbc.Label('Select X'), style={
                'text-align': 'right'}, width=5),
        dbc.Col(
            dbc.Select(
                id='selectX',
                options=select_x_options,
                value=DATA.columns[0]
            ), width=7
        )
    ], row=True
)
y_selection = dbc.FormGroup(
    [
        dbc.Col(dbc.Label('Select Y'), style={
                'text-align': 'right'}, width=5),
        dbc.Col(
            dbc.Select(
                id='selectY',
                options=select_x_options,
                value=DATA.columns[1]
            ), width=7
        )
    ], row=True
)
z_selection = dbc.FormGroup(
    [
        dbc.Col(dbc.Label('Select Z'), style={
                'text-align': 'right'}, width=5),
        dbc.Col(
            dbc.Select(
                id='selectZ',
                options=select_z_options,
                value='None'
            ), width=7
        )
    ], row=True
)

select_size_options = [
    {'label': column, 'value': column}
    for column in ['Fix value'] + list(DATA.columns)
]
size_selection = dbc.FormGroup(
    [
        dbc.Col(dbc.Label('Size'), style={
                'text-align': 'right'}, width=5),
        dbc.Col(
            dbc.Select(
                id='selectSize',
                options=select_size_options,
                value='Fix value'
            ), width=7
        )
    ], row=True
)
color_selection = dbc.FormGroup(
    [
        dbc.Col(dbc.Label('Select color'), style={
                'text-align': 'right'}, width=5),
        dbc.Col(
            dbc.Select(
                id='selectColor',
                options=select_x_options,
                value=DATA.columns[3]
            ), width=7
        )
    ], row=True
)
colormap_selection = dbc.FormGroup(
    [
        dbc.Col(dbc.Label('Color scale'), style={
                'text-align': 'right'}, width=5),
        dbc.Col(
            dbc.Select(
                id='selectColorMap',
                options=select_colormap_options,
                value='Viridis'
            ), width=7
        )
    ], row=True
)

app.layout = html.Div(children=[
    # html.Div(id='JmolDiv'),
    dcc.Store(id='memory'),
    html.H1(
        children='Molecular data explorer',
        style={
            'margin-left': '10pt',
            'margin-top': '10pt',
            'text-align': 'center',
            'margin-bottom': '20pt'
        }),
    html.Div([
        dbc.Row([
            dbc.Col([
                dbc.Button(
                    'Show data table', id='showTable', color='primary',
                    className='mr-1', style={'margin-bottom': '10pt'}),
                html.Div(
                    [
                        dbc.Button('Information', id='open',
                                   style={'margin-bottom': '10pt'}),
                        dbc.Modal(
                            [
                                dbc.ModalHeader('Information'),
                                dbc.ModalBody(info_text, style={'padding-left': '50pt',
                                                                'padding-right': '50pt'}),
                                dbc.ModalFooter(
                                    dbc.Button('Close', id='close',
                                               className='ml-auto')
                                ),
                            ],
                            id='modal', size='xl',
                        ),
                    ], style={'display': 'inline'}
                ),
                html.Div(
                    id='tableContainer',
                    style={
                        'overflow-x': 'auto',
                        'border': '1pt solid #d5d5d5',
                        'border-radius': '10pt',
                    }
                ),
            ], width=12)
        ], style={'padding': '40pt', 'padding-top': '0pt', 'padding-bottom': '0pt'}),
        dbc.Row([
            dbc.Col([
                dbc.Row(id='scatterContainer', style={
                        'padding': '40pt', 'padding-top': '0pt'}),
                dbc.Row([
                    dbc.Col(x_selection, width=4),
                    dbc.Col(y_selection, width=4),
                    dbc.Col(z_selection, width=4)
                ], style={'padding-left': '40pt', 'padding-right': '40pt'}),
                dbc.Row([
                    dbc.Col(size_selection, width=4),
                    dbc.Col(color_selection, width=4),
                    dbc.Col(colormap_selection, width=4)
                ], style={'padding-left': '40pt', 'padding-right': '40pt'}),
                dbc.Row([daq.Slider(
                    min=0.01,
                    max=2,
                    value=0.5,
                    handleLabel={'showCurrentValue': False, 'label': 'Scale'},
                    step=0.02,
                    id='size_scaler'
                )
                ], style={'padding-left': '80pt', 'padding-right': '40pt', 'padding-top': '10pt'})
            ], width=6, style={'padding': '10pt'}),
            dbc.Col([
                dbc.Row([
                    dbc.Card([
                        html.H4(
                            'Details',
                            className='card-title',
                            style={'margin': 'auto', 'margin-bottom': '5pt', 'margin-top': '5pt'}),
                        dbc.CardBody(
                            [
                                html.P(
                                    'Select a datapoint to display the details')
                            ],
                            id='detailsContainer',
                            className='mb-3',
                            style={
                                'max-height': '200pt',
                                'width': '800px',
                                'border': '1pt solid #d5d5d5',
                                'border-radius': '10pt',
                                'margin': 'auto',
                                'overflow-y': 'auto',
                                'overflow-x': 'auto',
                            }
                        ),
                        html.H4(
                            'Molecular view',
                            className='card-title',
                            style={'margin': 'auto', 'margin-bottom': '5pt'}),
                        dbc.CardBody(children=[
                            html.P(
                                'Select a datapoint to view a molecule')
                        ],
                            id='moleculeContainer',
                            className='mb-3',
                            style={
                                # 'max-width': '100%',
                                'overflow-x': 'auto',
                                'width': '494px',
                                'border': '1pt solid #d5d5d5',
                                'border-radius': '10pt',
                                'margin': 'auto'
                        }
                        ),
                    ], style={'width': '100%', 'margin-right': '60pt', 'border': '0pt'}),
                ]),
            ], width=6)
        ])
    ], id='container')
])


@app.callback(
    [
        Output('scatterContainer', 'children')
    ],
    [
        Input('selectX', 'value'),
        Input('selectY', 'value'),
        Input('selectZ', 'value'),
        Input('selectSize', 'value'),
        Input('selectColor', 'value'),
        Input('selectColorMap', 'value'),
        Input('size_scaler', 'value'),

    ],
    [
        State('selectX', 'value'),
        State('selectY', 'value'),
        State('selectZ', 'value'),
        State('selectSize', 'value'),
        State('selectColor', 'value'),
        State('selectColorMap', 'value'),
        State('size_scaler', 'value'),
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

    if button_id == 'selectX':
        x_column = input_x_column
    elif button_id == 'selectY':
        y_column = input_y_column
    elif button_id == 'selectZ':
        z_column = input_z_column
    elif button_id == 'selectSize':
        size_column = input_size_column
    elif button_id == 'selectColor':
        color_column = input_color_column
    elif button_id == 'selectColorMap':
        colormap_column = input_colormap
    elif button_id == 'size_scale':
        size_scale = input_size_scale
    if size_column != 'Fix value':
        size_data = DATA[size_column].values.copy()
        size_data -= np.min(size_data) - 0.4
        size_data = (size_data / size_data.max())
        size_data *= size_scale
    else:
        size_data = np.ones(len(DATA))

    if x_column and y_column and color_column:
        if z_column == 'None':
            fig = px.scatter(
                DATA, x=x_column, y=y_column, color=color_column,
                color_continuous_scale=colormap_column,
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
                style={'height': '500pt', 'width': '100%'},
                config={'scrollZoom': True},
            ),
        else:
            fig = px.scatter_3d(
                DATA, x=x_column, y=y_column, z=z_column, color=color_column,
                opacity=1,
                color_continuous_scale=colormap_column,
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
                style={'height': '500pt', 'width': '100%'},
                config={'scrollZoom': True},
            )
    else:
        scatter_children = html.P(
            'Select X, Y and color columns to display the plot')
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
        Output('detailsContainer', 'children'),
        Output('moleculeContainer', 'children'),
        Output('memory', 'data')
    ],
    [
        Input('mainScatter', 'clickData'),
    ]
)
def update_point_details(click_data):

    memory = {}
    details = html.P('Click on a point to see the details')
    vis_tabs = html.P('Click on a point to see the 3D molecular structure')
    if click_data:
        cdata = click_data['points'][0]
        if 'pointIndex' in cdata.keys():
            memory['point'] = click_data['points'][0]['pointIndex']
            pdata = DATA.iloc[click_data['points'][0]['pointIndex']].to_dict()
        else:
            memory['point'] = click_data['points'][0]['pointNumber']
            pdata = DATA.iloc[click_data['points'][0]['pointNumber']].to_dict()

        details = [
            html.P([html.Strong('{}: '.format(col)),
                    '{}'.format('{:.2E}'.format(pdata[col]) if isinstance(
                        pdata[col], float) else pdata[col])],
                   style={'margin-bottom': '0pt'})
            for col in DATA.columns
        ]
        # details = [dcc.Markdown(dangerously_allow_html=True,
        #                         children=[
        #                             '''**{}**: {}'''.format(col, ],
        #                         style={'margin-bottom': '-12pt',
        #                                'fontsize': '10%'}) for col in DATA.columns]
        # details = [
        #     html.P('{}: {}'.format(k, v), style={'margin-bottom': '0pt'})
        #     for k, v in DATA.iloc[
        #         click_data['points'][0]['pointIndex']].to_dict().items()]

        vis_tabs = dbc.Card(
            [
                dbc.CardHeader(
                    dbc.Tabs(
                        [
                            dbc.Tab(label='2D', tab_id='tab-2d'),
                            dbc.Tab(label='3D', tab_id='tab-speck'),
                        ],
                        id='vis-tabs',
                        card=True,
                        active_tab='tab-speck',
                    )
                ),
                dbc.CardBody(id='vis-content', style={'padding': 0}),
            ]
        )

    return [

        details, vis_tabs, memory
    ]


@app.callback(
    Output('vis-content', 'children'),
    [
        Input('vis-tabs', 'active_tab')
    ],
    State('memory', 'data')
)
def tab_content(active_tab, memory):
    if 'point' not in memory:
        raise PreventUpdate

    mol_name = DATA.iloc[memory['point']].xyz_file_name
    path = structures_dir + '/' + mol_name

    if active_tab == 'tab-speck':

        mol = aio.read(path)
        mol_data = [{'symbol': a, 'x': xyz[0],
                     'y': xyz[1], 'z': xyz[2]} for a, xyz in zip(
            mol.get_chemical_symbols(), mol.positions)]

        mol_plot = dashbio.Speck(
            id='my-speck', data=mol_data,
            view={
                'zoom': 0.06,
                'resolution': 450,
                # 'ao': 0.1,
                # 'outline': 1,
                'atomScale': 0.15,
                'relativeAtomScale': 0.33,
                'bonds': True
            },
            presetView='stickball',
        )
    else:
        mol = aio.read(path)
        model_data = xyz_to_json(mol)
        mol_plot = html.Div([
            dashbio.Molecule2dViewer(
                id='my-dashbio-molecule2d',
                modelData=model_data,
                width=450,
                height=400,
            ),
            # html.Hr(style={'padding-top': '0pt'}),
            html.Div(id='molecule2d-output')
        ])
    return mol_plot


@app.callback(
    Output('molecule2d-output', 'children'),
    [Input('my-dashbio-molecule2d', 'selectedAtomIds')],
    State('memory', 'data')
)
def update_selected_atoms(ids, memory):
    mol_name = DATA.iloc[memory['point']].xyz_file_name
    path = structures_dir + '/' + mol_name
    mol = aio.read(path)

    if ids is None or len(ids) == 0:
        return 'No atom has been selected. Select atoms by clicking on them.'
    elif len(ids) == 1:
        return 'Selected atom ID: {}.'.format(', '.join([str(i) for i in ids]))
    elif len(ids) == 2:
        dist = mol.get_distance(ids[0], ids[1])
        return 'Distance between atoms {}-{}: {:.3f}'.format(ids[0], ids[1], dist)
    elif len(ids) == 3:
        dist = mol.get_angle(ids[0], ids[1], ids[2])
        return 'Angle between atoms {}-{}-{}: {:.3f}'.format(ids[0], ids[1], ids[2], dist)
    elif len(ids) == 4:
        dist = mol.get_dihedral(ids[0], ids[1], ids[2], ids[3])
        return 'Dihedral angle between atoms {}-{}-{}-{}: {:.3f}'.format(ids[0], ids[1],
                                                                         ids[2], ids[
                                                                             3],
                                                                         dist)
    else:
        return 'Selected atom IDs: {}.'.format(', '.join([str(i) for i in ids]))

    # return 'Selected atom IDs: {}.'.format(', '.join([str(i) for i in ids]))


@app.callback(
    Output('modal', 'is_open'),
    [Input('open', 'n_clicks'), Input('close', 'n_clicks')],
    [State('modal', 'is_open')],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open


@app.callback(
    [
        Output('tableContainer', 'children'),
        Output('showTable', 'children')
    ],
    [
        Input('showTable', 'n_clicks'),
    ],
    [
        State('showTable', 'children')
    ]
)
def show_table(bt, text):
    if 'show' in text.lower() and bt:
        children = [
            dash_table.DataTable(
                id='table',
                columns=[{'name': i, 'id': i}
                         for i in DATA.columns if i != 'id'],
                data=DATA.head(10).to_dict('rows'),
                sort_action='custom',
                sort_mode='multi',
                row_selectable='single',
                page_action='custom',
                selected_rows=[],
                page_current=0,
                page_size=10,
                sort_by=[],
                style_cell={
                    'overflow': 'hidden',
                    'textOverflow': 'ellipsis',
                    'maxWidth': '300px',
                    'minWidth': '100px'
                },
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
        button_text = 'Hide data table'
    else:
        children = []
        button_text = 'Show data table'

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
        State('table', 'data')
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
            'name': i,
            'id': i,
            'clearable': True,
            'selectable': True,
            'renamable': False,
            'hideable': True,
            'deletable': False
        }
        for i in DATA.columns
    ]

    return [
        result_paginated.to_dict('records'),
        columns, True, False, page_count
    ]


app.run_server(port=8060, debug=True)
