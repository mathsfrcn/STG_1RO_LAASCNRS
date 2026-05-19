import json
import networkx as nx
import plotly.graph_objects as go
import tkinter as tk
from tkinter import filedialog

def visualize_benders_graph(json_file_path):
    with open(json_file_path, 'r') as f:
        data = json.load(f)
        
    G = nx.DiGraph()
    
    node_x = []
    node_y = []
    node_text = []
    
    for node in data['nodes']:
        G.add_node(node['id'], t=node['t'], b=node['b'], pi=node['pi'])
        node_x.append(node['t'])
        node_y.append(node['b'])
        node_text.append(f"t: {node['t']}<br>Budget: {node['b']}<br>Pi: {node['pi']:.2f}")

    edge_x = []
    edge_y = []
    worst_edge_x = []
    worst_edge_y = []
    
    for link in data['links']:
        source = link['source']
        target = link['target']
        G.add_edge(source, target, cost=link['cost'], is_worst=link['is_worst'])
        
        x0, y0 = G.nodes[source]['t'], G.nodes[source]['b']
        x1, y1 = G.nodes[target]['t'], G.nodes[target]['b']
        
        if link['is_worst'] == 1:
            worst_edge_x.extend([x0, x1, None])
            worst_edge_y.extend([y0, y1, None])
        else:
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

    edges_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines',
        name='Graphe Complet'
    )
    
    worst_edges_trace = go.Scatter(
        x=worst_edge_x, y=worst_edge_y,
        line=dict(width=3, color='red'),
        hoverinfo='none',
        mode='lines',
        name='Pire Scénario (Sous-graphe)'
    )

    nodes_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=node_text,
        marker=dict(
            showscale=True,
            colorscale='Viridis',
            color=[G.nodes[n]['pi'] for n in G.nodes()],
            size=10,
            colorbar=dict(thickness=15, title='Valeur Pi', xanchor='left')
        )
    )

    fig = go.Figure(data=[edges_trace, worst_edges_trace, nodes_trace],
             layout=go.Layout(
                title='Graphe de Budget Benders - Compilation de Connaissances',
                showlegend=True,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(title='Période (t)', showgrid=False, zeroline=False),
                yaxis=dict(title='Budget Consommé (b)', showgrid=True, zeroline=False)
             ))
    
    fig.show()

def main():
    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Sélectionnez l'instance Benders à visualiser",
        filetypes=[("Fichiers JSON", "*.json"), ("Tous les fichiers", "*.*")]
    )

    if not file_path:
        print("Aucun fichier sélectionné. Arrêt du programme.")
        return

    print(f"Génération du graphique pour : {file_path}")
    visualize_benders_graph(file_path)

if __name__ == "__main__":
    main()