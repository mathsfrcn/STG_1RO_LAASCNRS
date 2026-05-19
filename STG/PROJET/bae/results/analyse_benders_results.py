import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter as tk
from tkinter import filedialog

def analyze_benders_results(filepath):
    try:
        df = pd.read_csv(filepath, sep=r'[,\s]+', engine='python', 
                         names=['Method', 'Gamma', 'Tau', 'Iterations', 'Time'])
    except Exception as e:
        print("Erreur lors de la lecture du fichier :", e)
        return

    pivot_df = df.pivot_table(index=['Gamma', 'Tau'], 
                              columns='Method', 
                              values=['Iterations', 'Time']).reset_index()
    
    pivot_df.columns = [f"{col[0]}_{col[1]}" if col[1] else col[0] for col in pivot_df.columns]

    # If > 0, KC is the best
    pivot_df['Time_Saved_by_KC'] = pivot_df['Time_STANDARD'] - pivot_df['Time_KC']
    pivot_df['Iter_Saved_by_KC'] = pivot_df['Iterations_STANDARD'] - pivot_df['Iterations_KC']

    total_instances = len(pivot_df)
    kc_faster_count = (pivot_df['Time_Saved_by_KC'] > 0).sum()
    kc_fewer_iters_count = (pivot_df['Iter_Saved_by_KC'] > 0).sum()
    avg_speedup = (pivot_df['Time_STANDARD'] / pivot_df['Time_KC']).mean()

    print("\n" + "="*40)
    print("RAPPORT D'ANALYSE DE PERFORMANCE")
    print("="*40)
    print(f"Total de configurations (\u0393, \u03C4) testées : {total_instances}")
    print(f"KC plus rapide sur : {kc_faster_count}/{total_instances} instances ({(kc_faster_count/total_instances)*100:.1f}%)")
    print(f"KC fait moins d'itérations sur : {kc_fewer_iters_count}/{total_instances} instances ({(kc_fewer_iters_count/total_instances)*100:.1f}%)")
    print("-" * 40)
    print(f"Temps moyen STANDARD : {pivot_df['Time_STANDARD'].mean():.5f} s")
    print(f"Temps moyen KC       : {pivot_df['Time_KC'].mean():.5f} s")
    print(f"Speedup moyen (Std/KC): x{avg_speedup:.2f}")
    print("="*40 + "\n")

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Analyse Comparative : Décomposition Standard vs Heuristique KC', fontsize=18, fontweight='bold')

    ax = axes[0, 0]
    max_time = max(pivot_df['Time_STANDARD'].max(), pivot_df['Time_KC'].max()) * 1.1
    sns.scatterplot(data=pivot_df, x='Time_STANDARD', y='Time_KC', hue='Gamma', 
                    size='Tau', sizes=(20, 200), palette='viridis', ax=ax, alpha=0.8)
    
    ax.plot([0, max_time], [0, max_time], 'r--', label='Égalité (Temps Std = Temps KC)')
    ax.fill_between([0, max_time], [0, max_time], [0, 0], color='green', alpha=0.1, label='Zone de victoire KC')
    ax.set_title("Temps de Résolution Global (s)")
    ax.set_xlim(0, max_time)
    ax.set_ylim(0, max_time)
    ax.legend(title='Budget \u0393')

    ax = axes[0, 1]
    df_iter = df.groupby(['Tau', 'Method'])['Iterations'].mean().reset_index()
    sns.lineplot(data=df_iter, x='Tau', y='Iterations', hue='Method', marker='o', 
                 palette={'STANDARD': '#d62728', 'KC': '#2ca02c'}, ax=ax, linewidth=2.5, markersize=8)
    ax.set_title("Impact du coefficient d'approximation (\u03C4) sur les Itérations")
    ax.set_ylabel("Nombre moyen d'itérations")

    ax = axes[1, 0]
    heat_time = pivot_df.pivot(index="Gamma", columns="Tau", values="Time_Saved_by_KC")
    sns.heatmap(heat_time, cmap="RdYlGn", center=0, annot=True, fmt=".4f", 
                cbar_kws={'label': 'Temps économisé par KC (s)'}, ax=ax)
    ax.set_title("Où KC est-il plus rapide ? (Matrice Gain de Temps)")

    ax = axes[1, 1]
    heat_iter = pivot_df.pivot(index="Gamma", columns="Tau", values="Iter_Saved_by_KC")
    sns.heatmap(heat_iter, cmap="RdYlGn", center=0, annot=True, fmt=".1f", 
                cbar_kws={'label': 'Itérations économisées par KC'}, ax=ax)
    ax.set_title("Où KC réduit-il les coupes ? (Matrice Gain Itérations)")

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    
    file_path = filedialog.askopenfilename(
        title="Sélectionnez le fichier de résultats :",
        filetypes=[("Données", "*.csv *.txt"), ("Tous", "*.*")]
    )
    
    if file_path:
        analyze_benders_results(file_path)
    else:
        print("Aucun fichier sélectionné. Arrêt de l'analyse.")