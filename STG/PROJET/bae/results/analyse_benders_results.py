import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter as tk
import os
from tkinter import filedialog

def analyze_benders_results(file_path):
    try:
        df = pd.read_csv(file_path, sep=r'[,\s]+', engine='python', 
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

    file_name = os.path.basename(file_path)
    sns.set_theme(style='whitegrid')
    fig = plt.figure
    fig = plt.figure(figsize=(16, 12))
    fig.canvas.manager.set_window_title(f"Analyse Benders - {file_name}")
    axes = fig.subplots(2, 2)
    fig.suptitle(f'Analyse Comparative : {file_name}', fontsize=16, fontweight='bold')

    # Fig01
    ax = axes[0, 0]
    max_time = max(pivot_df['Time_STANDARD'].max(), pivot_df['Time_KC'].max()) * 1.1
    sns.scatterplot(data=pivot_df, x='Time_STANDARD', y='Time_KC', hue='Gamma', size='Tau', sizes=(20, 200), palette='viridis', ax=ax, alpha=0.8)
    ax.plot([0, max_time], [0, max_time], 'r--', label='Égalité')
    ax.fill_between([0, max_time], [0, max_time], [0, 0], color='green', alpha=0.1, label='Zone de victoire de KC')
    ax.set_title("Temps de Résolution Global (s)")
    ax.set_xlim(0, max_time)
    ax.set_ylim(0, max_time)
    ax.legend(title='Budget Gamma')

    # Fig02
    ax = axes[0, 1]
    df_iter = df.groupby(['Tau', 'Method'])['Iterations'].mean().reset_index()
    sns.lineplot(data=df_iter, x='Tau', y='Iterations', hue='Method', marker='o', palette={'STANDARD': '#d62728', 'KC': '#2ca02c'}, ax=ax)
    ax.set_title("Impact de tau sur les Itérations")

    # Fig03
    ax = axes[1, 0]
    heat_time = pivot_df.pivot(index="Gamma", columns="Tau", values="Time_Saved_by_KC")
    sns.heatmap(heat_time, cmap="RdYlGn", center=0, annot=True, fmt=".4f", ax=ax)
    ax.set_title("Matrice Gain de Temps (s)")

    # Fig04
    ax = axes[1, 1]
    heat_iter = pivot_df.pivot(index="Gamma", columns="Tau", values="Iter_Saved_by_KC")
    sns.heatmap(heat_iter, cmap="RdYlGn", center=0, annot=True, fmt=".1f", ax=ax)
    ax.set_title("Matrice Gain Itérations")

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    base_name_no_ext = os.path.splitext(file_name)[0]
    output_image_name = base_name_no_ext + "_analytics.png"
    target_directory = os.path.dirname(file_path)
    output_image_path = os.path.join(target_directory, output_image_name)

    plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
    plt.show()

def main():
    root = tk.Tk()
    root.withdraw()
    
    file_path = filedialog.askopenfilename(
        title="Sélectionnez le fichier de résultats",
        filetypes=[("Données", "*.csv *.txt"), ("Tous", "*.*")]
    )
    
    if file_path:
        analyze_benders_results(file_path)
    else:
        print("Aucun fichier sélectionné.")

if __name__ == "__main__":
    main()