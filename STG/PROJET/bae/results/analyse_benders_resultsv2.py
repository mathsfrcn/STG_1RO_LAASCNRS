import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def analyze_benders_results(file_path):
    try:
        df = pd.read_csv(file_path, sep=r'[,\s]+', engine='python', 
                         names=['Method', 'Gamma', 'Tau', 'Iterations', 'Time'])
    except Exception as e:
        print("Error reading file :", e)
        return

    methods = df["Method"].unique()

    if len(methods) != 2:
        print(f"Error : The file must contain at least 2 references.")
        return -1
    
    methods_sorted = sorted(methods, key=lambda x: 0 if x == 'STANDARD' else (1 if x == 'KC' else 2))
    ref_method = methods_sorted[0]
    alt_method = methods_sorted[1]

    pivot_df = df.pivot_table(index=['Gamma', 'Tau'], 
                              columns='Method', 
                              values=['Iterations', 'Time']).reset_index()
    
    pivot_df.columns = [f"{col[0]}_{col[1]}" if col[1] else col[0] for col in pivot_df.columns]

    pivot_df['Time_Saved'] = pivot_df[f'Time_{ref_method}'] - pivot_df[f'Time_{alt_method}']    # If > 0, KC is the best
    pivot_df['Iter_Saved'] = pivot_df[f'Iterations_{ref_method}'] - pivot_df[f'Iterations_{alt_method}']

    file_name = os.path.basename(file_path)
    sns.set_theme(style='whitegrid')
    fig = plt.figure
    fig = plt.figure(figsize=(16, 12))
    fig.canvas.manager.set_window_title(f"Benders analyse - {file_name}")
    axes = fig.subplots(2, 2)
    fig.suptitle(f'Comparative analysis : {ref_method} vs {alt_method}', fontsize=16, fontweight='bold')

    # Fig01
    ax = axes[0, 0]

    pivot_df['Speedup'] = pivot_df[f'Time_{ref_method}'] / pivot_df[f'Time_{alt_method}']
    
    sns.lineplot(data=pivot_df, x='Gamma', y='Speedup', hue='Tau', 
                 marker='o', palette='tab10', linewidth=2, ax=ax)
    
    ax.axhline(1.0, color='red', linestyle='--', label='Equality')
    
    ax.set_title(f"Speedup of {alt_method} relative to {ref_method}")
    ax.set_xlabel("Gamma")
    ax.set_ylabel(f"Ratio of Times ({ref_method} / {alt_method})")
    
    ax.fill_between(pivot_df['Gamma'].unique(), 1.0, pivot_df['Speedup'].max() * 1.1, 
                    color='green', alpha=0.05, label=f'Victory zone {alt_method}')
    
    ax.legend(title='Tau')

    # Fig02
    ax = axes[0, 1]
    df_iter = df.groupby(['Tau', 'Method'])['Iterations'].mean().reset_index()
    sns.lineplot(data=df_iter, x='Tau', y='Iterations', hue='Method', marker='o', palette={ref_method: '#d62728', alt_method: '#2ca02c'}, ax=ax)
    ax.set_title("Impact of tau on interations")

    # Fig03
    ax = axes[1, 0]
    heat_time = pivot_df.pivot(index="Gamma", columns="Tau", values="Time_Saved")
    sns.heatmap(heat_time, cmap="RdYlGn", center=0, annot=True, fmt=".4f", ax=ax)
    ax.set_title("Time saving matrix (s)")

    # Fig04
    ax = axes[1, 1]
    heat_iter = pivot_df.pivot(index="Gamma", columns="Tau", values="Iter_Saved")
    sns.heatmap(heat_iter, cmap="RdYlGn", center=0, annot=True, fmt=".1f", ax=ax)
    ax.set_title("Iteration gain matrix")

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)

    base_name_no_ext = os.path.splitext(file_name)[0]
    output_image_name = f"{base_name_no_ext}_{ref_method}_{alt_method}.png"
    target_directory = os.path.dirname(file_path)
    output_image_path = os.path.join(target_directory, output_image_name)

    plt.savefig(output_image_path, dpi=300, bbox_inches='tight')
    # plt.show()

method = False  # True, False

if method:  # Window with Tkinter
    import tkinter as tk
    from tkinter import filedialog

    def main():
        root = tk.Tk()
        root.withdraw()
    
        file_path = filedialog.askopenfilename(
            title="Select the results file",
            filetypes=[("Données", "*.csv *.txt"), ("Tous", "*.*")]
        )
        
        if file_path:
            analyze_benders_results(file_path)
        else:
            print("No file selected.")

    if __name__ == "__main__":
        main()
        print("===== Generation complete =====")
else:
    import argparse
    import os
    import matplotlib
    matplotlib.use('Agg')

    def main():
        
        parser = argparse.ArgumentParser(description="Analyse des résultats d'optimisation.")
        parser.add_argument("file_path", help="Chemin du fichier CSV/TXT à analyser", nargs="?")
        
        args = parser.parse_args()
        file_path = args.file_path

        if not file_path:
            file_path = input("Entrez le chemin du fichier de résultats : ").strip()

        if file_path and os.path.exists(file_path):
            analyze_benders_results(file_path)
        else:
            print(f"Erreur : Le fichier '{file_path}' n'existe pas ou est invalide.")

    if __name__ == "__main__":
        main()
        print("===== Generation complete =====")