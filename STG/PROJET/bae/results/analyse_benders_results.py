import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tkinter as tk
import os
from tkinter import filedialog

def analyze_benders_results(file_path):
    # Configuration
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
    max_time = max(pivot_df[f'Time_{ref_method}'].max(), pivot_df[f'Time_{alt_method}'].max()) * 1.1
    sns.scatterplot(data=pivot_df, x=f'Time_{ref_method}', y=f'Time_{alt_method}', hue='Gamma', size='Tau', sizes=(20, 200), palette='viridis', ax=ax, alpha=0.8)
    ax.plot([0, max_time], [0, max_time], 'r--', label='Equality')
    ax.fill_between([0, max_time], [0, max_time], [0, 0], color='green', alpha=0.1, label=f'{alt_method} victory zone')
    ax.set_title("Overall resolution time (s)")
    ax.set_xlim(0, max_time)
    ax.set_ylim(0, max_time)
    ax.legend(title='Gamma budget')

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
    plt.show()

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