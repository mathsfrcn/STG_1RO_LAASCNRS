import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def analyze_benders_results(file_path):
    try:
        df = pd.read_csv(file_path, sep=r'[,\s]+', engine='python', 
                         names=['Method', 'Gamma', 'Tau', 'nb_path_to_take', 'Iterations', 'Time', 'Time_Master', 'Time_Subproblem'])
    except Exception as e:
        print("Error: reading the file ", e)
        return

    methods = df['Method'].unique()
    if len(methods) != 2:
        print(f"Error: The file must contain exactly 2 methods. Found: {methods}")
        return
    
    methods_sorted = sorted(methods, key=lambda x: 0 if x == 'BA' else (1 if x == 'KC' else 2))
    ref_method = methods_sorted[0]
    alt_method = methods_sorted[1]
    
    df['Instance_ID'] = df.groupby(['Method', 'Gamma', 'Tau']).cumcount()   # Add id for variance

    pivot_df = df.pivot_table(index=['Instance_ID', 'Gamma', 'Tau'], 
                              columns='Method', 
                              values=['Iterations', 'Time']).reset_index()
    
    pivot_df.columns = [f"{col[0]}_{col[1]}" if col[1] else col[0] for col in pivot_df.columns]

    pivot_df['Time_Saved'] = pivot_df[f'Time_{ref_method}'] - pivot_df[f'Time_{alt_method}']
    pivot_df['Iter_Saved'] = pivot_df[f'Iterations_{ref_method}'] - pivot_df[f'Iterations_{alt_method}']
    
    pivot_df['Speedup'] = pivot_df[f'Time_{ref_method}'] / pivot_df[f'Time_{alt_method}']

    file_name = os.path.basename(file_path)
    base_name_no_ext = os.path.splitext(file_name)[0]
    target_directory = os.path.dirname(file_path)

    sns.set_theme(style='whitegrid')

    # Fig.01
    fig1 = plt.figure(figsize=(16, 12))
    fig1.canvas.manager.set_window_title(f"Benders decomposition - {file_name}")
    axes = fig1.subplots(2, 2)
    fig1.suptitle(f'Comparative Analysis: {ref_method} VS {alt_method} ({file_name})', fontsize=16, fontweight='bold')

    # Scatter Plot
    ax = axes[0, 0]

    pivot_df['Speedup'] = pivot_df[f'Time_{ref_method}'] / pivot_df[f'Time_{alt_method}']
    
    sns.lineplot(data=pivot_df, x='Gamma', y='Speedup', hue='Tau', 
                 marker='o', palette='tab10', linewidth=2, ax=ax)
    
    ax.axhline(1.0, color='red', linestyle='--', label='Equality')
    ax.set_title(f"Speedup of {alt_method} relative to {ref_method}")
    ax.set_xlabel("Gamma")
    ax.set_ylabel(f"Ratio of Times ({ref_method} / {alt_method})")
    ax.fill_between(pivot_df['Gamma'].unique(), 1.0, pivot_df['Speedup'].max() * 1.1, color='green', alpha=0.05, label=f'Victory zone {alt_method}')
    ax.legend(title='Tau')


    ax = axes[0, 1]
    sns.lineplot(data=df, x='Gamma', y='Iterations', hue='Method', marker='o', palette={ref_method: '#d62728', alt_method: '#2ca02c'}, ax=ax)
    ax.set_title("Impact of Gamma on averaged iterations")

    # Heatmaps
    pivot_means = pivot_df.groupby(['Gamma', 'Tau']).mean().reset_index()
    
    ax = axes[1, 0]
    heat_time = pivot_means.pivot(index="Gamma", columns="Tau", values="Time_Saved")
    sns.heatmap(heat_time, cmap="RdYlGn", center=0, annot=True, fmt=".4f", ax=ax)
    ax.set_title(f"Average Time Saved (s)")

    ax = axes[1, 1]
    heat_iter = pivot_means.pivot(index="Gamma", columns="Tau", values="Iter_Saved")
    sns.heatmap(heat_iter, cmap="RdYlGn", center=0, annot=True, fmt=".1f", ax=ax)
    ax.set_title(f"Average Iteration Gain")

    fig1.tight_layout()
    fig1.subplots_adjust(top=0.90)
    output_image_1 = os.path.join(target_directory, f"{base_name_no_ext}_{ref_method}_vs_{alt_method}_dashboard.png")
    fig1.savefig(output_image_1, dpi=300, bbox_inches='tight')

    # Fig.02
    fig2 = plt.figure(figsize=(12, 8))
    fig2.canvas.manager.set_window_title(f"Dispersion Analysis - {file_name}")
    ax2 = fig2.add_subplot(111)

    sns.lineplot(data=pivot_df, x='Gamma', y='Speedup', hue='Tau', 
                 estimator='median', errorbar=('pi', 50), 
                 marker='o', palette='tab10', linewidth=2, ax=ax2)

    ax2.axhline(1.0, color='red', linestyle='--', label='Egalite (Speedup = 1)')
    ax2.set_title(f"Speedup Dispersion: {ref_method} vs {alt_method}\n(Line = Media | Shaded area = Interquartile range Q1–Q3)", fontsize=14)
    ax2.set_xlabel("Uncertainty budget", fontsize=12)
    ax2.set_ylabel(f"Ratio Speedup (T_ref / T_alt)", fontsize=12)
    
    ax2.legend(title='Tolerance', bbox_to_anchor=(1.05, 1), loc='upper left')

    fig2.tight_layout()
    output_image_2 = os.path.join(target_directory, f"{base_name_no_ext}_{ref_method}_vs_{alt_method}_dispersion.png")
    fig2.savefig(output_image_2, dpi=300, bbox_inches='tight')

    # Fig.03
    fig3 = plt.figure(figsize=(16, 8))
    fig3.canvas.manager.set_window_title(f"Boxplot Analysis - {file_name}")
    ax3 = fig3.add_subplot(111)

    sns.boxplot(data=pivot_df, x='Gamma', y='Speedup', hue='Tau', 
                palette='tab10', width=0.7, ax=ax3, fliersize=4)

    ax3.axhline(1.0, color='red', linestyle='--', linewidth=2, label='Egalite (Speedup = 1)')
    ax3.set_title(f"Seepup Distribution: {ref_method} vs {alt_method}", fontsize=14, fontweight='bold')
    ax3.set_xlabel("Uncertainty budget", fontsize=12)
    ax3.set_ylabel(f"Ratio Speedup (T_{ref_method} / T_{alt_method})", fontsize=12)
    
    ax3.legend(title='Tolerance', bbox_to_anchor=(1.01, 1), loc='upper left')

    fig3.tight_layout()
    output_image_3 = os.path.join(target_directory, f"{base_name_no_ext}_{ref_method}_vs_{alt_method}_boxplots.png")
    fig3.savefig(output_image_3, dpi=300, bbox_inches='tight')

    print(f"Dashboards saved under:\n- {output_image_1}\n- {output_image_2}\n- {output_image_3}")
    print("Succes: Generation complete")

def main():
    root = tk.Tk()
    root.withdraw()
    
    file_path = filedialog.askopenfilename(
        title="Select the results file",
        filetypes=[("Data", "*.csv *.txt"), ("All Files", "*.*")]
    )
    
    if file_path:
        analyze_benders_results(file_path)
    else:
        print("Error: No file selected")

method = False  # True if you want to use the GUI, False if you want to use the command line

if method:  # Window with Tkinter
    import tkinter as tk
    from tkinter import filedialog

    def main():
        root = tk.Tk()
        root.withdraw()
    
        file_path = filedialog.askopenfilename(
            title="Select the results file",
            filetypes=[("Data", "*.csv *.txt"), ("All Files", "*.*")]
        )
        
        if file_path:
            analyze_benders_results(file_path)
        else:
            print("Error: No file selected")

    if __name__ == "__main__":
        main()
else:
    import argparse
    import os
    import matplotlib
    matplotlib.use('Agg')

    def main():
        
        parser = argparse.ArgumentParser(description="Analyze optimization results.")
        parser.add_argument("file_path", help="Path to the CSV/TXT file to analyze", nargs="?")
        
        args = parser.parse_args()
        file_path = args.file_path

        if not file_path:
            file_path = input("Enter the path to the results file: ").strip()

        if file_path and os.path.exists(file_path):
            analyze_benders_results(file_path)
        else:
            print(f"Error: The file '{file_path}' does not exist or is invalid")

    if __name__ == "__main__":
        main()