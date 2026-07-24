
from argparse import Namespace
import pandas as pd
from pathlib import Path

def generate_summary(file_path):
    path_in = Path(file_path)
    
    df = pd.read_csv(path_in, sep=',', names=['method', 'gamma', 'tau', 'time', 'iter'])
    df.columns = df.columns.str.strip()

    gamma_values: list[int] = [1, 11, 21, 31, 41, 51, 61, 61, 81, 91, 101]
    df_filtered = df[df['gamma'].isin(gamma_values)]
    
    summary_table: DataFrame = pd.pivot_table(
        df_filtered,
        index='gamma',
        columns='method',
        values=['time', 'iter'],
        aggfunc='mean'
    )
    
    methods: list[str] = sorted(df_filtered['method'].unique())
    columns_filtered: list[tuple[str, str]] = [(metric, method) for method in methods for metric in ['time', 'iter']]
    summary_table = summary_table[columns_filtered]
    summary_table.columns = [f"{metric}_{method}" for metric, method in summary_table.columns]
    summary_table = summary_table.reset_index()
    
    filename_out: str = f"{path_in.stem}_summary.csv"
    path_out: Path = path_in.parent / filename_out
    
    summary_table.to_csv(path_out, index=False)
    print(f"Succes: File exported to\n-> {path_out}")
    
    return summary_table

method = False  # True if you want to use the GUI, False if you want to use the command line

if method:  # Window with Tkinter
    import tkinter as tk
    from tkinter import filedialog

    def main():
        root = tk.Tk()
        root.withdraw()
    
        file_path = filedialog.askopenfilename(
            title="Select the results file",
            filetypes=[("Data", "*.csv"), ("All", "*.*")]
        )
        
        if file_path:
            df_final = generate_summary(file_path)
            print("\nSummary table:")
            print(df_final)
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
        parser.add_argument("file_path", help="Path to the CSV file to analyze", nargs="?")
        
        args: Namespace = parser.parse_args()
        file_path = args.file_path

        if not file_path:
            file_path: str = input("Enter the path to the results file: ").strip()

        if file_path and os.path.exists(file_path):
            df_final = generate_summary(file_path)
            print("\nSummary table:")
            print(df_final)
        else:
            print(f"Error: The file '{file_path}' does not exist or is invalid.")

    if __name__ == "__main__":
        main()