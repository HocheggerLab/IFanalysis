import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ifanalysis._helper_functions import save_fig



def count_per_cond(df: pd.DataFrame) -> pd.DataFrame:
    # Combine groupby, count, and merge operations into a single step
    df_count = (
        df.groupby(["cell_line", "gwli", "palb", "well", "plate_id"])
        ["experiment"]
        .count()
        .reset_index()
        .rename(columns={"experiment": "abs cell count"})
    )

    # Calculate the mean cell count for gwli where palb is 0.0 within the same DataFrame
    mean_cell_count = df_count[df_count['palb'] == 0.0].groupby('gwli')['abs cell count'].mean()

    # Join the mean cell count with the original DataFrame and compute the normalized cell count
    df_merged = df_count.join(mean_cell_count, on='gwli', rsuffix='_mean')
    df_merged['normalized_cell_count'] = (df_merged['abs cell count'] / df_merged['abs cell count_mean']) * 100

    # Drop the auxiliary mean column
    df_merged.drop(columns=['abs cell count_mean'], inplace=True)

    return df_merged
   

def count_plots(df_count, title, count_type, path):
    # Create the grouped bar plot
    fig, ax = plt.subplots(figsize=(12, 8))
    ax = sns.barplot(x="gwli", y=count_type, hue="palb", data=df_count, palette="Blues_d")
    # Add labels and title
    plt.xlabel("GWLI Concentration", fontsize=14)
    plt.ylabel("Normalized Cell Count (%)", fontsize=14)
    plt.title(title, fontsize=16)
    plt.legend(title='PALB Concentration', fontsize='large', title_fontsize='15')

    # Add annotations for better readability
    for p in ax.patches:
        ax.annotate(f"{p.get_height():.2f}", (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha="center", va="bottom", fontsize=10)

    save_fig(fig, path, title)

def intensityplot(df, measurement, title, path):
    # Assuming 'df' is your DataFrame, and you have 'gwli', 'area_cell', and 'palb' columns
    conditions = df['gwli'].unique()
    hue_levels = df['palb'].unique()
    hue_count = len(hue_levels)
    width_per_hue = 0.8 / hue_count  # Default width of violin plots in seaborn

    # Create the plot
    fig, ax = plt.subplots(figsize=(5, 5))

    # Create the violin 
    sns.violinplot(x="gwli", y=measurement, hue="palb", data=df, ax=ax, palette="Blues_d", inner=None, density_norm='width', dodge=True, linewidth=0)

    # Overlay box plots
    for i, condition in enumerate(conditions):
        for j, level in enumerate(hue_levels):
            subset = df[(df["gwli"] == condition) & (df["palb"] == level)]
            box_center = i - 0.4 + width_per_hue/2 + j*width_per_hue
            sns.boxplot(y=subset[measurement], ax=ax,
                        showfliers=False, 
                        #color='white', 
                        saturation=0.5, 
                        width=0.1, 
                        boxprops={'facecolor': (1, 1, 1, 0.5), 'zorder': 2}, # Adjusted alpha here
                        positions=[box_center])

    # Set plot properties
    ax.set_title(title, fontsize=16)
    ax.set_xlabel('GWLI Concentration', fontsize=14)
    ax.set_ylabel('Cell Size (area_cell)', fontsize=14)
   # Set y-axis limits based on 1st and 99th percentiles
    y_min = df[measurement].quantile(0.001)
    y_max = df[measurement].quantile(0.99)
    ax.set_ylim(y_min, y_max)


    # Customizing the legend to only show one legend (from the violin plot)
    handles, labels = ax.get_legend_handles_labels()
    n = len(df['palb'].unique())
    ax.legend(handles[:n], labels[:n], title='PALB Concentration', fontsize='medium', title_fontsize='medium')

    save_fig(fig, path, title)  