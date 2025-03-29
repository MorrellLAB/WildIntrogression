#!/usr/bin/env python3

"""
Written with Claude 3.7 Sonnnet through Copilot
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read the interval data
print("Reading interval genetic distances file...")
df = pd.read_csv('interval_genetic_distances.tsv', sep='\t')

# Create a copy of the dataframe for processing
df_orig = df.copy()

# Function to extract numeric ID from sample names (e.g., WBDC123 -> 123)

def extract_numeric_id(sample_name):
    try:
        return int(sample_name.replace('WBDC', ''))
    except (ValueError, AttributeError):
        return 0  # Default value if extraction fails

# Convert 'NA' values to NaN for proper handling
df['Genetic_Distance(cM)'] = pd.to_numeric(df['Genetic_Distance(cM)'], errors='coerce')

# Get counts of intervals per sample
sample_counts = df.groupby('Sample').size().reset_index(name='Count')

# Sort samples by their numeric ID
sample_counts['numeric_id'] = sample_counts['Sample'].apply(extract_numeric_id)
sample_counts = sample_counts.sort_values('numeric_id')
sorted_samples = sample_counts['Sample'].tolist()

# Use the same sample order for both plots (sorted by sample ID)
genetic_order = sorted_samples
physical_order = sorted_samples

# Print basic statistics
print(f"Total intervals: {len(df)}")
print(f"Intervals with genetic distance: {df['Genetic_Distance(cM)'].notna().sum()}")
print(f"Unique samples: {df['Sample'].nunique()}")
print("\nNumber of intervals per sample:")
for _, row in sample_counts.sort_values('Count', ascending=False).iterrows():
    print(f"{row['Sample']}: {row['Count']} intervals")
# Set the aesthetic style of the plots
sns.set_theme(style="whitegrid")

# Reshape the data into a long format for paired plotting
# First, let's create a dataframe with just the columns we need
plot_df = df[['Sample', 'Genetic_Distance(cM)', 'Physical_Distance(Mbp)']].copy()

# Reshape to long format
plot_df = pd.melt(plot_df,
                  id_vars=['Sample'],
                  value_vars=['Genetic_Distance(cM)', 'Physical_Distance(Mbp)'],
                  var_name='Distance_Type', 
                  value_name='Distance')

# Clean up the distance type names for the legend
plot_df['Distance_Type'] = plot_df['Distance_Type'].replace({
    'Genetic_Distance(cM)': 'Genetic (cM)',
    'Physical_Distance(Mbp)': 'Physical (Mbp)'
})

# Set figure height to accommodate all samples (about 0.3 inch per sample)
fig_height = max(12, df['Sample'].nunique() * 0.3)
fig_width = 14

# Create figure
plt.figure(figsize=(fig_width, fig_height))

# Define colors for the two distance types
colors = {"Genetic (cM)": "red", "Physical (Mbp)": "blue"}

# Create paired violin plots
ax = sns.violinplot(y='Sample', x='Distance', 
                   hue='Distance_Type', 
                   data=plot_df,
                   order=sorted_samples,
                   inner='box',
                   orient='h',
                   palette=colors,
                   split=False,
                   cut=0)  # Don't extend the violin past the observed data

# Add swarm plots
sns.swarmplot(y='Sample', x='Distance', 
              hue='Distance_Type',
              data=plot_df, 
              order=sorted_samples,
              size=2, 
              alpha=0.6, 
              orient='h',
              palette=colors,
              dodge=True)  # Dodge to align with violin plots

# Set x-axis to start at 0
plt.xlim(0, None)

# Add titles and labels
plt.title('Interval Lengths by Sample', fontsize=18)
plt.xlabel('Distance', fontsize=14)
plt.ylabel('Sample', fontsize=14)

# Add interval counts to y-axis labels
sample_count_dict = dict(zip(sample_counts['Sample'], sample_counts['Count']))
new_labels = [f"{sample} (n={sample_count_dict.get(sample, 0)})" for sample in sorted_samples]
plt.yticks(range(len(sorted_samples)), new_labels)

# Customize the legend to make it more clear
handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], title="Distance Type", loc="upper right")

# Adjust layout
plt.tight_layout()

# Save the figure
plt.savefig('interval_lengths_by_sample.png', dpi=300)
print("Saved plot to interval_lengths_by_sample.png")
print("Plotting complete!")

