import pandas as pd
import matplotlib.pyplot as plt

wd = "../../results/04_varcalls"
df = pd.read_csv(f"{wd}/snpdev.genic.m_z.txt", sep='\t')

#Histogram of p-values
plt.figure(figsize=(8,4))
plt.hist(df['pval'], bins=50, color='skyblue', edgecolor='black')
plt.title("Distribution of chi-square p-values")
plt.xlabel("p-value")
plt.ylabel("Count")
plt.tight_layout()
plt.savefig("pval_histogram.png")   # Save file
plt.show()

#Boxplot of p-values by coverage bins
df['coverage_bin'] = (df['COV'] // 500) * 500
bins = sorted(df['coverage_bin'].unique())
data_to_plot = [df.loc[df['coverage_bin'] == b, 'pval'].values for b in bins]

plt.figure(figsize=(12,6))
plt.boxplot(data_to_plot, positions=bins, widths=400, patch_artist=True,
            boxprops=dict(facecolor='lightgreen', color='green'),
            medianprops=dict(color='red'))
plt.xlabel("Coverage bin (total reads)")
plt.ylabel("p-value")
plt.title("Boxplot of p-values by coverage bins")
plt.xticks(bins, rotation=45)
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig("pval_boxplot_by_coverage.png")
plt.show()
