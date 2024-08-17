import numpy as np
import scipy.stats as stats

A1=27
A2=13
B1=423
B2=3227 
C1=51  
C2=168
D1=1224 
D2=9012

# Step 1: Define your contingency tables
# Immune gene loss table (example values, replace with your data)
immune_gene_table = np.array([[A1, B1], [A2, B2]])

# Background gene loss table (example values, replace with your data)
background_gene_table = np.array([[C1, D1], [C2, D2]])

# Step 2: Perform Fisher's exact test
odds_ratio_immune, p_value_immune = stats.fisher_exact(immune_gene_table)
odds_ratio_background, p_value_background = stats.fisher_exact(background_gene_table)

corrected_odds_ratio = odds_ratio_immune / odds_ratio_background

# Step 4: Output results
print(f"Immune Gene Loss Fisher's Test: Odds Ratio = {odds_ratio_immune}, p-value = {p_value_immune}")
print(f"Background Gene Loss Fisher's Test: Odds Ratio = {odds_ratio_background}, p-value = {p_value_background}")
print(f"Corrected Odds Ratio for Immune Gene Loss: {corrected_odds_ratio}")

n_perm = 10000
perm_corrected_or = []

# Step 5: Combined data for permutation
combined_data = np.array(
	[['Group1', 'Immune', 'Loss']] * A1 + [['Group1', 'Immune', 'Retained']] * B1 +
    [['Group2', 'Immune', 'Loss']] * A2 + [['Group2', 'Immune', 'Retained']] * B2 +
    [['Group1', 'Background', 'Loss']] * C1 + [['Group1', 'Background', 'Retained']] * D1 +
    [['Group2', 'Background', 'Loss']] * C2 + [['Group2', 'Background', 'Retained']] * D2)


for _ in range(n_perm):
    np.random.shuffle(combined_data[:, 0])  # Permute group labels
    immune_perm_table = np.array([
        [np.sum((combined_data[:, 0] == 'Group1') & (combined_data[:, 1] == 'Immune') & (combined_data[:, 2] == 'Loss')),
         np.sum((combined_data[:, 0] == 'Group1') & (combined_data[:, 1] == 'Immune') & (combined_data[:, 2] == 'Retained'))],
        [np.sum((combined_data[:, 0] == 'Group2') & (combined_data[:, 1] == 'Immune') & (combined_data[:, 2] == 'Loss')),
         np.sum((combined_data[:, 0] == 'Group2') & (combined_data[:, 1] == 'Immune') & (combined_data[:, 2] == 'Retained'))]
    ])
    background_perm_table = np.array([
        [np.sum((combined_data[:, 0] == 'Group1') & (combined_data[:, 1] == 'Background') & (combined_data[:, 2] == 'Loss')),
         np.sum((combined_data[:, 0] == 'Group1') & (combined_data[:, 1] == 'Background') & (combined_data[:, 2] == 'Retained'))],
        [np.sum((combined_data[:, 0] == 'Group2') & (combined_data[:, 1] == 'Background') & (combined_data[:, 2] == 'Loss')),
         np.sum((combined_data[:, 0] == 'Group2') & (combined_data[:, 1] == 'Background') & (combined_data[:, 2] == 'Retained'))]
    ])
    or_immune_perm, _ = stats.fisher_exact(immune_perm_table)
    or_background_perm, _ = stats.fisher_exact(background_perm_table)
    perm_corrected_or.append(or_immune_perm / or_background_perm)

# Calculate the empirical p-value
empirical_p_value = np.mean(np.array(perm_corrected_or) >= corrected_odds_ratio)
print(f"Empirical p-value: {empirical_p_value}")