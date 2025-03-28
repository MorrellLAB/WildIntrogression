import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 12))
ax.set_aspect('equal')

# Define dataset sizes
bopa_total = 2963
gbs_total = 619179
iselect_total = 7622

# Use same radius for all circles
radius = 0.25

# Create three equal-sized circles
circle1 = plt.Circle((0.4, 0.6), radius, alpha=0.4, fc='#ff9999', ec='black', label='BOPA')
circle2 = plt.Circle((0.6, 0.6), radius, alpha=0.4, fc='#99ff99', ec='black', label='GBS')
circle3 = plt.Circle((0.5, 0.4), radius, alpha=0.4, fc='#99ccff', ec='black', label='iSelect_9k')

# Add circles to plot
ax.add_patch(circle1)
ax.add_patch(circle2)
ax.add_patch(circle3)

# Add text annotations for intersection values
plt.text(0.5, 0.65, '78', ha='center', va='center', fontsize=12, fontweight='bold')  # BOPA-GBS
plt.text(0.35, 0.45, '2745', ha='center', va='center', fontsize=12, fontweight='bold')  # BOPA-9k
plt.text(0.65, 0.45, '196', ha='center', va='center', fontsize=12, fontweight='bold')  # GBS-9k
plt.text(0.5, 0.5, '69', ha='center', va='center', fontsize=12, fontweight='bold', color='darkred')  # Three-way intersection

# Add dataset labels with total variants
plt.text(0.25, 0.7, f'BOPA\n({bopa_total:,})', ha='center', va='center', fontsize=14, fontweight='bold')
plt.text(0.75, 0.7, f'GBS\n({gbs_total:,})', ha='center', va='center', fontsize=14, fontweight='bold')
plt.text(0.5, 0.2, f'iSelect_9k\n({iselect_total:,})', ha='center', va='center', fontsize=14, fontweight='bold')

# Set title
plt.title('Shared Variants Between Datasets', pad=20, fontsize=16)

# Set plot limits
plt.xlim(0, 1)
plt.ylim(0, 1)

# Remove axes
ax.set_xticks([])
ax.set_yticks([])

# Save the plot
plt.savefig('equal_sized_venn_with_three_way.png', dpi=300, bbox_inches='tight')
plt.show()