import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read CSV file
df = pd.read_csv('file_sizes.csv')

# Print column names to debug
print("Column names in CSV:")
print(df.columns.tolist())
print("\nFirst few rows:")
print(df.head())

# Clean column names (remove extra spaces)
df.columns = df.columns.str.strip()

# Create subplots for each error rate
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('File Size Comparison: FastGA PAF vs Tracepoint Encoded PAF', fontsize=16, fontweight='bold')

error_rates = [0.001, 0.01, 0.05, 0.1]
positions = [(0, 0), (0, 1), (1, 0), (1, 1)]

colors = {'Standard': '#1f77b4', 'Mixed': '#ff7f0e', 'Variable': '#2ca02c', 'FastGA': '#d62728'}

for idx, error_rate in enumerate(error_rates):
    ax = axes[positions[idx]]
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    # Get unique sequence lengths and tracepoint types
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    x_positions = np.arange(len(seq_lengths))
    
    # Plot FastGA PAF and Tracepoint Encoded PAF for each type
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            fastga_sizes = []
            encoded_sizes = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    fastga_sizes.append(sl_data['FastGA PAF (MB)'].values[0])
                    encoded_sizes.append(sl_data['Tracepoint Encoded PAF (MB)'].values[0])
                else:
                    fastga_sizes.append(0)
                    encoded_sizes.append(0)
            
            # Plot FastGA PAF (circles, dashed line)
            ax.plot(x_positions, fastga_sizes, 'o--', color=colors[tp_type], 
                    label=f'{tp_type} FastGA PAF', markersize=6, alpha=0.7, linewidth=1)
            
            # Plot Tracepoint Encoded PAF (squares, solid line)
            ax.plot(x_positions, encoded_sizes, 's-', color=colors[tp_type], 
                    label=f'{tp_type} Tracepoint', markersize=6, alpha=0.9, linewidth=2)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('File Size (MB)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    # Set log scale for y-axis if there are large differences
    if error_rate in [0.001, 0.01]:  # For datasets with FastGA's large 10000bp files
        ax.set_yscale('log')

plt.tight_layout()
plt.show()

# Create a comprehensive bar chart comparison
plt.figure(figsize=(20, 10))

# Calculate global min and max for consistent Y scale
all_ratios = []
for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    fastga_size = sl_data['FastGA PAF (MB)'].values[0]
                    encoded_size = sl_data['Tracepoint Encoded PAF (MB)'].values[0]
                    if fastga_size > 0:  # Avoid division by zero
                        ratio = encoded_size / fastga_size
                        all_ratios.append(ratio)

# Set Y scale limits with some padding
y_min = min(all_ratios) * 0.9
y_max = max(all_ratios) * 1.1

# Store handles and labels for single legend
handles = []
labels = []

for idx, error_rate in enumerate(error_rates):
    plt.subplot(2, 2, idx + 1)
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    # Get unique sequence lengths and tracepoint types
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    # Create grouped bar chart
    x = np.arange(len(seq_lengths))
    width = 0.08
    
    for i, tp_type in enumerate(tracepoint_types):
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            fastga_sizes = []
            encoded_sizes = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    fastga_sizes.append(sl_data['FastGA PAF (MB)'].values[0])
                    encoded_sizes.append(sl_data['Tracepoint Encoded PAF (MB)'].values[0])
                else:
                    fastga_sizes.append(0)
                    encoded_sizes.append(0)
            
            # Convert to numpy arrays for division
            fastga_array = np.array(fastga_sizes)
            encoded_array = np.array(encoded_sizes)
            ratios = encoded_array / fastga_array
            
            # Only collect legend info from first subplot
            if idx == 0:
                bar = plt.bar(x + i*width*2 + width, ratios, width, 
                           label=f'{tp_type} Tracepoint', color=colors[tp_type], alpha=1.0, hatch='//')
                handles.append(bar)
                labels.append(f'{tp_type} Tracepoint')
            else:
                plt.bar(x + i*width*2 + width, ratios, width, 
                       color=colors[tp_type], alpha=1.0, hatch='//')
    
    plt.xlabel('Sequence Length (bp)')
    plt.ylabel('PAF with Tracepoints / PAF with Cigar')
    plt.title(f'Error Rate: {error_rate}')
    plt.xticks(x + width*3.5, [f'{sl}' for sl in seq_lengths])
    plt.grid(True, alpha=0.3)
    
    # Set consistent Y scale for all subplots
    plt.ylim(y_min, y_max)
    
    # Add horizontal line at y=1 for reference
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.5, linewidth=1)

# Add single legend outside the subplots
plt.figlegend(handles, labels, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)

plt.tight_layout()
plt.suptitle('Compression Ratios: Tracepoint Encoded PAF vs FastGA PAF', 
             fontsize=16, fontweight='bold', y=0.98)
plt.subplots_adjust(right=0.85)  # Make room for legend
plt.show()

# Create compression ratio analysis
print("\nCompression Ratios (Tracepoint Encoded / FastGA PAF):")
print("=" * 80)
print(f"{'Error':<6} {'Type':<8} {'Length':<7} {'FastGA':<8} {'Encoded':<8} {'Ratio':<6} {'Status'}")
print("-" * 80)

for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    for tp_type in ['Standard', 'Mixed', 'Variable', 'FastGA']:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in sorted(tp_data['Sequence Length'].unique()):
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    fastga_size = sl_data['FastGA PAF (MB)'].values[0]
                    encoded_size = sl_data['Tracepoint Encoded PAF (MB)'].values[0]
                    ratio = encoded_size / fastga_size
                    status = "Compression" if ratio < 1.0 else "Expansion"
                    print(f"{error_rate:<6.3f} {tp_type:<8} {sl:<7} {fastga_size:<8.2f} {encoded_size:<8.2f} {ratio:<6.3f} {status}")