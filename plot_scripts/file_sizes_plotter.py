import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read CSV file
df = pd.read_csv('csv/file_sizes.csv')

# Print column names to debug
print("Column names in CSV:")
print(df.columns.tolist())
print("\nFirst few rows:")
print(df.head())
print("\nUnique values in key columns:")
print("Error values:", df['Error'].unique() if 'Error' in df.columns else "Error column not found")
print("Tracepoint Type values:", df['Tracepoint Type'].unique() if 'Tracepoint Type' in df.columns else "Tracepoint Type column not found")
print("Sequence Length values:", df['Sequence Length'].unique() if 'Sequence Length' in df.columns else "Sequence Length column not found")

# Clean column names (remove extra spaces)
df.columns = df.columns.str.strip()

# Check if required columns exist
required_columns = ['Error', 'Tracepoint Type', 'Sequence Length']
missing_columns = [col for col in required_columns if col not in df.columns]

if missing_columns:
    print(f"\nMissing columns: {missing_columns}")
    print("Available columns:", df.columns.tolist())
    exit(1)

# Check for file size columns - try different possible names
file_size_columns = [col for col in df.columns if 'paf' in col.lower() or 'mb' in col.lower()]
print(f"\nFile size related columns: {file_size_columns}")

# Try to identify the correct columns for FastGA PAF and Tracepoint PAF
fastga_col = None
tracepoint_col = None

for col in df.columns:
    col_lower = col.lower()
    if 'fastga' in col_lower and ('paf' in col_lower or 'mb' in col_lower):
        fastga_col = col
    elif 'tp' in col_lower or ('tracepoint' in col_lower and 'paf' in col_lower):
        tracepoint_col = col

if not fastga_col or not tracepoint_col:
    print(f"\nCould not find appropriate columns:")
    print(f"FastGA column: {fastga_col}")
    print(f"Tracepoint column: {tracepoint_col}")
    print("Please check your CSV file structure")
    exit(1)

print(f"\nUsing columns:")
print(f"FastGA PAF: {fastga_col}")
print(f"Tracepoint PAF: {tracepoint_col}")

# Create subplots for each error rate
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('File Size Comparison: FastGA PAF vs Tracepoint Encoded PAF', fontsize=16, fontweight='bold')

# Get actual error rates from data
error_rates = sorted(df['Error'].unique())
print(f"\nActual error rates in data: {error_rates}")

# Ensure we have exactly 4 error rates for the 2x2 grid
if len(error_rates) != 4:
    print(f"Warning: Expected 4 error rates, found {len(error_rates)}")
    # Pad or truncate as needed
    while len(error_rates) < 4:
        error_rates.append(error_rates[-1])  # Duplicate last value
    error_rates = error_rates[:4]  # Take only first 4

positions = [(0, 0), (0, 1), (1, 0), (1, 1)]

# Get actual tracepoint types from data
tracepoint_types = sorted(df['Tracepoint Type'].unique())
print(f"Actual tracepoint types: {tracepoint_types}")

# Create color map for actual types
colors = {}
color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
for i, tp_type in enumerate(tracepoint_types):
    colors[tp_type] = color_list[i % len(color_list)]

all_ratios = []  # For calculating global min/max

for idx, error_rate in enumerate(error_rates):
    ax = axes[positions[idx]]
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    if error_data.empty:
        ax.text(0.5, 0.5, f'No data for error rate {error_rate}', 
                transform=ax.transAxes, ha='center', va='center')
        ax.set_title(f'Error Rate: {error_rate}')
        continue
    
    # Get unique sequence lengths
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    
    if not seq_lengths:
        ax.text(0.5, 0.5, f'No sequence lengths for error rate {error_rate}', 
                transform=ax.transAxes, ha='center', va='center')
        ax.set_title(f'Error Rate: {error_rate}')
        continue
    
    x_positions = np.arange(len(seq_lengths))
    
    # Plot FastGA PAF and Tracepoint Encoded PAF for each type
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            fastga_sizes = []
            encoded_sizes = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty and len(sl_data) > 0:
                    fastga_size = sl_data[fastga_col].iloc[0]
                    encoded_size = sl_data[tracepoint_col].iloc[0]
                    fastga_sizes.append(fastga_size)
                    encoded_sizes.append(encoded_size)
                    
                    # Store ratio for global scaling
                    if fastga_size > 0:
                        all_ratios.append(encoded_size / fastga_size)
                else:
                    fastga_sizes.append(0)
                    encoded_sizes.append(0)
            
            # Only plot if we have data
            if any(size > 0 for size in fastga_sizes + encoded_sizes):
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
    
    # Only add legend if there are labeled artists
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    # Set log scale for y-axis if there are large differences
    if len(seq_lengths) > 0 and any(error_data[fastga_col] > 10):
        ax.set_yscale('log')

plt.tight_layout()
plt.show()

# Create compression ratio analysis only if we have data
if all_ratios:
    # Set Y scale limits with some padding
    y_min = min(all_ratios) * 0.9
    y_max = max(all_ratios) * 1.1

    # Create a comprehensive bar chart comparison
    plt.figure(figsize=(20, 10))
    
    # Store handles and labels for single legend
    handles = []
    labels = []

    for idx, error_rate in enumerate(error_rates):
        plt.subplot(2, 2, idx + 1)
        
        # Filter data for current error rate
        error_data = df[df['Error'] == error_rate]
        
        if error_data.empty:
            continue
        
        # Get unique sequence lengths
        seq_lengths = sorted(error_data['Sequence Length'].unique())
        
        if not seq_lengths:
            continue
        
        # Create grouped bar chart
        x = np.arange(len(seq_lengths))
        width = 0.8 / len(tracepoint_types)  # Adjust width based on number of types
        
        for i, tp_type in enumerate(tracepoint_types):
            tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
            
            if not tp_data.empty:
                ratios = []
                
                for sl in seq_lengths:
                    sl_data = tp_data[tp_data['Sequence Length'] == sl]
                    if not sl_data.empty and len(sl_data) > 0:
                        fastga_size = sl_data[fastga_col].iloc[0]
                        encoded_size = sl_data[tracepoint_col].iloc[0]
                        if fastga_size > 0:
                            ratios.append(encoded_size / fastga_size)
                        else:
                            ratios.append(0)
                    else:
                        ratios.append(0)
                
                # Only plot if we have data
                if any(ratio > 0 for ratio in ratios):
                    # Only collect legend info from first subplot
                    if idx == 0:
                        bar = plt.bar(x + i*width - width*(len(tracepoint_types)-1)/2, ratios, width, 
                                   label=f'{tp_type} Tracepoint', color=colors[tp_type], alpha=0.8)
                        handles.append(bar)
                        labels.append(f'{tp_type} Tracepoint')
                    else:
                        plt.bar(x + i*width - width*(len(tracepoint_types)-1)/2, ratios, width, 
                               color=colors[tp_type], alpha=0.8)
        
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Tracepoint PAF / FastGA PAF')
        plt.title(f'Error Rate: {error_rate}')
        plt.xticks(x, [f'{sl}' for sl in seq_lengths])
        plt.grid(True, alpha=0.3)
        
        # Set consistent Y scale for all subplots
        plt.ylim(y_min, y_max)
        
        # Add horizontal line at y=1 for reference
        plt.axhline(y=1, color='black', linestyle='--', alpha=0.5, linewidth=1)

    # Add single legend outside the subplots
    if handles:
        plt.figlegend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.95), fontsize=14, ncol=len(tracepoint_types))

    plt.tight_layout()
    plt.suptitle('Compression Ratios: Tracepoint Encoded PAF vs FastGA PAF', 
                 fontsize=16, fontweight='bold', y=0.98)
    plt.subplots_adjust(top=0.85)  # Make room for legend
    plt.show()

    # Print compression ratio analysis
    print("\nCompression Ratios (Tracepoint Encoded / FastGA PAF):")
    print("=" * 90)
    print(f"{'Error':<8} {'Type':<12} {'Length':<8} {'FastGA':<10} {'Encoded':<10} {'Ratio':<8} {'Status'}")
    print("-" * 90)

    for error_rate in error_rates:
        error_data = df[df['Error'] == error_rate]
        for tp_type in tracepoint_types:
            tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
            if not tp_data.empty:
                for sl in sorted(tp_data['Sequence Length'].unique()):
                    sl_data = tp_data[tp_data['Sequence Length'] == sl]
                    if not sl_data.empty and len(sl_data) > 0:
                        fastga_size = sl_data[fastga_col].iloc[0]
                        encoded_size = sl_data[tracepoint_col].iloc[0]
                        if fastga_size > 0:
                            ratio = encoded_size / fastga_size
                            status = "Compression" if ratio < 1.0 else "Expansion"
                            print(f"{error_rate:<8.3f} {tp_type:<12} {sl:<8} {fastga_size:<10.2f} {encoded_size:<10.2f} {ratio:<8.3f} {status}")
else:
    print("\nNo valid data found for compression ratio analysis")
    print("Please check your CSV file structure and data")