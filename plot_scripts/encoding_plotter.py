import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read CSV file
df = pd.read_csv('csv/encoding.csv')

# Clean column names (remove extra spaces)
df.columns = df.columns.str.strip()

print("Column names in CSV:")
print(df.columns.tolist())

# Create subplots for each error rate
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Number of Tracepoints by Tracepoint Type', fontsize=16, fontweight='bold')

error_rates = [0.001, 0.01, 0.05, 0.1]
positions = [(0, 0), (0, 1), (1, 0), (1, 1)]

colors = {'Standard': '#1f77b4', 'Mixed': '#ff7f0e', 'Variable': '#2ca02c', 'FastGA': '#d62728'}
markers = {'Standard': 'o', 'Mixed': 's', 'Variable': '^', 'FastGA': 'd'}

# Calculate global min and max for consistent Y scale
all_tracepoints = []
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
                    tracepoints = sl_data['Total Tracepoints'].values[0] / 1000000  # Convert to millions
                    all_tracepoints.append(tracepoints)
                else:
                    all_tracepoints.append(0)

# Set Y scale limits with some padding
y_min = 0
y_max = max(all_tracepoints) * 1.1

# Store handles and labels for single legend
handles = []
labels = []

for idx, error_rate in enumerate(error_rates):
    ax = axes[positions[idx]]
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    # Get unique sequence lengths and tracepoint types
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    x_positions = np.arange(len(seq_lengths))
    width = 0.2
    
    for i, tp_type in enumerate(tracepoint_types):
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            tracepoint_counts = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    tracepoints = sl_data['Total Tracepoints'].values[0] / 1000000  # Convert to millions
                    tracepoint_counts.append(tracepoints)
                else:
                    tracepoint_counts.append(0)
            
            # Only collect legend info from first subplot
            if idx == 0:
                bar = ax.bar(x_positions + i*width, tracepoint_counts, width, 
                           label=f'{tp_type}', color=colors[tp_type], alpha=0.8)
                handles.append(bar)
                labels.append(f'{tp_type}')
            else:
                ax.bar(x_positions + i*width, tracepoint_counts, width, 
                       color=colors[tp_type], alpha=0.8)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Number of Tracepoints (Millions)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions + width*1.5)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    
    # Set consistent Y scale for all subplots
    ax.set_ylim(y_min, y_max)

# Add single legend outside the subplots
fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)

plt.tight_layout()
plt.subplots_adjust(right=0.85)  # Make room for legend
plt.show()

# Create analysis table
print("\nNumber of Tracepoints Analysis:")
print("=" * 80)
print(f"{'Error':<6} {'Type':<8} {'Length':<7} {'Tracepoints':<12} {'Bytes':<12}")
print("-" * 80)

for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    for tp_type in ['Standard', 'Mixed', 'Variable', 'FastGA']:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in sorted(tp_data['Sequence Length'].unique()):
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    tracepoints = sl_data['Total Tracepoints'].values[0]
                    bytes_count = sl_data['Tracepoints in Bytes'].values[0]
                    print(f"{error_rate:<6.3f} {tp_type:<8} {sl:<7} {tracepoints:<12} {bytes_count:<12}")

# Create a second plot showing tracepoints in bytes
plt.figure(figsize=(16, 12))

# Store handles and labels for single legend (second plot)
handles2 = []
labels2 = []

for idx, error_rate in enumerate(error_rates):
    plt.subplot(2, 2, idx + 1)
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    # Get unique sequence lengths and tracepoint types
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    x_positions = np.arange(len(seq_lengths))
    width = 0.2
    
    for i, tp_type in enumerate(tracepoint_types):
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            bytes_counts = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    bytes_count = sl_data['Tracepoints in Bytes'].values[0] / 1024 / 1024  # Convert to MB
                    bytes_counts.append(bytes_count)
                else:
                    bytes_counts.append(0)
            
            # Only collect legend info from first subplot
            if idx == 0:
                bar = plt.bar(x_positions + i*width, bytes_counts, width, 
                           label=f'{tp_type}', color=colors[tp_type], alpha=0.8)
                handles2.append(bar)
                labels2.append(f'{tp_type}')
            else:
                plt.bar(x_positions + i*width, bytes_counts, width, 
                       color=colors[tp_type], alpha=0.8)
    
    plt.xlabel('Sequence Length (bp)')
    plt.ylabel('Tracepoints Size (MB)')
    plt.title(f'Error Rate: {error_rate}')
    plt.xticks(x_positions + width*1.5, [f'{sl}' for sl in seq_lengths])
    plt.grid(True, alpha=0.3)

# Add single legend outside the subplots for second plot
plt.figlegend(handles2, labels2, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)

plt.tight_layout()
plt.suptitle('Tracepoints Storage Size by Tracepoint Type', fontsize=16, fontweight='bold', y=0.98)
plt.subplots_adjust(right=0.85)  # Make room for legend
plt.show()

# Create a third plot showing average CPU time per alignment
fig3, axes3 = plt.subplots(2, 2, figsize=(16, 12))
fig3.suptitle('Average CPU Time per Alignment by Tracepoint Type', fontsize=16, fontweight='bold')

# Calculate global min and max for consistent Y scale
all_cpu_times = []
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
                    cpu_time = sl_data['Average CPU Time (Per Aignment) (ms)'].values[0]
                    all_cpu_times.append(cpu_time)

# Set Y scale limits with some padding
y_min_cpu = 0
y_max_cpu = max(all_cpu_times) * 1.1

# Store handles and labels for single legend
handles3 = []
labels3 = []

for idx, error_rate in enumerate(error_rates):
    ax = axes3[positions[idx]]
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    # Get unique sequence lengths and tracepoint types
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    x_positions = np.arange(len(seq_lengths))
    
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            cpu_times = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    cpu_time = sl_data['Average CPU Time (Per Aignment) (ms)'].values[0]
                    cpu_times.append(cpu_time)
                else:
                    cpu_times.append(0)
            
            # Only collect legend info from first subplot
            if idx == 0:
                line = ax.plot(x_positions, cpu_times, marker=markers[tp_type], 
                             label=f'{tp_type}', color=colors[tp_type], 
                             linewidth=2, markersize=8, alpha=0.8)
                handles3.append(line[0])
                labels3.append(f'{tp_type}')
            else:
                ax.plot(x_positions, cpu_times, marker=markers[tp_type], 
                       color=colors[tp_type], linewidth=2, markersize=8, alpha=0.8)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Average CPU Time (ms)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    
    # Set consistent Y scale for all subplots
    ax.set_ylim(y_min_cpu, y_max_cpu)

# Add single legend outside the subplots
fig3.legend(handles3, labels3, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)

plt.tight_layout()
plt.subplots_adjust(right=0.85)  # Make room for legend
plt.show()

# Create analysis table
print("\nAverage CPU Time per Alignment Analysis:")
print("=" * 80)
print(f"{'Error':<6} {'Type':<8} {'Length':<7} {'CPU Time (ms)':<12} {'Total CPU (s)':<12}")
print("-" * 80)

for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    for tp_type in ['Standard', 'Mixed', 'Variable', 'FastGA']:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in sorted(tp_data['Sequence Length'].unique()):
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    avg_cpu = sl_data['Average CPU Time (Per Aignment) (ms)'].values[0]
                    total_cpu = sl_data['CPU Time (s)'].values[0]
                    print(f"{error_rate:<6.3f} {tp_type:<8} {sl:<7} {avg_cpu:<12.3f} {total_cpu:<12.2f}")

# Create a forth plot showing total CPU time
plt.figure(figsize=(16, 12))

# Store handles and labels for single legend (second plot)
handles2 = []
labels2 = []

for idx, error_rate in enumerate(error_rates):
    plt.subplot(2, 2, idx + 1)
    
    # Filter data for current error rate
    error_data = df[df['Error'] == error_rate]
    
    # Get unique sequence lengths and tracepoint types
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    
    x_positions = np.arange(len(seq_lengths))
    
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            total_cpu_times = []
            
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    total_cpu = sl_data['CPU Time (s)'].values[0]
                    total_cpu_times.append(total_cpu)
                else:
                    total_cpu_times.append(0)
            
            # Only collect legend info from first subplot
            if idx == 0:
                line = plt.plot(x_positions, total_cpu_times, marker=markers[tp_type], 
                              label=f'{tp_type}', color=colors[tp_type], 
                              linewidth=2, markersize=8, alpha=0.8)
                handles2.append(line[0])
                labels2.append(f'{tp_type}')
            else:
                plt.plot(x_positions, total_cpu_times, marker=markers[tp_type], 
                        color=colors[tp_type], linewidth=2, markersize=8, alpha=0.8)
    
    plt.xlabel('Sequence Length (bp)')
    plt.ylabel('Total CPU Time (s)')
    plt.title(f'Error Rate: {error_rate}')
    plt.xticks(x_positions, [f'{sl}' for sl in seq_lengths])
    plt.grid(True, alpha=0.3)

# Add single legend outside the subplots for second plot
plt.figlegend(handles2, labels2, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)

plt.tight_layout()
plt.suptitle('Total CPU Time by Tracepoint Type', fontsize=16, fontweight='bold', y=0.98)
plt.subplots_adjust(right=0.85)  # Make room for legend
plt.show()