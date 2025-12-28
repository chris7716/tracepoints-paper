import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read CSV file
df = pd.read_csv('csv/decoding.csv')

# Clean column names (remove extra spaces)
df.columns = df.columns.str.strip()

print("Column names in CSV:")
print(df.columns.tolist())
print("\nUnique values in data:")
print("Error rates:", sorted(df['Error'].unique()))
print("Tracepoint types:", sorted(df['Tracepoint Type'].unique()))
print("Sequence lengths:", sorted(df['Sequence Length'].unique()))

# Use actual data values
error_rates = sorted(df['Error'].unique())  # [0.01, 0.05, 0.1, 0.2]
tracepoint_types = sorted(df['Tracepoint Type'].unique())  # ['fastga', 'mixed', 'standard', 'variable']

# Create proper subplot layout
if len(error_rates) == 4:
    positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
    fig_size = (16, 12)
    subplot_shape = (2, 2)
else:
    # Adjust for different number of error rates
    rows = (len(error_rates) + 1) // 2
    positions = [(i // 2, i % 2) for i in range(len(error_rates))]
    fig_size = (16, 6 * rows)
    subplot_shape = (rows, 2)

# Color and marker mapping for actual types
colors = {
    'fastga': '#d62728',     # Red
    'mixed': '#ff7f0e',      # Orange  
    'standard': '#1f77b4',   # Blue
    'variable': '#2ca02c'    # Green
}
markers = {
    'fastga': 'd',           # Diamond
    'mixed': 's',            # Square
    'standard': 'o',         # Circle
    'variable': '^'          # Triangle
}

# Convert Peak Memory from KB to MB
df['Peak Memory (MB)'] = df['Peak Memory (KB)'] / 1024

# 1. Average CPU Time per Alignment (Decoding)
fig1, axes1 = plt.subplots(*subplot_shape, figsize=fig_size)
if len(error_rates) == 1:
    axes1 = [axes1]
elif len(error_rates) <= 2:
    axes1 = axes1.flatten()
fig1.suptitle('Average CPU Time per Alignment - Decoding', fontsize=16, fontweight='bold')

# Calculate global min and max for consistent Y scale
all_cpu_times = []
for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    cpu_time = sl_data['Average CPU Time (Per Aignment) (ms)'].values[0]
                    all_cpu_times.append(cpu_time)

y_min_cpu = 0
y_max_cpu = max(all_cpu_times) * 1.1 if all_cpu_times else 100
handles1 = []
labels1 = []

for idx, error_rate in enumerate(error_rates):
    if len(error_rates) <= 2:
        ax = axes1[idx] if len(error_rates) > 1 else axes1[0]
    else:
        ax = axes1[positions[idx]]
    
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
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
            
            if idx == 0:
                line = ax.plot(x_positions, cpu_times, marker=markers[tp_type], 
                             label=tp_type.capitalize(), color=colors[tp_type], 
                             linewidth=2, markersize=8, alpha=0.8)
                handles1.append(line[0])
                labels1.append(tp_type.capitalize())
            else:
                ax.plot(x_positions, cpu_times, marker=markers[tp_type], 
                       color=colors[tp_type], linewidth=2, markersize=8, alpha=0.8)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Average CPU Time (ms)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    ax.set_ylim(y_min_cpu, y_max_cpu)

fig1.legend(handles1, labels1, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.show()

# 2. Average Runtime per Alignment (Decoding)
fig2, axes2 = plt.subplots(*subplot_shape, figsize=fig_size)
if len(error_rates) == 1:
    axes2 = [axes2]
elif len(error_rates) <= 2:
    axes2 = axes2.flatten()
fig2.suptitle('Average Runtime per Alignment - Decoding', fontsize=16, fontweight='bold')

all_runtimes = []
for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    runtime = sl_data['Average Runtime (per aignment) (ms)'].values[0]
                    all_runtimes.append(runtime)

y_min_runtime = 0
y_max_runtime = max(all_runtimes) * 1.1 if all_runtimes else 10
handles2 = []
labels2 = []

for idx, error_rate in enumerate(error_rates):
    if len(error_rates) <= 2:
        ax = axes2[idx] if len(error_rates) > 1 else axes2[0]
    else:
        ax = axes2[positions[idx]]
    
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    x_positions = np.arange(len(seq_lengths))
    width = 0.8 / len(tracepoint_types)  # Adjust width based on number of types
    
    for i, tp_type in enumerate(tracepoint_types):
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            runtimes = []
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    runtime = sl_data['Average Runtime (per aignment) (ms)'].values[0]
                    runtimes.append(runtime)
                else:
                    runtimes.append(0)
            
            offset = i * width - width * (len(tracepoint_types) - 1) / 2
            
            if idx == 0:
                bar = ax.bar(x_positions + offset, runtimes, width, 
                           label=tp_type.capitalize(), color=colors[tp_type], alpha=0.8)
                handles2.append(bar)
                labels2.append(tp_type.capitalize())
            else:
                ax.bar(x_positions + offset, runtimes, width, 
                       color=colors[tp_type], alpha=0.8)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Average Runtime (ms)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    ax.set_ylim(y_min_runtime, y_max_runtime)

fig2.legend(handles2, labels2, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.show()

# 3. Peak Memory Usage (Decoding)
fig3, axes3 = plt.subplots(*subplot_shape, figsize=fig_size)
if len(error_rates) == 1:
    axes3 = [axes3]
elif len(error_rates) <= 2:
    axes3 = axes3.flatten()
fig3.suptitle('Peak Memory Usage - Decoding', fontsize=16, fontweight='bold')

all_memory = []
for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    memory = sl_data['Peak Memory (MB)'].values[0]
                    all_memory.append(memory)

y_min_memory = 0
y_max_memory = max(all_memory) * 1.1 if all_memory else 50
handles3 = []
labels3 = []

for idx, error_rate in enumerate(error_rates):
    if len(error_rates) <= 2:
        ax = axes3[idx] if len(error_rates) > 1 else axes3[0]
    else:
        ax = axes3[positions[idx]]
    
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    x_positions = np.arange(len(seq_lengths))
    width = 0.8 / len(tracepoint_types)
    
    for i, tp_type in enumerate(tracepoint_types):
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        
        if not tp_data.empty:
            memory_usage = []
            for sl in seq_lengths:
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    memory = sl_data['Peak Memory (MB)'].values[0]
                    memory_usage.append(memory)
                else:
                    memory_usage.append(0)
            
            offset = i * width - width * (len(tracepoint_types) - 1) / 2
            
            if idx == 0:
                bar = ax.bar(x_positions + offset, memory_usage, width, 
                           label=tp_type.capitalize(), color=colors[tp_type], alpha=0.8)
                handles3.append(bar)
                labels3.append(tp_type.capitalize())
            else:
                ax.bar(x_positions + offset, memory_usage, width, 
                       color=colors[tp_type], alpha=0.8)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Peak Memory (MB)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    ax.set_ylim(y_min_memory, y_max_memory)

fig3.legend(handles3, labels3, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.show()

# Print analysis table
print("\nDecoding Performance Analysis:")
print("=" * 100)
print(f"{'Error':<6} {'Type':<10} {'Length':<7} {'Avg CPU (ms)':<12} {'Avg Runtime (ms)':<16} {'Peak Memory (MB)':<16}")
print("-" * 100)

for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    for tp_type in tracepoint_types:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in sorted(tp_data['Sequence Length'].unique()):
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    avg_cpu = sl_data['Average CPU Time (Per Aignment) (ms)'].values[0]
                    avg_runtime = sl_data['Average Runtime (per aignment) (ms)'].values[0]
                    peak_memory = sl_data['Peak Memory (MB)'].values[0]
                    print(f"{error_rate:<6.2f} {tp_type:<10} {sl:<7} {avg_cpu:<12.2f} {avg_runtime:<16.2f} {peak_memory:<16.1f}")