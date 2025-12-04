import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Read CSV file
df = pd.read_csv('csv/decoding.csv')

# Clean column names (remove extra spaces)
df.columns = df.columns.str.strip()

print("Column names in CSV:")
print(df.columns.tolist())

error_rates = [0.001, 0.01, 0.05, 0.1]
positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
colors = {'Standard': '#1f77b4', 'Mixed': '#ff7f0e', 'Variable': '#2ca02c', 'FastGA': '#d62728'}
markers = {'Standard': 'o', 'Mixed': 's', 'Variable': '^', 'FastGA': 'd'}

# 1. Average CPU Time per Alignment (Decoding)
fig1, axes1 = plt.subplots(2, 2, figsize=(16, 12))
fig1.suptitle('Average CPU Time per Alignment - Decoding', fontsize=16, fontweight='bold')

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

y_min_cpu = 0
y_max_cpu = max(all_cpu_times) * 1.1
handles1 = []
labels1 = []

for idx, error_rate in enumerate(error_rates):
    ax = axes1[positions[idx]]
    error_data = df[df['Error'] == error_rate]
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
            
            if idx == 0:
                line = ax.plot(x_positions, cpu_times, marker=markers[tp_type], 
                             label=f'{tp_type}', color=colors[tp_type], 
                             linewidth=2, markersize=8, alpha=0.8)
                handles1.append(line[0])
                labels1.append(f'{tp_type}')
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
fig2, axes2 = plt.subplots(2, 2, figsize=(16, 12))
fig2.suptitle('Average Runtime per Alignment - Decoding', fontsize=16, fontweight='bold')

all_runtimes = []
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
                    runtime = sl_data['Average Runtime (per aignment) (ms)'].values[0]
                    all_runtimes.append(runtime)

y_min_runtime = 0
y_max_runtime = max(all_runtimes) * 1.1
handles2 = []
labels2 = []

for idx, error_rate in enumerate(error_rates):
    ax = axes2[positions[idx]]
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    x_positions = np.arange(len(seq_lengths))
    
    for tp_type in tracepoint_types:
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
            
            if idx == 0:
                line = ax.plot(x_positions, runtimes, marker=markers[tp_type], 
                             label=f'{tp_type}', color=colors[tp_type], 
                             linewidth=2, markersize=8, alpha=0.8)
                handles2.append(line[0])
                labels2.append(f'{tp_type}')
            else:
                ax.plot(x_positions, runtimes, marker=markers[tp_type], 
                       color=colors[tp_type], linewidth=2, markersize=8, alpha=0.8)
    
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
fig3, axes3 = plt.subplots(2, 2, figsize=(16, 12))
fig3.suptitle('Peak Memory Usage - Decoding', fontsize=16, fontweight='bold')

all_memory = []
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
                    memory = sl_data['Peak Memory (MB)'].values[0]
                    all_memory.append(memory)

y_min_memory = 0
y_max_memory = max(all_memory) * 1.1
handles3 = []
labels3 = []

for idx, error_rate in enumerate(error_rates):
    ax = axes3[positions[idx]]
    error_data = df[df['Error'] == error_rate]
    seq_lengths = sorted(error_data['Sequence Length'].unique())
    tracepoint_types = ['Standard', 'Mixed', 'Variable', 'FastGA']
    x_positions = np.arange(len(seq_lengths))
    
    for tp_type in tracepoint_types:
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
            
            if idx == 0:
                line = ax.plot(x_positions, memory_usage, marker=markers[tp_type], 
                             label=f'{tp_type}', color=colors[tp_type], 
                             linewidth=2, markersize=8, alpha=0.8)
                handles3.append(line[0])
                labels3.append(f'{tp_type}')
            else:
                ax.plot(x_positions, memory_usage, marker=markers[tp_type], 
                       color=colors[tp_type], linewidth=2, markersize=8, alpha=0.8)
    
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

# 4. Total CPU Time (Decoding)
fig4, axes4 = plt.subplots(2, 2, figsize=(16, 12))
fig4.suptitle('Total CPU Time - Decoding', fontsize=16, fontweight='bold')

all_total_cpu = []
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
                    total_cpu = sl_data['CPU Time (s)'].values[0]
                    all_total_cpu.append(total_cpu)

y_min_total = 0
y_max_total = max(all_total_cpu) * 1.1
handles4 = []
labels4 = []

for idx, error_rate in enumerate(error_rates):
    ax = axes4[positions[idx]]
    error_data = df[df['Error'] == error_rate]
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
            
            if idx == 0:
                line = ax.plot(x_positions, total_cpu_times, marker=markers[tp_type], 
                             label=f'{tp_type}', color=colors[tp_type], 
                             linewidth=2, markersize=8, alpha=0.8)
                handles4.append(line[0])
                labels4.append(f'{tp_type}')
            else:
                ax.plot(x_positions, total_cpu_times, marker=markers[tp_type], 
                       color=colors[tp_type], linewidth=2, markersize=8, alpha=0.8)
    
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Total CPU Time (s)')
    ax.set_title(f'Error Rate: {error_rate}')
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f'{sl}' for sl in seq_lengths])
    ax.grid(True, alpha=0.3)
    ax.set_ylim(y_min_total, y_max_total)

fig4.legend(handles4, labels4, loc='center right', bbox_to_anchor=(1.0, 0.5), fontsize=10)
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.show()

# Print analysis table
print("\nDecoding Performance Analysis:")
print("=" * 100)
print(f"{'Error':<6} {'Type':<8} {'Length':<7} {'Avg CPU (ms)':<12} {'Avg Runtime (ms)':<16} {'Peak Memory (MB)':<16}")
print("-" * 100)

for error_rate in error_rates:
    error_data = df[df['Error'] == error_rate]
    for tp_type in ['Standard', 'Mixed', 'Variable', 'FastGA']:
        tp_data = error_data[error_data['Tracepoint Type'] == tp_type]
        if not tp_data.empty:
            for sl in sorted(tp_data['Sequence Length'].unique()):
                sl_data = tp_data[tp_data['Sequence Length'] == sl]
                if not sl_data.empty:
                    avg_cpu = sl_data['Average CPU Time (Per Aignment) (ms)'].values[0]
                    avg_runtime = sl_data['Average Runtime (per aignment) (ms)'].values[0]
                    peak_memory = sl_data['Peak Memory (MB)'].values[0]
                    print(f"{error_rate:<6.3f} {tp_type:<8} {sl:<7} {avg_cpu:<12.3f} {avg_runtime:<16.3f} {peak_memory:<16.1f}")