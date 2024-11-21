import matplotlib.pyplot as plt

# Define the points
points = {
    '4_files': (0.14530159017874744, 0.7338338027501042),
    '7_files': (0.14486317846326538, 0.7346534460067548),
    '16_files': (0.13304170843951246, 0.7253669848138146),
    'clustal_omega': (0.25356285576640764, 0.7221922998445067),
    'mafft': (0.23532388170161972, 0.7522102321758275),
    'muscle5_divvy': (0.0417322150362372, 0.521403981175085),
    'muscle5_partial_filtering': (0.031842732935118745, 0.47128731235615734),
    'muscle5': (0.21938855618981626, 0.7698268362364417),
    'mergealign': (0.21941503174888144, 0.7682606860798054),
    'TrimAl': (0.21938855618981626, 0.775),
    'Ground Truth': (0.0, 1.0),
    'hmmcleaner':()
}



# Define symbols for each point
symbols = ['o', 's', 'D', 'v', 'P', '*', 'X', 'h', '^', '+', 'p']

# Extract x and y coordinates
x_coords = [point[0] for point in points.values()]
y_coords = [point[1] for point in points.values()]

# Plot the points with symbols and labels
for i, (name, point) in enumerate(points.items()):
    if name in ['4_files', '7_files']:
        offset = 0.0025 if name == '4_files' else -0.0025
        plt.plot(point[0] + offset, point[1], marker=symbols[i], markersize=12, label=name)
    elif name == 'Ground Truth':
        plt.plot(point[0], point[1], marker='s', markersize=15, label=name, alpha=0.7)  # Adjusted marker style and size for 'Ground Truth'
        plt.annotate('Ground Truth', xy=(point[0], point[1]), xytext=(point[0]+0.05, point[1]-0.15),
                     arrowprops=dict(facecolor='black', arrowstyle='->'), fontsize=10)
    else:
        plt.plot(point[0], point[1], marker=symbols[i], markersize=12, label=name)

# Set the axis limits
plt.xlim(0, 0.3)
plt.ylim(0.4, 1)

# Set the axis labels
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')

# Create a legend
plt.legend(loc='upper right', fontsize='small')  # Adjusted legend position to 'upper right' and font size to 'small'

# Set the title
plt.title('')

# Display the plot
plt.show()
