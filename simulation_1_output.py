# -*- coding: utf-8 -*-
from odbAccess import openOdb
import csv

Z_LENGTH = 10.0  # <-- set your model's Z_LENGTH here

# Open the ODB file
odb = openOdb(path='HeatJob-1.odb')

# Define output CSV file
csv_file = open('pure_data.csv', mode='w')
csv_writer = csv.writer(csv_file)

# Write header row
csv_writer.writerow(['Time (s)', 'Node Label', 'X', 'Y', 'Z', 'Temperature (°C)'])

# Access the instance
assembly = odb.rootAssembly
instance = assembly.instances['COMPOSITE_INSTANCE']

line_count = 0

# Loop over each frame (time step) in the step
step = odb.steps['heat_step']
for frame in step.frames:
    time = frame.frameValue
    temp_field = frame.fieldOutputs['NT11']

    # Loop through temperature values
    for value in temp_field.values:
        node_label = value.nodeLabel
        coords = instance.getNodeFromLabel(node_label).coordinates
        z = coords[2]
        # Only include nodes near z=0
        if 0 <= z <= Z_LENGTH / 8000000:
            temp = value.data
            csv_writer.writerow([time, node_label, coords[0], coords[1], z, temp])
	    line_count += 1

	    if line_count % 20000 == 0:
                
		print "✔️  %d lines written..." % line_count

# Clean up
csv_file.close()
odb.close()
print("✅ Temperature data written to 'pure_data.csv'")
