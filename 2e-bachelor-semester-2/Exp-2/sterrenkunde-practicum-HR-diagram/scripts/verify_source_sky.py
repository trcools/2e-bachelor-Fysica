#!/usr/bin/env python3
import csv

# Check if Source-Sky values in tabel_rodesterren_1 match export data

print("Verification of Source-Sky values from tabel_rodesterren_1:\n")

# Read tabel_rodesterren
source_sky = {}
with open('data/raw_tables/tabel_rodesterren-1.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        label = row.get('Label', '')
        band = 'B' if 'B-band' in label else 'I' if 'I-Band' in label else None
        if band:
            for i in range(2, 10):  # Check a few stars
                star = f'T{i}'
                val = row.get(f'Source-Sky_{star}', '')
                if val and val != '0':
                    source_sky[(band, star)] = float(val)

print("Source-Sky values from tabel_rodesterren_1.csv:")
for (band, star), val in sorted(source_sky.items()):
    print(f"  {star} ({band}): {val}")

# Read export file
print("\nExport data (flux values):")
rename_map = {f'T{i}': f'T{128+i}' for i in range(2, 10)}

with open('data/exports/belangrijke_data_long.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        star_id = row['star_id']
        if star_id.startswith('T12') or star_id.startswith('T13'):
            orig_num = int(star_id[1:]) - 127
            if 2 <= orig_num <= 9:
                print(f"  {star_id}: flux_B={float(row['flux_B']):.4f}, flux_I={float(row['flux_I']):.4f}")

print("\n✓ If values match Source-Sky from above, data extraction is correct!")
