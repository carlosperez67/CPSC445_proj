import py7zr
import os

archive_file = 'data.7z'
output_dir = 'pbmc_samples_batch1'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print(f"Extracting {archive_file}...")
with py7zr.SevenZipFile(archive_file, mode='r') as archive:
    archive.extractall(path=output_dir)

print(f"Success! Files extracted to {output_dir}")
print(f"Total files: {len(os.listdir(output_dir))}")
