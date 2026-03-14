import os

source_dir = '/Users/brandonlukas/Library/CloudStorage/Box-Box/data/timkoh/macrophages'
target_dir = './box'

os.makedirs(target_dir, exist_ok=True)

for name in os.listdir(source_dir):
    src_path = os.path.join(source_dir, name)
    dst_path = os.path.join(target_dir, name)
    if os.path.isdir(src_path):
        if not os.path.exists(dst_path):
            os.symlink(src_path, dst_path)
            print(f"Symlinked: {dst_path} -> {src_path}")
        else:
            print(f"Skipped (already exists): {dst_path}")