import os
import re

def replace_in_file(file_path):
    with open(file_path, 'r+') as file:
        content = file.read()

        patterns = [
            (r"(?<!self)\.A ", ".toarray() "),
            (r"(?<!self)\.A\n", ".toarray()\n"),
            (r"(?<!self)\.A\)", ".toarray())"),
            (r"(?<!self)\.A\.", "toarray()."),
            (r"(?<!self)\.A1 ", ".toarray().flatten() "),
            (r"(?<!self)\.A1\n", ".toarray().flatten()\n"),
            (r"(?<!self)\.A1\)", ".toarray().flatten())"),
            (r"(?<!self)\.A1\.", ".toarray().flatten()."),
            (r"\\\.", "."),
            ]
        
        for (pattern, replacment) in patterns:
            content = re.sub(pattern, replacment, content)

        file.seek(0)
        file.write(content)
        f.truncate()

def process_directory(directory_path):
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for filename in filenames:
            if filename.endswith(".py"):
                file_path = os.path.join(dirpath, filename)
                replace_in_file(file_path)

# Example usage:
directory_path = '/home/jd438446/miniconda3/envs/scvelo/lib/python3.10/site-packages/scvelo/'
process_directory(directory_path)
