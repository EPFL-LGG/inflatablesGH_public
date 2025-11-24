import os
import nbformat as nbf

# Walk through the current directory
for dirpath, dirnames, filenames in os.walk('.'):
    for filename in filenames:
        if filename.endswith('.ipynb'):
            # Construct the full file path
            filepath = os.path.join(dirpath, filename)
            
            # Load the notebook
            with open(filepath, 'r') as f:
                nb = nbf.read(f, as_version=4)

            # Clear the outputs
            for cell in nb.cells:
                if cell['cell_type'] == 'code':
                    cell['outputs'].clear()

            # Save the notebook
            with open(filepath, 'w') as f:
                nbf.write(nb, f)