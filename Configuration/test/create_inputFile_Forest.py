import os

# Directory path
directory = '/eos/cms/store/group/phys_heavyions/clemahie/RECO/pythiaRECO/MinBias_PbPb_5p36TeV_Hydjet_RECO_Run2/240805_111038/0000'

# Output file path
output_file = 'inputFile_Forest.txt'

# Prefix
prefix = '/store/group/phys_heavyions/clemahie/RECO/pythiaRECO/MinBias_PbPb_5p36TeV_Hydjet_RECO_Run2/240805_111038/0000/'

# Open the output file in write mode
with open(output_file, 'w') as f:
    # Loop through each file in the directory
    for filename in os.listdir(directory):
        # Check if the 'filename' is a file (not a directory)
        if os.path.isfile(os.path.join(directory, filename)):
            # Write the prefixed filename to the output file, followed by a new line
            f.write(prefix + filename + '\n')

