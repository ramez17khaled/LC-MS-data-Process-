import subprocess
import os

# Path to your batch file
batch_file_path = r"C:\Users\mo-lipidomique.i2mc\Desktop\MSDIAL ver.4.9.221218 Windowsx64\MsdialConsoleApp.exe"

##POS
# Arguments for your batch file
argumentsPOS = [
    "lcmsdda",
    "-i", r".\01-mzmlNEG-thermo\",
    "-o", r".\1.2-NICOresultsNEG-thermo\",
    "-m", r".\01-mzmlNEG-thermo\parametre.txt"
]

# Execute the batch file with arguments
subprocess.call([batch_file_path] + argumentsPOS)

##NEG
# Arguments for your batch file
argumentsNEG = [
    "lcmsdda",
    "-i", r".\01-mzmlNEG-thermo\",
    "-o", r".\1.2-NICOresultsNEG-thermo\",
    "-m", r".\01-mzmlNEG-thermo\parametre.txt"
]

# Execute the batch file with arguments
subprocess.call([batch_file_path] + argumentsNEG)