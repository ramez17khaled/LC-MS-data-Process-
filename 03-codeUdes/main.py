from tools import *

# INPUT FILES

inputPOSFile = "D:/data/MSDial/01.1-TermoData/1.2-NICOresultsPOS-thermo/AlignPOSData.csv"
inputNEGFile = "D:/data/MSDial/01.1-TermoData/1.2-NICOresultsNEG-thermo/AlignNEGData.csv"
metadataPOS = "D:/data/MSDial/02-Nico_metadata/TermoMetaData/thermoPOS.csv"
metadataNEG = "D:/data/MSDial/02-Nico_metadata/TermoMetaData/thermoNEG.csv"
outputPath = "D:/data/MSDial/05-codeOutput/Thermo_results/"
outputPOSFile = 'POS-manipulated'
outputNEGFile = 'NEG-manipulated'


## POS

#data import
dataPOS = read_file(inputPOSFile, sep = '\t')
metadataPOSData = read_file(metadataPOS, sep = ';', encoding='latin-1' )
#creating metabolite column start with M
dataPOS = add_metabolite_column(dataPOS)
excludDataPOS = ['Rt(min)', 'Mz']
fdataPOS = filter_column(dataPOS, exclude = excludDataPOS)
#transposing => metabolites as columns name
tdataPOS = transpose_data(fdataPOS, 'metabolite')
#exclud machine blc
exludRaw = ["sample_name == 'blc'","sample_name == 'blc_20240403164953'","sample_name == 'blc_20240404121923'", "sample_name == 'blc_20240404124839'", "sample_name == 'istd_ode'"]
ftdataPOS = filter_rows (tdataPOS,exclude=exludRaw)
#metadata merging 
excludDataPOS = ['id natif', 'class', 'injectionOrder']
metadataPOSData = filter_column(metadataPOSData,exclude=excludDataPOS)
metadataPOSData.set_index('sample_name', inplace = True)
start_dataPOS = merge_data(ftdataPOS, metadataPOSData)
# blank filtering 0 besed
blank_filtered_dataPOS= blank_filter(start_dataPOS)
# QC filtering: if 0 in 3/4 QC elimined
QC_blank_filterPOS = QC_filter_with_zeros (blank_filtered_dataPOS, 0.75)
# QCDil filtering CV <=10 elimined
raw_list = ["sample_name == '240326NCE_Globale_neg_QC3-DIL8'", "sample_name == '240326NCE_Globale_neg_QC3'" , "sample_name == '240326NCE_Globale_neg_QC3-DIL2'"]
QCDil_QC_blank_filterPOS = cv_filter (QC_blank_filterPOS, raw_list, 10, "<=")
# metadata merging for filtered data
QCDil_QC_blank_filterPOS = merge_data(QCDil_QC_blank_filterPOS, metadataPOSData)
# caving in a csv
output_pathPOS = save_as_csv(QCDil_QC_blank_filterPOS, output_dir=outputPath, output_file=outputPOSFile, file_conflict="replace")
print("CSV file saved at:", output_pathPOS)




## NEG

#data import
dataNEG = read_file(inputNEGFile, sep = ';')
metadataNEGData = read_file(metadataNEG, sep = ';', encoding='latin-1' )
#creating metabolite column start with M
dataNEG = add_metabolite_column(dataNEG)
excludDataNEG = ['Rt(min)', 'Mz']
fdataNEG = filter_column(dataNEG, exclude=excludDataNEG)
#transposing => metabolites as columns name
tdataNEG = transpose_data(fdataNEG, 'metabolite')
#exclud machine blc
exludRawNEG = ["sample_name == 'blc'","sample_name == 'blc_20240403164953'","sample_name == 'blc_20240404121923'", "sample_name == 'blc_20240404124839'", "sample_name == 'istd_ode'"]
ftdataNEG = filter_rows (tdataNEG,exclude=exludRawNEG)
#metadata merging 
excludDataNEG = ['id natif', 'class', 'injectionOrder']
metadataNEGData = filter_column(metadataNEGData,exclude=excludDataNEG)
metadataNEGData.set_index('sample_name', inplace = True)
start_dataNEG = merge_data(ftdataNEG, metadataNEGData)
# blank filtering 0 besed
blank_filtered_dataNEG= blank_filter(start_dataNEG)
# QC filtering: if 0 in 3/4 QC elimined
QC_blank_filterNEG = QC_filter_with_zeros (blank_filtered_dataNEG, 0.75)
# QCDil filtering CV <=10 elimined
raw_list = ["sample_name == '240326NCE_Globale_neg_QC3-DIL8'", "sample_name == '240326NCE_Globale_neg_QC3'" , "sample_name == '240326NCE_Globale_neg_QC3-DIL2'"]
QCDil_QC_blank_filterNEG = cv_filter (QC_blank_filterNEG, raw_list, 10, "<=")
# metadata merging for filtered data
QCDil_QC_blank_filterNEG = merge_data(QCDil_QC_blank_filterNEG, metadataNEGData)
# caving in a csv
output_pathNEG = save_as_csv(QCDil_QC_blank_filterNEG, output_dir=outputPath, output_file=outputNEGFile, file_conflict="replace")
print("CSV file saved at:", output_pathNEG)