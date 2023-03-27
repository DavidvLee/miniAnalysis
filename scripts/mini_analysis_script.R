### Install and load the libraries needed (DO NOT EDIT)
if("openxlsx" %in% rownames(installed.packages()) == FALSE) {
  install.packages("openxlsx")
}
if("readxl" %in% rownames(installed.packages()) == FALSE) {
  install.packages("readxl")
}
if("devtools" %in% rownames(installed.packages()) == FALSE) {
  install.packages("devtools")
}
library(devtools)
if("miniAnalysis" %in% rownames(installed.packages()) == FALSE) {
  install_github("DavidvLee/miniAnalysis")
}

library(openxlsx)
library(readxl)
library(miniAnalysis)


##### Set Parameters (EDIT)
# files and folders
file_combined = 'events_combined.xlsx' # the file where all the atf files data will be saved
new_file_plots = 'plots.xlsx' # the file where the combined data + plots will be saved
directory_EventsToXlsx = "C:/Users/David/OneDrive - UvA/Jaar 3/4-6 internship/analysis David van Lee updated/3_events/all (used for flow analysis)/" # folder with all the atf files
file_overview = "C:/Users/David/OneDrive - UvA/Jaar 3/4-6 internship/analysis David van Lee updated/analysis_event_statistics.xlsm" # the file with the general overview of all the recordings (event statistics file)

# bins amplitude
size_bins_ampl = 1 # size of bins
n_bins_ampl = 100 + size_bins_ampl # number of bins (e.g., 100 means 0 - 100)
col_name_ampl = 'Peak.Amp..pA.' # name of the column of the amplitudes in the combined file

# bins frequency
size_bins_hz = 0.1
n_bins_hz = 20 + size_bins_hz
col_name_hz = 'Peak.to.Peak.Frequency..Hz.'


##### DO NOT EDIT

# combines the atf files if needed.
data_retrieval(directory_EventsToXlsx, file_combined)

# load the sheet names of the file
sheets = readxl::excel_sheets(file_combined)

# dataframes and bins amplitude and frequencies
data_ampl = bins_dataframe(n_bins_ampl, size_bins_ampl)
data_freqs_ampl = data_ampl$dataframe
bins_ampl = data_ampl$bins

data_hz = bins_dataframe(n_bins_hz, size_bins_hz)
data_freqs_hz = data_hz$dataframe
bins_hz = data_hz$bins

# Load the overview file
data_overview = read_excel(file_overview, sheet = 'Overview')
data_overview = data_overview[complete.cases(data_overview$FileName),]


##### Creates the plots per cell and the means
plot_histcdf_per_cell(data_overview = data_overview,
                      file_combined = file_combined,
                      sheets = sheets,

                      new_file_plots = new_file_plots,

                      col_name_ampl = col_name_ampl,
                      bins_ampl = bins_ampl,

                      col_name_hz = col_name_hz,
                      bins_hz = bins_hz,

                      break_point = -1)
