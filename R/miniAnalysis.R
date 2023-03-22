### functions

# calc peak-peak freq
calc_peakpeak_freq <- function(table, col_name = 'Peak.to.Peak.Frequency..Hz.') {
  
  col_num = which(colnames(table) == col_name)
  start_times = table[,col_num]
  
  previous = NA
  instant_freqs = list()
  
  for (start_time in start_times) {
    
    start_time = as.numeric(start_time)
    
    if (!is.na(previous) & !is.nan(previous) & previous != 'N/A'){
      
      instant_freq = 1 / ((start_time - previous) / 1000)
      
    } else {
      
      instant_freq = NA
      
    }
    
    instant_freqs = append(instant_freqs, instant_freq)
    previous = start_time
    
  }
  
  return(instant_freqs)
  
}
# data retrieval and saving
data_retrieval <- function(directory_ret, save_file) {
  wb = createWorkbook()
  addWorksheet(wb, 'Amplitudes')
  addWorksheet(wb, 'Frequencies')
  
  files <- list.files(directory_ret, pattern="*.atf",full.names=TRUE)
  
  for (file in files) {
    
    txt_file = readLines(file)
    txt_file = txt_file[-c(1,2)]
    tabs = read.table(text=txt_file, sep="\t", header=TRUE)
    
    file_name = gsub("^.*/","",file)
    file_name = gsub("_events.atf","",file_name)
    
    
    addWorksheet(wb, file_name)
    
    writeDataTable(
      wb=wb, sheet=file_name,
      x=tabs[,1:33],
      startCol = 1,
      startRow = 1
    )
    
  }
  
  saveWorkbook(wb, file = save_file, overwrite = TRUE)
  
  print("Data retrieved and saved.")
  
}

# create empty dataframe with bins
bins_dataframe <- function(n_bins, size_bin) {
  
  info_cells = c('Cell', 'Genotype', 'Pharmacology')
  
  bins = seq(0,n_bins,by=size_bin)
  ranges = paste(head(bins,-1), bins[-1], sep=" - ")
  
  data_bins = data.frame(matrix(ncol = length(bins)+length(info_cells)-1, nrow = 0))
  
  names_cols = append(ranges, info_cells,after=0)
  names_cols[length(names_cols)] = sprintf('> %f', bins[length(bins)]-size_bin)
  colnames(data_bins) <- names_cols
  
  return(list("dataframe" = data_bins, "bins" = bins))
  
}

# return frequency table
freq_table <- function(file, sheet, col_name, n_bins, breaks) {
  
  # Load the data in the sheet
  data_sing_cell = read_excel(file, sheet = sheet)
  
  # create the frequency table and create an overflow bin
  data_sing_overflow = data_sing_cell
  col_num = which(colnames(data_sing_overflow) == col_name)
  
  data_sing_overflow_na_removed = data_sing_overflow[data_sing_overflow[,col_num] != 'N/A',]
  
  data_sing_overflow = data_sing_overflow_na_removed
  
  if (as.numeric(unlist(data_sing_overflow[1,col_num])) < 0) {
    j = -1
  } else {
    j = 1
  }
  
  data_sing_overflow = ifelse(as.numeric(unlist(data_sing_overflow[,col_num]))*j>n_bins-1, n_bins, as.numeric(unlist(data_sing_overflow[,col_num]))*j) # n_bins is overflow
  
  # create the frequency table in the freq variable.
  freq = hist(data_sing_overflow, breaks=breaks, include.lowest=TRUE, plot=FALSE)
  
  output = as.numeric(freq$counts)
  
  return(output)
  
}

# histogram
create_hist <- function(sheet, freqs, breaks, name = NULL) {  
  # do not overlap the histograms
  par(new = FALSE)
  par(mar = c(5, 5, 5, 5))
  # the counts variable is the numeric version of the frequencies
  counts = as.numeric(freqs)
  
  # create the histogram of the current recording and plot it
  newhist = list(breaks=breaks, counts=counts, density=counts/diff(breaks),
                 xname=paste(sheet,name))
  class(newhist) = "histogram"
  plot(newhist)
  
}

# cdf
calc_cdf <- function(counts) {
  # add cumulative distribution
  cdf = c(0)
  
  # The total frequencies needs to be known to use it as a percentage
  total_freqs = sum(counts)
  
  # calculate the cumulative percentage
  for (bin in counts) {
    cdf = append(cdf,cdf[length(cdf)]+bin/total_freqs * 100)
  }
  cdf = cdf[-1]
  return(cdf)
  
}

# cdf plot
plot_cdf <- function(cdf, breaks) {
  
  # make sure it overlaps the histogram
  par(new = TRUE)
  par(mar = c(5, 5, 5, 5))
  
  # plot the cumulative distribution in red and add the second y-axis
  plot(x = cdf, type = 'l', col = "red",
       axes = FALSE, xlab = "", ylab = "")
  axis(side = 4, at = seq(0,100,10))
  mtext("Cumulative distribution", side = 4, line = 3, outer=FALSE)
  
}

# create final plot 
plot_hist_cdf <- function(sheet, output, breaks, name = NULL) {
  
  create_hist(sheet, output, breaks, name)
  
  cdf = calc_cdf(output)
  
  plot_cdf(cdf, breaks)
  
}

# create mean row for bins

create_mean_bins <- function(data, genotype, pharmacology) {
  
  mean_info = c('Mean', genotype, pharmacology)
  
  btwn_stp = data[data$Genotype == genotype & data$Pharmacology == pharmacology, 4:ncol(data)]
  
  mean_bins = colMeans(btwn_stp)    
  mean_bins = append(mean_bins, mean_info, after=0)
  
  return(mean_bins)
  
}

plot_histcdf_per_cell <- function(data_overview, file_combined, sheets, new_file_plots, col_name_ampl, bins_ampl, col_name_hz, bins_hz, break_point = -1) {
  
  wb = loadWorkbook(file_combined) # load workbook
  
  n_bins_ampl = bins_ampl[length(bins_ampl)]
  n_bins_hz = bins_hz[length(bins_hz)]
  
  i = 1
  for (sheet in sheets) {
    
    if (i==break_point){
      break
    }
    
    if (sheet == 'Amplitudes' || sheet == 'Frequencies') {
      i = 0
      
    } else {
      
      overview_row = data_overview[data_overview$FileName == sheet,]
      genotype = overview_row$Genotype
      pharma = overview_row$Pharmacology
      
      if (is.na(pharma)) {
        pharma = 'ACSF'
      }
      
      
      ##### ampl
      output = freq_table(file_combined, sheet, col_name_ampl, n_bins_ampl, bins_ampl)
      
      info_cells = c(sheet, genotype, pharma)
      output = append(output, info_cells, after=0)
      
      data_freqs_ampl[i,] = output
      output = as.numeric(output[-c(1:length(info_cells))])
      
      plot_hist_cdf(sheet, output, bins_ampl, name = 'Amplitudes')
      
      # insert the plot in de correct sheet in cell AI2 (35,1) (that is next to the tables)
      insertPlot(wb, sheet, width = 6, height = 4, xy = NULL, startRow = 1, 
                 startCol = 1, fileType = "png", units = "in", dpi = 300)
      
      ##### Hz
      
      output = freq_table(file_combined, sheet, col_name_hz, n_bins_hz, bins_hz)
      
      info_cells = c(sheet, genotype, pharma)
      
      output = append(output, info_cells, after=0)
      
      data_freqs_hz[i,] = output
      output = as.numeric(output[-c(1:length(info_cells))])
      
      plot_hist_cdf(sheet, output, bins_hz, name = 'Frequencies')
      
      insertPlot(wb, sheet, width = 6, height = 4, xy = NULL, startRow = 1, 
                 startCol = 10, fileType = "png", units = "in", dpi = 300)
      
    }
    
    i = i + 1
    
  }
  
  # add mean row per condition
  # Conditions:
  # - Genotype
  # - Pharmacology
  
  genotypes = unique(data_freqs_ampl$Genotype)
  pharmas = unique(data_freqs_ampl$Pharmacology)
  # pharmacologies = c('ACSF','DMSO','ISO + IBMX')
  
  starting_col = 1
  for (genotype in genotypes) {
    
    starting_row = 1
    
    for (pharma in pharmas) {
      
      pharma = toString(pharma)
      data_freqs_ampl = type.convert(data_freqs_ampl, as.is = TRUE)
      data_freqs_hz = type.convert(data_freqs_hz, as.is = TRUE)
      
      mean_ampl = create_mean_bins(data_freqs_ampl, genotype, pharma)
      mean_hz = create_mean_bins(data_freqs_hz, genotype, pharma)
      
      if(!any(mean_ampl == 'NaN')){
        
        data_freqs_ampl[nrow(data_freqs_ampl)+1,] = mean_ampl
        
        plot_hist_cdf('Amplitudes', as.numeric(mean_ampl[-c(1,2,3)]), bins_ampl, name = paste(genotype,pharma))
        insertPlot(wb, 'Amplitudes', width = 6, height = 4, xy = NULL, startRow = starting_row, 
                   startCol = starting_col, fileType = "png", units = "in", dpi = 300)
        
      }
      
      if(!any(mean_hz == 'NaN')) {
        
        data_freqs_hz[nrow(data_freqs_hz)+1,] = mean_hz
        
        plot_hist_cdf('Frequencies', as.numeric(mean_hz[-c(1,2,3)]), bins_hz, name = paste(genotype,pharma))
        insertPlot(wb, 'Frequencies', width = 6, height = 4, xy = NULL, startRow = starting_row, 
                   startCol = starting_col, fileType = "png", units = "in", dpi = 300)
        
      }
      
      starting_row = starting_row + 21
      
    }
    
    starting_col = starting_col + 8
  }
  
  
  # add means to main sheets
  writeDataTable(
    wb=wb, sheet='Amplitudes',
    x=data_freqs_ampl,
    startCol = 1,
    startRow = 1
  )
  
  writeDataTable(
    wb=wb, sheet='Frequencies',
    x=data_freqs_hz,
    startCol = 1,
    startRow = 1
  )
  
  # save workbook
  saveWorkbook(wb, new_file_plots, overwrite = TRUE)
  
  print("Analysis is done")
  
}