### variables needed
version_miniAnalysis = 0.4.1

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
data_retrieval <- function(directory_ret, save_file, col_name_ampl = 'Peak.Amp..pA.', col_name_hz = 'Peak.to.Peak.Frequency..Hz.') {
  wb = createWorkbook()
  addWorksheet(wb, 'Settings')
  addWorksheet(wb, 'Amplitudes')
  addWorksheet(wb, 'Frequencies')

  files <- list.files(directory_ret, pattern="*.atf",full.names=TRUE)

  for (file in files) {

    txt_file = readLines(file)
    txt_file = txt_file[-c(1,2)]
    tabs = read.table(text=txt_file, sep="\t", header=TRUE)

    file_name = gsub("^.*/","",file)
    file_name = gsub("_events.atf","",file_name)

    tabs = type.convert(tabs, as.is = TRUE)

    addWorksheet(wb, file_name)

    col_num = which(colnames(tabs) == col_name_hz)
    tabs[,col_num] = as.numeric(tabs[,col_num])

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
freq_table <- function(file, sheet, col_name, n_bins, breaks, cut_off = FALSE) {

  # Load the data in the sheet
  data_sing_cell = read_excel(file, sheet = sheet)

  # create the frequency table and create an overflow bin
  data_sing_overflow = data_sing_cell
  col_num = which(colnames(data_sing_overflow) == col_name)

  data_sing_overflow_na_removed = data_sing_overflow[!is.na(data_sing_overflow[,col_num]),]

  data_sing_overflow = data_sing_overflow_na_removed

  if (as.numeric(unlist(data_sing_overflow[1,col_num])) < 0) {
    j = -1
  } else {
    j = 1
  }

  if (cut_off) {
    i = 1
    data_updated_cutoff = data.frame(data.frame(matrix(ncol = 1, nrow = 0)))
    for (row in 1:nrow(data_sing_overflow)) {
      if (!data_sing_overflow[row,col_num]*j > cut_off) {
        data_updated_cutoff[i,] = as.numeric(data_sing_overflow[row,col_num])
        i = i + 1
      }

    }

    data_sing_overflow = data_updated_cutoff[,1]*j

  } else {
    data_sing_overflow = ifelse(as.numeric(unlist(data_sing_overflow[,col_num]))*j>n_bins-1, n_bins, as.numeric(unlist(data_sing_overflow[,col_num]))*j) # n_bins is overflow
  }

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
plot_cdf <- function(cdf, plot_over_plot = TRUE, second_axis = TRUE) {

  if (plot_over_plot) {
    # make sure it overlaps the histogram
    par(new = TRUE)
  }
  par(mar = c(5, 5, 5, 5))

  # plot the cumulative distribution in red and add the second y-axis
  plot(x = cdf, type = 'l', col = "red",
       axes = !second_axis, xlab = "", ylab = "")

  if (second_axis) {
    axis(side = 4, at = seq(0,100,10))
    mtext("Cumulative distribution", side = 4, line = 3, outer=FALSE)
  }

}

# create final plot
plot_hist_cdf <- function(sheet, output, breaks, name = NULL) {

  create_hist(sheet, output, breaks, name)

  cdf = calc_cdf(output)

  plot_cdf(cdf)

}

# create mean row for bins

create_mean_bins <- function(data, genotype, pharmacology) {

  mean_info = c('Mean', genotype, pharmacology)

  btwn_stp = data[data$Genotype == genotype & data$Pharmacology == pharmacology, 4:ncol(data)]

  mean_bins = colMeans(btwn_stp)
  mean_bins = append(mean_bins, mean_info, after=0)

  return(mean_bins)

}


# plot all cumulative distributions
plot_cumdistr <- function(data, sheet, n_bins, size_bins) {

  mean_cums = type.convert(data[data$Cell == 'Mean',], as.is = TRUE)

  order = paste(mean_cums$Genotype, mean_cums$Pharmacology)
  color = c('black','blue','red','green')

  if (sheet == 'Amplitudes') {
    xlabs = 'Amplitude (pA)'
    legend_x = 55
    legend_y = 30
  } else if (sheet == 'Frequencies') {
    xlabs = 'Frequency (Hz)'
    legend_x = 11
    legend_y = 30
  } else {
    stop('"sheet" needs to be either "Amplitudes" or "Frequencies"')
  }

  plot.new()
  # par(mar = rep(1,4))
  i = 1
  for (row in seq(1,nrow(mean_cums))) {

    cdf = calc_cdf(mean_cums[row,4:ncol(mean_cums)])

    if (i==1) {
      plot(seq(0,n_bins-size_bins,size_bins), cdf, xlab = xlabs, ylab='Cumulative percentage',
           type = 'l', col = color[i],
           xaxp = c(0,100,5/size_bins))
    } else {
      lines(seq(0,n_bins-size_bins,size_bins), cdf,col = color[i])
    }

    i = i + 1

  }

  legend(legend_x, legend_y, legend=order,
         col=color, lty=1)

}

# plot the amplitude and freqency histograms/cdfs per cell, condition and overlapping plots
plot_histcdf_per_cell <- function(data_overview, file_combined, sheets, new_file_plots, col_name_ampl, bins_ampl, col_name_hz, bins_hz, cut_off_ampl = FALSE, cut_off_hz = FALSE, exclude_manual = FALSE, exclude_auto = FALSE, break_point = -1) {

  wb = loadWorkbook(file_combined) # load workbook

  n_bins_ampl = bins_ampl[length(bins_ampl)]
  n_bins_hz = bins_hz[length(bins_hz)]

  i = 1
  for (sheet in sheets) {

    if (i==break_point){
      break
    }

    if (sheet == 'Amplitudes' || sheet == 'Frequencies' || sheet == 'Settings') {
      i = 0

    } else {

      overview_row = data_overview[data_overview$FileName == sheet,]
      genotype = overview_row$Genotype
      pharma = overview_row$Pharmacology

      if (is.na(pharma)) {
        pharma = 'ACSF'
      }

      ##### ampl
      output = freq_table(file_combined, sheet, col_name_ampl, n_bins_ampl, bins_ampl, cut_off = cut_off_ampl)

      info_cells = c(sheet, genotype, pharma)
      output = append(output, info_cells, after=0)

      data_freqs_ampl[i,] = output
      output = as.numeric(output[-c(1:length(info_cells))])

      plot_hist_cdf(sheet, output, bins_ampl, name = 'Amplitudes')

      # insert the plot in de correct sheet in cell AI2 (35,1) (that is next to the tables)
      insertPlot(wb, sheet, width = 6, height = 4, xy = NULL, startRow = 1,
                 startCol = 1, fileType = "png", units = "in", dpi = 300)

      ##### Hz

      output = freq_table(file_combined, sheet, col_name_hz, n_bins_hz, bins_hz, cut_off = cut_off_hz)

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

  # exclusion
  if (exclude_manual & exclude_auto) {
    exclusion_manual = data_overview[which(data_overview$Exclude == 'x'),'FileName']
    exclusion_auto = data_overview[which(data_overview$`Automatic Exclusion` == 'x'),'FileName']
    data_freqs_ampl = data_freqs_ampl[!(data_freqs_ampl$Cell %in% as.list(exclusion_manual)$FileName) &
                                        !(data_freqs_ampl$Cell %in% as.list(exclusion_auto)$FileName),]
    data_freqs_hz = data_freqs_hz[!(data_freqs_hz$Cell %in% as.list(exclusion_manual)$FileName) &
                                    !(data_freqs_hz$Cell %in% as.list(exclusion_auto)$FileName),]

  } else {

    if (exclude_manual) {
      exclusion_manual = data_overview[which(data_overview$Exclude == 'x'),'FileName']
      data_freqs_ampl = data_freqs_ampl[!(data_freqs_ampl$Cell %in% as.list(exclusion_manual)$FileName),]
      data_freqs_hz = data_freqs_hz[!(data_freqs_hz$Cell %in% as.list(exclusion_manual)$FileName),]
    }
    if (exclude_auto) {
      exclusion_auto = data_overview[which(data_overview$`Automatic Exclusion` == 'x'),'FileName']
      data_freqs_ampl = data_freqs_ampl[!(data_freqs_ampl$Cell %in% as.list(exclusion_auto)$FileName),]
      data_freqs_hz = data_freqs_hz[!(data_freqs_hz$Cell %in% as.list(exclusion_auto)$FileName),]
    }

  }


  starting_col = 1
  for (genotype in genotypes) {

    starting_row = 32

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

  n_bins_ampl = bins_ampl[length(bins_ampl)]
  n_bins_hz = bins_hz[length(bins_hz)]

  size_bins_ampl = bins_ampl[2]-bins_ampl[1]
  size_bins_hz = bins_hz[2]-bins_hz[1]

  plot_cumdistr(data_freqs_ampl, 'Amplitudes', n_bins_ampl, size_bins_ampl)
  insertPlot(wb, 'Amplitudes', width = 8, height = 6, xy = NULL, startRow = 1,
             startCol = 1, fileType = "png", units = "in", dpi = 300)

  plot_cumdistr(data_freqs_hz, 'Frequencies', n_bins_hz, size_bins_hz)
  insertPlot(wb, 'Frequencies', width = 8, height = 6, xy = NULL, startRow = 1,
             startCol = 1, fileType = "png", units = "in", dpi = 300)

  data_freqs_ampl = type.convert(data_freqs_ampl, as.is = TRUE)
  data_freqs_hz = type.convert(data_freqs_hz, as.is = TRUE)

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

  # add settings information
  date = Sys.Date()
  version_miniAnalysis = version_miniAnalysis
  cut_off_frequency = cut_off_hz
  cut_off_amplitudes = cut_off_ampl

  writeData(wb, 'Settings', 'Date', 1, 1)
  writeData(wb, 'Settings', date, 2, 1)

  writeData(wb, 'Settings', 'Version', 1, 2)
  writeData(wb, 'Settings', version_miniAnalysis, 2, 2)

  writeData(wb, 'Settings', 'Amplitudes', 2, 4)
  writeData(wb, 'Settings', 'Frequencies', 3, 4)

  writeData(wb, 'Settings', 'Number of bins', 1, 5)
  writeData(wb, 'Settings', n_bins_ampl, 2, 5)
  writeData(wb, 'Settings', n_bins_hz, 3, 5)

  writeData(wb, 'Settings', 'Size of bins', 1, 6)
  writeData(wb, 'Settings', size_bins_ampl, 2, 6)
  writeData(wb, 'Settings', size_bins_hz, 3, 6)

  writeData(wb, 'Settings', 'Cut off', 1, 7)
  writeData(wb, 'Settings', cut_off_ampl, 2, 7)
  writeData(wb, 'Settings', cut_off_hz, 3, 7)

  writeData(wb, 'Settings', 'Manual exclusion', 1, 9)
  writeData(wb, 'Settings', 'Automated exclusion', 1, 10)
  writeData(wb, 'Settings', exclude_manual, 2, 9)
  writeData(wb, 'Settings', exclude_auto, 2, 10)

  writeData(wb, 'Settings', 'Dates of patching', 1, 12)
  writeData(wb, 'Settings', 'First day', 1, 13)
  writeData(wb, 'Settings', as.numeric(min(data_overview$Date)), 2, 13)
  writeData(wb, 'Settings', 'Last day', 1, 14)
  writeData(wb, 'Settings', as.numeric(max(data_overview$Date)), 2, 14)

  # save workbook
  saveWorkbook(wb, new_file_plots, overwrite = TRUE)

  print("Analysis is done")
  # return(data_freqs_ampl)

}

# add conditional formatting to insanely high frequencies
add_conditional_format <- function(file_plots, new_file, cut_off = 20) {

  wb = loadWorkbook(file_plots)
  sheets = readxl::excel_sheets(file_plots)

  negStyle <- createStyle(fontColour = "#ff0000", bgFill = "#eb7a7a")

  for (sheet in sheets) {

    if (!sheet == 'Amplitudes' & !sheet == 'Frequencies' & !sheet == 'Settings') {

      data_sheet = read.xlsx(file_plots, sheet)
      print(sheet)
      rows_cut_off = which(data_sheet[,29] > as.numeric(cut_off))

      print(data_sheet)
      print(1:ncol(data_sheet))

      addStyle(wb, sheet, negStyle, rows_cut_off + 1, cols=1:ncol(data_sheet), gridExpand = TRUE)

    }

  }

  saveWorkbook(wb, new_file, TRUE)

  print("Formatting completed")

}
