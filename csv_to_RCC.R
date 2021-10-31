# this script aims to convert a csv file into a mnumber of RCC files. The number 
# of columns in the csv file is the same as that of the RCC files.
#
# Method: 
# 1. read the header information into a several arrays, each of which represents 
#    an attribute, e.g., file name, sample date, gene RLF, and Lane ID, etc
# 2. meanwhile, load the counts into a matrix and save the regular lines, e.g., <header> 
#    and </header>
# 3. as long as all counts and attributes are saved in the arrays and matrix, we generate 
#    a number of RCC files each of which corresponds to a column in the original csv file


# all parameters that need to be adjusted in the future
setwd("/fs/ess/PCON0022/liyang/Dr_Almaani/")


# load libraries
library(dplyr)
library(magrittr)
library(data.table)


# the function used to convert a csv file into a number of RCC files
# note that the last line should end with a "\n" symbol
convert_file <- function(csv.file) {
    
    # csv.file : the csv file
    
    count.data <- readLines(paste0("csv_files/", csv.file))
    # read all lines into the program
    
    # file attributes
    file.name <- grep('^File name', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the file names
    cat ("There are in total", length(file.name), "files.\n")
    
    sample.ID <- grep('^Sample ID', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the sample IDs
    cat ("There are in total", length(sample.ID), "samples.\n")
    
    if (length(file.name) != length(sample.ID)) {
        stop ("Error: the number of files is unequal with that of samples.\n")
    }
    
    file.version <- grep('^File Version', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the file versions
    # cat ("There are in total", length(file.version), "samples.\n")
    
    if (length(file.name) != length(file.version)) {
        stop ("Error: the number of files is unequal with that of file versions.\n")
    }
    
    gene.RLF <- grep('^GeneRLF', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the gene RLF
    cat ("There are in total", length(gene.RLF), "gene RLFs.\n")
    
    if (length(file.name) != length(gene.RLF)) {
        stop ("Error: the number of files is unequal with that of gene RLFs.\n")
    }
    
    
    # lane attributes
    lane.ID <- grep('^Lane ID', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the lane IDs
    cat ("There are in total", length(lane.ID), "lane IDs.\n")
    
    if (length(file.name) != length(lane.ID)) {
        stop ("Error: the number of files is unequal with that of lane IDs.\n")
    }
    
    FOV.count <- grep('^FOV Count', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the FOV count
    cat ("There are in total", length(lane.ID), "lane IDs.\n")
    
    if (length(file.name) != length(FOV.count)) {
        stop ("Error: the number of files is unequal with that of FOV counts.\n")
    }
    
    FOV.counted <- grep('^FOV Counted', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the FOV counted
    cat ("There are in total", length(FOV.counted), "FOV counted.\n")
    
    if (length(file.name) != length(FOV.counted)) {
        stop ("Error: the number of files is unequal with that of FOV counted.\n")
    }
    
    scanner.ID <- grep('^Scanner ID', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the scanner IDs
    cat ("There are in total", length(scanner.ID), "scanner IDs.\n")
    
    if (length(file.name) != length(scanner.ID)) {
        stop ("Error: the number of files is unequal with that of scanner IDs.\n")
    }
    
    stage.position <- grep('^StagePosition', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the stage position
    cat (length(stage.position), "\n")
    
    binding.density <- grep('^Binding Density', x = count.data, value = T) %>% `[[` (1) %>%
        strsplit(., split = ',') %>% `[[` (1) %>%
        gsub(., pattern = ' ', replacement = '-') %>%
        `[` (-(1:3)) # get the bidning density
    cat (length(binding.density), "\n")
    
    
    # Reporter Counts
    reporter.count <- count.data[-(1:(grep('^Reporter Counts', x = count.data, 
                                           value = F) + 1))] %>% 
        .[!grepl('^,', x = .)] %>% strsplit(., split = ",")
    # retain the lines for Reporter Counts
    
    cutted.count <- lapply(reporter.count, "[", (1:(length(file.name) + 3)))
    count.mm <- do.call(rbind, cutted.count)
    
    
    # generate the RCC files
    RCC.dir <- paste0("RCC_files/", gsub(".csv", "_RCC/", csv.file))
    dir.create(RCC.dir)
    
    for (i in seq_along(file.name)) {
        tmp.ll <- c("<Header>", 
                    paste0("FileVersion,", file.version[[i]]), 
                    "SoftwareVersion,2.1.2.3", 
                    "</Header>", 
                    "", 
                    "<Sample_Attributes>", 
                    paste0("ID,", sample.ID[[i]]),
                    "Owener,", 
                    "Comments,",
                    "Date,",
                    paste0("GeneRLF,", gene.RLF[[i]]), 
                    "SystemAPF,", 
                    "</Sample_Attributes>", 
                    "", 
                    "<Lane_Attributes>",
                    paste0("ID,", lane.ID[[i]]),
                    paste0("FovCount,", FOV.count[[i]]), 
                    paste0("FovCounted,", FOV.counted[[i]]),
                    paste0("ScannerID,", scanner.ID[[i]]), 
                    paste0("StagePosition,", stage.position[[i]]), 
                    paste0("BindingDensity,", binding.density[[i]]), 
                    "CartridgeID,", 
                    "</Lane_Attributes>", 
                    "", 
                    "<Code_Summary>", 
                    "CodeClass,Name,Accession,Count"
                    ) # an empty list
        
        for (j in 1:nrow(count.mm)) {
            line.ll <- c(count.mm[j, 1:3], count.mm[j, 3 + i])
            tmp.ll <- c(tmp.ll, paste(line.ll, collapse = ","))
        }
        
        RCC.file <- paste0(RCC.dir, paste0(file.name[[i]], ".RCC"))
        file.create(RCC.file)
        FD <- file(RCC.file)
        writeLines(tmp.ll, FD)
        close(FD)
    }
}


# main function
convert_file("Glom.csv")
convert_file("Tubular.csv")
