get_data <- function(fileLocation) {
    # load every file in fileLocation. 
    # Expected to find  1 RData containing the dataframe
    #                   1 ShpFile containing the shape 
    
    allFiles = dir()
    shpFile = allFiles[endsWith(allFiles,".shp")]
}