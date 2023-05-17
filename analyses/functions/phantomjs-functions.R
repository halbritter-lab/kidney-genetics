#### This file holds analyses functions for downloading web urls using phantomjs

#' Download Web URL Files Using PhantomJS
#'
#' This function uses PhantomJS to download files from a specified URL. The downloaded
#' file is saved with a specified basename, creation date, and extension. The path to save
#' the file can also be specified.
#' 
#' @param input_url A character string. The URL to download file from.
#' @param output_basename A character string. The base name for the output file.
#' @param output_extension A character string. The extension for the output file.
#' Default: "html"
#' @param save_path A character string. The path to save the downloaded file.
#' Default: "data/downloads"
#'
#' @return A character string. The path of the downloaded file.
#' 
#' @examples
#' download_url_by_phantomjs("https://www.example.com", "example_file")
#'
#' @export
download_url_by_phantomjs <- function(input_url,
		output_basename,
		output_extension = "html",
		save_path = "data/downloads") {
	# based on: https://ladal.edu.au/webcrawling.html

	pjs_instance <- run_phantomjs(timeout = 5000)
	pjs_session <- Session$new(port = pjs_instance$port)

	# go to input URL
	pjs_session$go(input_url)
	# render page
	rendered_source <- pjs_session$getSource()

	# generate file name
	creation_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")
	output_filename <- paste0(save_path, "/", output_basename, ".", creation_date, ".", output_extension)
	
	cat(rendered_source, file = output_filename)
	
	return(output_filename)
}