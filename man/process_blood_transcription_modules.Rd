% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blood_transcription_modules.R
\name{process_blood_transcription_modules}
\alias{process_blood_transcription_modules}
\title{Create a metadata table for blood transcription modules}
\usage{
process_blood_transcription_modules()
}
\value{
data.table
}
\description{
Updates the gene symbols with Hugo in blood transcription modules table from
Li S, et al. This table is used to display metadata about the selected module
in the app app.
}
\examples{
\dontrun{
btm <- process_data("blood_transcription_modules")
validate(btm)
res <- update_table(btm)
}
}
\references{
Li S, Rouphael N, Duraisingham S, et al.
Molecular signatures of antibody responses derived from a systems biology
study of five human vaccines.
Nat Immunol. 2014;15(2):195-204. doi:10.1038/ni.2789
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3946932/bin/NIHMS540680-supplement-26.zip}
}
