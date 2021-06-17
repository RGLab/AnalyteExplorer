if (!dir.exists("data")) {
  dir.create("data")
}

if (!file.exists("data/gene_expression.rds")) {
  download.file("https://fh-pi-gottardo-r-eco-public.s3-us-west-2.amazonaws.com/AnalyteExplorer/gene_expression.rds", "gene_expression.rds")
}
if (!file.exists("data/blood_transcription_modules.csv")) {
  download.file("https://fh-pi-gottardo-r-eco-public.s3-us-west-2.amazonaws.com/AnalyteExplorer/blood_transcription_modules.csv", "blood_transcription_modules.csv")
}
if (!file.exists("data/gene_signatures.csv")) {
  download.file("https://fh-pi-gottardo-r-eco-public.s3-us-west-2.amazonaws.com/AnalyteExplorer/gene_signatures.csv", "gene_signatures.csv")
}

AnalyteExplorer:::log_message("Loading data...")

ge <- readRDS("data/gene_expression.rds")
ge_gene <- gene_expression[gene_expression$analyte_type == "gene", ]
ge_btm <- gene_expression[gene_expression$analyte_type == "blood transcription module", ]
ge_gs <- gene_expression[gene_expression$analyte_type == "gene signature", ]

btm <- read.csv("data/blood_transcription_modules.csv")
gs <- read.csv("data/gene_signatures.csv")
