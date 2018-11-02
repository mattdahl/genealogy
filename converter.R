library(haven)

caseType <- c('Privacy')
caseTypesLabel <- caseType

stata_citation_data <- read_dta(file = paste('Data/', caseType, '/', caseType, 'CitationData.dta',sep=""))
stata_case_data <- read_dta(file = paste('Data/', caseType, '/', caseType, 'Cases.dta', sep = ''))
saveRDS(stata_citation_data, file = paste('Data/', caseType, '/', caseType, 'CitationData.rds', sep = ''))
saveRDS(stata_case_data, file = paste('Data/', caseType, '/', caseType, 'Cases.rds', sep = ''))