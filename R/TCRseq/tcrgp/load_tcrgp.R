
## Load in the TCRGP summary files
helsinki_predictions    <- fread("results/tcrgp/summary/meta/helsinki_predictions.txt")   %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
tumeh_predictions       <- fread("results/tcrgp/summary/meta/tumeh_predictions.txt")      %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
riaz_predictions        <- fread("results/tcrgp/summary/meta/riaz_predictions.txt")       %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
yusko_mnc_predictions   <- fread("results/tcrgp/summary/meta/yusko_mnc_predictions.txt")  %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget 
yusko_cd8_predictions   <- fread("results/tcrgp/summary/meta/yusko_cd8_predictions.txt")  %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
robert_predictions      <- fread("results/tcrgp/summary/meta/robert_predictions.txt")     %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget


## Correct riaz; a mislabeled based on prior IO status (they were the opposite)
riaz_correction <- fread("data/unselected_TCRseq/riaz_conversion.txt") %>% dplyr::rename(filename = oldfiles) %>% mutate(filename = extractFileName(filename)) # %>% mutate(filename = paste0("data/unselected_TCRseq/Riaz//", filename))

riaz_predictions <- riaz_predictions %>% left_join(riaz_correction)
riaz_predictions$io_stat[riaz_predictions$io_stat == "IO.naive"] <- "temp"
riaz_predictions$io_stat[riaz_predictions$io_stat == "Prior.IO"] <- "IO.naive"
riaz_predictions$io_stat[riaz_predictions$io_stat == "temp"] <- "Prior.IO"
riaz_predictions$filename <- extractFileName(riaz_predictions$newfiles)
riaz_predictions$newfiles <- NULL




## Emerson
emerson_predictions1   <- fread("results/tcrgp/summary/meta/emerson_predictions_1.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions2   <- fread("results/tcrgp/summary/meta/emerson_predictions_2.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions3   <- fread("results/tcrgp/summary/meta/emerson_predictions_3.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions4   <- fread("results/tcrgp/summary/meta/emerson_predictions_4.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions5   <- fread("results/tcrgp/summary/meta/emerson_predictions_5.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions6   <- fread("results/tcrgp/summary/meta/emerson_predictions_6.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions7   <- fread("results/tcrgp/summary/meta/emerson_predictions_7.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions8   <- fread("results/tcrgp/summary/meta/emerson_predictions_8.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions9   <- fread("results/tcrgp/summary/meta/emerson_predictions_9.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
emerson_predictions10  <- fread("results/tcrgp/summary/meta/emerson_predictions_10.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget

emerson_predictions     <- rbind(emerson_predictions1, emerson_predictions2, emerson_predictions3, emerson_predictions4, emerson_predictions5, emerson_predictions6, emerson_predictions7, emerson_predictions8, emerson_predictions9, emerson_predictions10)

act_til_predictions <- fread("results/tcrgp/summary/meta/act_til_predictions.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget
act_predictions     <- fread("results/tcrgp/summary/meta/act_predictions.txt") %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget


## TCGA
tcga_predictions      <- fread("results/tcrgp/summary/meta/tcga_beta.txt")     %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% correctMelanomaTarget


helsinki_predictions    <- fread("results/tcrgp/summary/meta/helsinki_predictions_new.txt") #  %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies()
tumeh_predictions       <- fread("results/tcrgp/summary/meta/tumeh_predictions.txt") 

