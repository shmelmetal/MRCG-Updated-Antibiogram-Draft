library(tidyverse)

#read data
micro.raw <- read.csv("raw data/micro-longJun26.csv")

#vector of strings that represent a culture negative sample
neg.samples <- c(
  "No growth", "No growth after 48 hours",
  "No pathogens isolated", "No significant growth",
  "Salmonella and Shigella NOT isolated",
  "Salmonella and shigella not isolated",
  "Heavy growth", "mixed growth", "Bactec Negative",
  "NO ORGANISM SEEN", "NO ORGANISM SEEN WBC (+++) MAINLY POLYMPHS",
  
  # non-organism entries found in data review:
  "5", "AZOOSPERMIA", "CONTAMINATS", "ESBL Positive)", "Gram negative rods", "5.00",
  "NO FUNGAL ELEMENT ISOLATED", "contaminant)", "growing fungus",
  "no fungus isolated", "positive rods", "positive cocci")

#infection postive or negative
micro.raw <- micro.raw |> mutate(infection = ifelse(culture %in% neg.samples,
                                                    "Negative", "Positive")) |>
  relocate(infection, .after = "culture")

#change sampletype col to be less specific, fix specimen_class col, add MRSA col
micro.raw <- micro.raw |> 
  mutate(
    sampletype = 
      case_when(str_detect(sampletype, "swab") ~
                  "Swab",
                str_detect(sampletype, "aspirate") ~
                  "Aspirate",
                !str_detect(sampletype, "swab") & 
                  !str_detect(sampletype, "aspirate") ~
                  sampletype),
    
    MRSA = ifelse(culture_clean == 'Staphylococcus aureus' &
                    Cefoxitin == "R",
                  T, F)) |>
  mutate(specimen_class = ifelse(sampletype %in% 
                                   c("Aspirate", "Blood Culture",
                                     "CSF"), "Invasive", "Non-Invasive"))|>
  relocate(MRSA, .after = "ESBL")

#rename Amoxicillin to Amoxicillin-Clavulanate
micro.raw <- micro.raw |> rename(Amoxicillin_Clavulanate = Amoxicilin)

#new df with all positve samples
micro <- micro.raw[, c(1:4,23,24:26,29,33,34:57,59,60:64,67, 71)] |> 
  subset(infection == "Positive")

################################
# PART 1: Descriptive Analysis #
################################

#Sec A.Summary statistics for patient demographics 
#create more robust age grouping
micro <- micro |> mutate(age_group = 
                           case_when(Age <= 1 ~ "Neonates/Infants",
                                     Age <= 9 & 
                                       Age > 1 ~ "Children",
                                     Age <= 19 & 
                                       Age > 9 ~ "Adolescents",
                                     Age <= 59 & 
                                       Age > 19 ~ "Adults",
                                     Age > 60 ~ "Elders" ),
                         Sex = case_when(Sex == "F" ~ "Female",
                                         Sex == "M" ~ "Male",
                                         Sex == "" ~ "Unknown")) |>
  relocate(age_group, .after = "Age")

#age and sex summaries
pat.demo <- micro |> count(Sex, age_group) |>
  `colnames<-`(c("Sex", "Age Group", "n"))

pat.demo |> write_csv("tables/Patient_Demographics.csv")

#more summaries, bc I am not exactly sure what abdulrahman wanted
micro |> group_by(Sex) |> 
  summarise_at(vars(Age),
               list(min = min, median = median,
                    max = max, mean = mean,
                    sd = sd)) |>
  na.omit() |> write_csv("tables/Patient_Demographics_Summary.csv")


#barchart for patient count
pat.demo.bp <- pat.demo |> na.omit()|> 
  ggplot(aes(x = factor(`Age Group`,
                        levels = c("Neonates/Infants", "Children",
                                   "Adolescents", "Adults", "Elders")),
             y = n, fill = Sex)) +
  geom_col(position = "dodge", width = 0.75, color = "black") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("#FF8899","#5599FF")) + 
  theme_classic() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(vjust = -1),
        axis.title.x = element_text(vjust = -1),
        panel.grid.major.y = element_line(color  = "grey60")) +
  labs(x = "Age Group", y = "n")

#save
ggsave("figures/DescAnalysis-AgeSex.png",
       plot = pat.demo.bp,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)

#Sec B. Specimen Positivity Rates Table

#count infection by sample type from raw data
micro.raw |> group_by(sampletype) |> count(infection) |> 
  filter_out(sampletype == "")|> mutate(pct = round(n/sum(n)*100, 1)) |>
  `colnames<-`(c("Specimen Type", "Infection", "n", "%")) |>
  write_csv("tables/Specimen_Positivity_Rates.csv")

#Sec C. Specimen Positivity Rates Bar Chart
#first make DF of sampletypes and their respective totals for labellng
sampletypes.labs <- micro.raw |> count(sampletype) |> 
  mutate(label = paste0(sampletype, " ", "(n = ", n, ")"),
         labby = paste0( "(n = ", n, ")"))

#now prep for plot
spec.ispos.bp <- micro.raw |> group_by(sampletype) |> 
  count(infection, specimen_class) |> 
  left_join(sampletypes.labs[,-2], by= c("sampletype")) |>
  filter_out(sampletype == "") |> 
#plot
  ggplot(aes(x = sampletype, y = n, fill = infection)) + 
  geom_col(position = "fill", width = 0.75) +
  geom_text(aes(y = 1.025, label = labby)) + 
  theme_classic()+
  facet_wrap(~specimen_class, scale = "free_x") +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0,0.025)),
                     minor_breaks = seq(0,1,0.125)) + 
  scale_fill_manual(values = c("#32a852", "#a83232")) +
  labs(x = "Specimen Type", y = "Percent of Positive/Negative Specimens",
       fill = "Infection") + 
  theme(panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.y = element_line(color = "grey80"),)

spec.ispos.bp

#save
ggsave("figures/DescAnalysis-SpecimenIsPos2.png",
       plot = spec.ispos.bp,
       width = 10,
       height = 5,
       units = "in",
       device = "png",
       create.dir = T)

#Sec D. Pathogen Distribution by invasive and non-invasive
micro |> group_by(culture_clean) |> count(specimen_class) |> 
  filter_out(specimen_class == "") |>
  mutate(pct = round(n/sum(n)*100, 1)) |>
  write_csv("tables/DescAnalysis-InvasivePathogenDistribution.csv")

#Sec F, Count of all pathogens isolated
patho.count <- micro |> group_by(sampletype) |> 
  count(culture_clean,  specimen_class) |> 
  arrange(desc(n)) |>
  mutate(pct = n/sum(n)*100) 

#write into CSV for table
patho.count |> `colnames<-`(c("Culture", "Sample Type", "Specimen Class",
                              "n", "% (n/Sample Type Total)")) |>
  write_csv("tables/DescAnalysis-PathogenSampleType.csv")

#make an order vector, it will also serve as the significant paothgen vector
patho.signif <- micro |> count(culture_clean) |> 
  arrange(desc(n)) |> subset(n >= 20) |> select(culture_clean)|>
  unlist()|> unique() |> as.vector()


#pathogen count bar chart
patho.count.bp <- patho.count |> aggregate(n ~ culture_clean, sum) |>
  subset(culture_clean %in% patho.signif) |>
  ggplot(aes(x = factor(culture_clean, levels = patho.signif), y = n)) + 
  geom_col(fill = "white", color = "black") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) + 
  labs(x = "Culture")

#save
ggsave("figures/DescAnalysis-PathogenCount.png",
       plot = patho.count.bp,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)

#pathogens by sample type
patho.sampletype.bp <- patho.count |>
  mutate(culture_clean = 
           ifelse(culture_clean %in% patho.signif[-c(9:18)],
                  culture_clean, "Other"))|>
  filter_out(sampletype == "") |>
  aggregate(n ~ sampletype + culture_clean, sum) |>
  ggplot(aes(x = sampletype, y = n , 
             fill = factor(culture_clean, 
                           levels = c(patho.signif[-c(9:18)], "Other")))) +
  geom_col(position = "fill", width = 0.75) +
  scale_x_discrete(labels = sampletypes.labs[-1]) + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Sample Type", y = "% of Culture", fill = "Culture Isolated") +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 25, vjust = 1),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.y = element_line(color = "grey80"),
        legend.text = element_text(face = "italic"))

#save
ggsave("figures/DescAnalysis-PathogenSampleType.png",
       plot = patho.sampletype.bp,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)

################################
# PART 2: Antibiogram Analysis #
################################

#Sec A. Antibiogram of Key pathogens

#we aim for resistance rates of key pathogens
#start by pivoting long the micro df
micro.long <- micro[,c(1:3,10,13:41,42,43)] |> 
  pivot_longer(cols = Penicillin:imipenem,
               names_to = "Antibiotic",
               values_to = "Effectiveness") |>
  filter_out(Effectiveness == "" | is.na(Effectiveness)) |>
  #fix capitalization error in data
  mutate(Antibiotic = tools::toTitleCase(Antibiotic))

#Antibiotic classes col
micro.long <- micro.long |> 
  mutate(
    Antibiotic_Class = 
      case_when(Antibiotic %in% c("Ampicillin", "Amoxicilin", 
                                  "Amoxicillin_Clavulanate",
                                  "Penicillin","Cloxacillin", 
                                  "Oxacillin", "Methicillin") 
                ~ "Penicillins",
                Antibiotic %in% c("Cotrimoxazole") 
                ~ "Sulfonamides",
                Antibiotic %in% c("Gentamicin") 
                ~ "Aminoglycosides",
                Antibiotic %in% c("Chloramphenicol") 
                ~ "Peptidyl transferases",
                Antibiotic %in% c("Tetracycline") 
                ~ "Tetracycline",
                Antibiotic %in% c("Nitrofirantoin") 
                ~ "Nitrofurans",
                Antibiotic %in% c("Ciprofloxacin") 
                ~ "Fluoroquinolones",
                Antibiotic %in% c("Cefoxitin", "Cefotaxime", "Ceftazidime", 
                                  "Cefuroxime", "Ceftriaxone") 
                ~ "Cephalosporins",
                Antibiotic %in% c("Erythromycin")
                ~ "Macrolides",
                Antibiotic %in% c("Vancomycin") 
                ~ "Glycopeptides",
                Antibiotic %in% c("Meropenem", "Imipenem") 
                ~ "Carbapenems",
                Antibiotic %in% c("PolymixinB") 
                ~ "Polymixins"))

#now antibiogram
abx.overall <- micro.long |> 
  group_by(culture_clean, Antibiotic, Antibiotic_Class) |> 
  count(Effectiveness) |> 
  #remove blank culture 
  filter_out(culture_clean == "") |>
  #sum for year and calculate percent resistance
  mutate(abx_total = sum(n))  

#plot it on barchart
abx.overall.pp <- abx.overall |> 
  subset(culture_clean %in% patho.signif[-c(6:18)]) |>
  ggplot(aes(x = factor(culture_clean, levels = patho.signif[-c(6:18)]), 
             fill = Effectiveness, y = n)) +
  geom_col(position = "fill", width = 0.75) +
  facet_wrap(~Antibiotic) +
  geom_text(aes(y = 0.5, label = abx_total))+
  scale_y_continuous(labels = scales::percent,
                     expand = c(0,0)) +
  scale_fill_manual(values = c("#E1756F", "lightblue"),
                    labels = c("Resistant", "Susceptible")) +
  labs(x = "Culture", y = "n (%)")+ 
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90,
                                   vjust = 0.5,
                                   face = "italic"),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.y = element_line(color = "grey80"))

abx.overall.pp

ggsave("figures/antib-SPKE1.png",
       plot = abx.overall.pp,
       width = 12,
       height = 6,
       units = "in",
       device = "png",
       create.dir = T)

#sec B: rates of ESBL/AmpC/MRSA
micro |> select(culture_clean, year, MRSA, ESBL, AmpC)|>
  subset(culture_clean %in% 
           c("Staphylococcus aureus", "Pseudomonas aeruginosa",
             "Klebsiella pneumoniae", "Escherichia coli")) |>
  mutate(ESBL = ifelse(ESBL == T, 1,0),
         AmpC = ifelse(AmpC == T, 1,0),
         MRSA = ifelse(MRSA == T, 1,0)) |>
  filter_out(ESBL == 0 & AmpC == 0 & MRSA  == 0 | culture_clean == "") |>
  group_by(culture_clean, year) |>
  mutate(ESBL = sum(ESBL),
         AmpC = sum(AmpC),
         MRSA = sum(MRSA)) |> unique() |> 
  left_join(micro |> count(culture_clean, year), 
            by = c("culture_clean", "year")) |>
  mutate(`ESBL (%)` = round((ESBL / n)*100, 1),
         `AmpC (%)` = round((AmpC / n)*100, 1),
         `MRSA (%)` = round((MRSA / n)*100, 1)) |>
  relocate(`ESBL (%)`, .after = "ESBL") |>
  relocate(`AmpC (%)`, .after = "AmpC") |>
  relocate(`MRSA (%)`, .after = "MRSA") |>
  arrange(culture_clean, year) |>
  write_csv("tables/EAM_rates_SPKE_year.csv")

###################################
# PART 3: Temporal Trend Analysis #
###################################

#make antibiotic class resistance/sensitivity rate df by year for trends
abx.year <- micro.long |> 
  group_by(culture_clean, Antibiotic_Class, year) |> 
  count(Effectiveness) |> 
  #remove blank culture 
  filter_out(culture_clean == "") |>
  #sum for year and calculate percent resistance
  mutate(year_total = sum(n),
         pct = n / sum(n))  |>
  #filter for resistance and remove effectiness col
  filter(Effectiveness == "R") |> mutate(Effectiveness = NULL)

#point plot
abx.year.pp <- abx.year |> 
  #subset for specific  pathogens
  subset(culture_clean %in% c("Staphylococcus aureus", "Pseudomonas aeruginosa",
                              "Klebsiella pneumoniae", "Escherichia coli")) |>
  ggplot(aes(x = year, y = pct, color = culture_clean)) +
  geom_point() + geom_line() + 
  facet_wrap(~Antibiotic_Class) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Year", y = "% Resistant", color = "Culture") +
  theme_classic()

abx.year.pp

ggsave("figures/antib-SPKE-year.png",
       plot = abx.year.pp,
       width = 12,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)

#since year totals will really influence the data here I think it's important to 
#include a table with the plot to provide context, probably best as a supplement
abx.year |> 
  subset(culture_clean %in% 
           c("Staphylococcus aureus", "Pseudomonas aeruginosa",
             "Klebsiella pneumoniae", "Escherichia coli")) |>
  arrange(culture_clean, year) |> 
  mutate(pct = paste0(round(pct*100, 1), "%")) |>
  relocate(year, 1) |> 
  `colnames<-`(c("Year", "Culture", "Antibiotic Class", "n Resistant", 
                 "Year Total", "% Resistant")) |>
  write_csv("tables/AMR_rates_SPKE.csv")

###############################################
# PART 4: Multidrug Resistance (MDR) Analysis #
###############################################

#Sec A: MDR identification overall
#create df to determine a sample as MDR
#pivot long df into wide
micro.mdr <- micro.long |> 
  #change effectiveness to num format for antibiotic class aggregation
  mutate(Effectiveness = ifelse(Effectiveness == "R",
                                1,0)) |>
  aggregate(Effectiveness ~ opdno + serial + culture_clean + sampletype +
              Antibiotic_Class + MRSA + ESBL + AmpC + year, sum) |> 
  #now pivot wide
  pivot_wider(id_cols = c(opdno, serial, culture_clean, sampletype, 
                          MRSA, ESBL, AmpC, year),
              names_from = Antibiotic_Class, 
              values_from = Effectiveness)

#vector for the antibiotic classes for later selection
abx.classes <- micro.long$Antibiotic_Class |> unique()

#later is now, replace NA with 0 so the next mutate() works
micro.mdr[, abx.classes][is.na(micro.mdr[, abx.classes])] <- 0

#mutate in MDR alhamdullilah
micro.mdr <- micro.mdr |> 
  mutate(MDR = ifelse(rowSums(across(Aminoglycosides:Tetracycline) > 0) >= 3,
                      1, 0))

#count MDR and insances
micro.mdr.count <- micro.mdr |> aggregate(MDR ~ culture_clean, sum) |>
  filter_out(culture_clean == "") |>
  #now we join to general pathogen count for MDR %s
  left_join(aggregate(patho.count, n ~ culture_clean, sum)) |>
  #mutate percent cols and move them for the table later  
  mutate(`MDR %` = round(MDR/n*100, 1)) |>
  relocate(`MDR %`, .after = MDR) 

#write into CSV
micro.mdr.count |>
  arrange(desc(n)) |>
  `colnames<-`(c("Isolated Culture", "n of MDR Isolates",
                 "% of MDR Isolates", "Total n of Isolates")) |>
  write_csv("tables/MDR_rates_overall.csv")


#rates of MDR for signif pathogens
#pivot back to long for the plot
mdr.overall.bp <- micro.mdr.count[-c(3,5)] |>
  #cols for non MDR/XDR, to become their own set of rows later
  mutate(`non MDR` = n - MDR)|>
  pivot_longer(cols = c("MDR","non MDR"),
               names_to = "status",
               values_to = "count")|>
  #subset to include only significant pathogens
  subset(culture_clean %in% patho.signif[-c(6:18)]) |>
  #and plot
  ggplot(aes(x = factor(culture_clean, levels = patho.signif[-c(6:18)]), 
             y = count,
             fill = status)) +
  geom_col(position = "fill", width = 0.75) + 
  geom_text(aes(y = 0.5, label = n)) +
  scale_y_continuous(labels = scales::percent,
                     expand = c(0,0)) +
  scale_fill_manual(values = c("#E1756F", "lightblue"),
                    labels = c("MDR", "Non-MDR")) + 
  labs(y = "n (%)", x = "Isolated Culture", fill = "") + 
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 35, 
                                   face = "italic"),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.y = element_line(color = "grey80"))

mdr.overall.bp

#save
ggsave("figures/MDR/MDR-SPKE.png",
       plot = mdr.overall.bp,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)

#Sec B: MDR rates over time

micro.mdr.year <- micro.mdr |> aggregate(MDR ~ culture_clean + year, sum) |>
  filter_out(culture_clean == "") |>
  #now we join to general pathogen count for MDR and XDR %s
  left_join(micro |> count(year, culture_clean), 
            by = c("year", "culture_clean")) |>
  group_by(year, culture_clean) |>
  #mutate percent cols and move them for the table later  
  mutate(`MDR %` = round(MDR/n*100, 1)) |>
  relocate(`MDR %`, .after = MDR) |>
  relocate(year, 1) 

#write into CSV
micro.mdr.year |> 
  arrange(culture_clean) |>
  `colnames<-`(c("Year", "Isolated Culture", "n of MDR Isolates",
                 "% of MDR Isolates", "Yearly Total n of Isolates")) |>
  write_csv("tables/MDR_rates_year.csv")

#some manipulations for plotting
mdr.year.bp <- micro.mdr.year |>
  #mutate the percents to be decimals without rounding
  mutate(`MDR %` = MDR/n) |>
  #pivot longer for the plot to be more workable
  pivot_longer(cols = c(`MDR %`),
               names_to = "res_type", 
               values_to = "pct") |>
  #mutate out the % signs for better aesthetics
  mutate(res_type = str_remove_all(res_type, " %")) |>
  #subset to include only relevent pathogens
  subset(culture_clean %in% patho.signif[-c(6:18)]) |>
  #plot time
  ggplot(aes(x = year, y = pct, color = culture_clean))  +
  geom_line() + 
  geom_point()+
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,1)) +
  labs(x = "year", y = "% of MDR isolates", color = "Isolated Culture") + 
  theme_classic() + 
  theme(legend.text = element_text(face = "italic"),
        panel.grid.major.y = element_line(color = "grey50"),
        panel.grid.minor.y = element_line(color = "grey80"))

mdr.year.bp

#save
ggsave("figures/MDR/MDR-SPKE-year.png",
       plot = mdr.year.bp,
       width = 8,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)

##################################
# PART 5: Heat map Visualisation #
##################################

abx.overall.hm <- abx.overall |>
  #take only resistant isolates and those of a significant pathogen  
  subset(Effectiveness == "R" & 
           culture_clean %in% patho.signif[-c(6:18)]) |>
  #mutate a pct col for the fil, shorten Abx classes, col for geom_text() label
  mutate(pct = n/abx_total, 
         Antibiotic = (str_sub(Antibiotic, 1, 5)),
         Antibiotic_Class = str_sub(Antibiotic_Class, 1, 5),
         labby = paste0(round(pct*100, 1), "%" , "\n(", abx_total, ")")) |>
  #and we plot  
  ggplot((aes(x = Antibiotic, y = culture_clean, fill = pct))) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = labby),
            color = "black") +
  facet_grid(~Antibiotic_Class, scales = "free_x", space = "free_x") + 
  scale_fill_gradient(low = "lightblue", high = "#E1756F",
                      limits = c(0,1), labels = scales::percent) +
  scale_x_discrete(position = "top", expand = c(0.025,0.025)) +
  scale_y_discrete(expand = c(0.125,0.125)) + 
  labs(y = "Isolated Culture", x = "Antibiotic Class", 
       fill = "% Resistant") +
  theme_classic() +
  theme(axis.line = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(face = "italic"))

abx.overall.hm

#save
ggsave("figures/antib-SPKE-hm1.png",
       plot = abx.overall.hm,
       width = 16,
       height = 4,
       units = "in",
       device = "png",
       create.dir = T)
