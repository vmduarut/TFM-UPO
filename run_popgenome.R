library("PopGenome")
# Cargar librería de PopGenome

setwd("D://tfm/data")
# Ajustar el directorio de trabajo a la carpeta local donde tenemos los datos

files_names <- vector()
# Creamos un vector vacío donde meteremos los nombres de los ficheros vcf que serán cargados posteriormente por el paquete

for (i in 1:22) {
  files_names <- c(files_names, paste("ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.split.AF_0.05.vcf.gz", sep = ""))
}
# Vamos metiendo los nombres de los ficheros con el bucle for

popIBS <- as.character(read.table("IBS_population.panel")[[1]])
popTSI <- as.character(read.table("TSI_population.panel")[[1]])
# Creamos dos vectores para las dos poblaciones que vamos a estudiar, a partir de los ficheros con los individuos de cada subpoblación que ya teníamos

j <- "1"
# Iniciamos una variable a 1, que será el número de cromosoma que irá aumentando en cada iteración del bucle

for (i in files_names) {
  # Cada iteración de este bucle será un análisis completo con cada fichero de cromosoma, así hasta completar los 22

  GENOME.class <- readVCF(i, numcols=10000, tid = j, from=1, to= 300000000, approx=FALSE, out="", parallel=FALSE, gffpath=FALSE)
  # Leemos el fichero VCF

  GENOME.class <- set.populations(GENOME.class, list(popIBS, popTSI), diploid = TRUE)
  # Ajustamos las dos poblaciones dentro del objeto GENOME

  GENOME.class.slide.300.50 <- sliding.window.transform(GENOME.class,300,50, type=2)
  GENOME.class.slide.600.50 <- sliding.window.transform(GENOME.class,600,50, type=2)
  GENOME.class.slide.1000.50 <- sliding.window.transform(GENOME.class,1000,50, type=2)
  GENOME.class.slide.2000.200 <- sliding.window.transform(GENOME.class,2000,200, type=2)
  GENOME.class.slide.3000.200 <- sliding.window.transform(GENOME.class,3000,200, type=2)
  GENOME.class.slide.5000.200 <- sliding.window.transform(GENOME.class,5000,200, type=2)
  # Creamos nuevos objetos para ir dividiendo el cromosoma original en fragmentos para el análisis en sliding windows, cada uno con diferentes parámetros

  GENOME.class.slide.300.50 <- neutrality.stats(GENOME.class.slide.300.50, FAST=TRUE)
  GENOME.class.slide.600.50 <- neutrality.stats(GENOME.class.slide.600.50, FAST=TRUE)
  GENOME.class.slide.1000.50 <- neutrality.stats(GENOME.class.slide.1000.50, FAST=TRUE)
  GENOME.class.slide.2000.200 <- neutrality.stats(GENOME.class.slide.2000.200, FAST=TRUE)
  GENOME.class.slide.3000.200 <- neutrality.stats(GENOME.class.slide.3000.200, FAST=TRUE)
  GENOME.class.slide.5000.200 <- neutrality.stats(GENOME.class.slide.5000.200, FAST=TRUE)
  # Llevamos a cabo el análisis de neutralidad para cada uno de los objetos de sliding windows creados previamente

  write.table(get.neutrality(GENOME.class.slide.300.50)[[1]], file = paste("chr", j, "_neutrality_stats_sliding_300_50_IBS.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.300.50)[[2]], file = paste("chr", j, "_neutrality_stats_sliding_300_50_TSI.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.600.50)[[1]], file = paste("chr", j, "_neutrality_stats_sliding_600_50_IBS.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.600.50)[[2]], file = paste("chr", j, "_neutrality_stats_sliding_600_50_TSI.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.1000.50)[[1]], file = paste("chr", j, "_neutrality_stats_sliding_1000_50_IBS.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.1000.50)[[2]], file = paste("chr", j, "_neutrality_stats_sliding_1000_50_TSI.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.2000.200)[[1]], file = paste("chr", j, "_neutrality_stats_sliding_2000_200_IBS.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.2000.200)[[2]], file = paste("chr", j, "_neutrality_stats_sliding_2000_200_TSI.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.3000.200)[[1]], file = paste("chr", j, "_neutrality_stats_sliding_3000_200_IBS.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.3000.200)[[2]], file = paste("chr", j, "_neutrality_stats_sliding_3000_200_TSI.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.5000.200)[[1]], file = paste("chr", j, "_neutrality_stats_sliding_5000_200_IBS.tsv", sep = ""), sep = "\t")
  write.table(get.neutrality(GENOME.class.slide.5000.200)[[2]], file = paste("chr", j, "_neutrality_stats_sliding_5000_200_TSI.tsv", sep = ""), sep = "\t")
  # Exportamos a ficheros .tsv los resultados de los análisis, en diferentes ficheros para la población IBS y para TSI

  j <- as.character((as.integer(j)) + 1)
  # Para ir aumentando el valor de esta variable tenemos que transformarla a entero, para luego volver a convertirlo en string, ya que los argumentos de readVCF sólo permiten strings

}



