# Revgadgets

###############################
# ver prios del Yule
###############################

#Cargar librerías necesarias
library(RevGadgets)
library(ggplot2)

# Leer los datos del MCMC
mcmc_trace <- readTrace("/home/daniela/GSM/Revbayes/Rev_out/98/98_2701_diversification_Yule.log")

# Visualizar la distribución posterior de birth_rate
plotTrace(mcmc_trace, vars="birth_rate") # vars="birth_rate"

# Calcular la media posterior y el HPD del 95%
summary_stats <- summarizeTrace(mcmc_trace, vars="birth_rate")
print(summary_stats)


# Leer los datos del MCMC combinando las dos cadenas
posterior_trace <- readTrace(c("/home/daniela/GSM/Revbayes/Rev_out/98/98_2701_diversification_Yule_run_1.log", "/home/daniela/GSM/Revbayes/Rev_out/98/98_2701_diversification_Yule_run_2.log"))

# Extraer el primer conjunto de muestras (lista de data frames)
yule_posterior <- posterior_trace[[1]]

# Simular 10,000 valores de la distribución previa
yule_prior <- data.frame(birth_rate = runif(10000, min=0, max=100))

# Agregar la columna de la distribución previa en el posterior
yule_posterior$birth_rate_prior <- sample(yule_prior$birth_rate, size = nrow(yule_posterior), replace = TRUE)

# Graficar la comparación
plotTrace(list(yule_posterior), vars = c("birth_rate", "birth_rate_prior"))[[1]] +
  theme(legend.position = c(0.80, 0.80),  # Ubicación de la leyenda
        legend.text = element_text(size =20),  # Tamaño del texto de la leyenda
        legend.title = element_text(size = 22)) +  # Tamaño del título de la leyenda
  xlim(0, 1)  # Ajusta el límite según los valores observados

################################################################################
#
# Plot estimates of the birth-death model
#
# authors: Sebastian Höhna
#
################################################################################

library(RevGadgets)
library(ggplot2)

# read the posterior and prior output
bd_posterior <- readTrace("/home/daniela/GSM/Revbayes/Rev_out/98/98_2701_diversification_BD.log")

# plot the prior vs the posterior
plot <- plotTrace(bd_posterior, vars=c("birth_rate", "death_rate"))[[1]]  +
  # modify legend location using ggplot2
  theme(legend.position = c(0.80,0.80))

plot



##########################
# Ver EBD
##########################


# Especificar la ruta de los archivos de salida 
speciation_time_file <- "/home/daniela/GSM/Revbayes/Rev_out/98/EBD/98_2701_EBD_speciation_times.log"
speciation_rate_file <- "/home/daniela/GSM/Revbayes/Rev_out/98/EBD/98_2701_EBD_speciation_rates.log"
extinction_time_file <- "/home/daniela/GSM/Revbayes/Rev_out/98/EBD/98_2701_EBD_extinction_times.log"
extinction_rate_file <- "/home/daniela/GSM/Revbayes/Rev_out/98/EBD/98_2701_EBD_extinction_rates.log"

# Leer y procesar los datos
rates <- processDivRates(speciation_time_log = speciation_time_file,
                         speciation_rate_log = speciation_rate_file,
                         extinction_time_log = extinction_time_file,
                         extinction_rate_log = extinction_rate_file,
                         burnin = 0.25, # Se descarta el 25% inicial de las iteraciones como burn-in
                         summary = "median")  # Se toman las tasas medianas como resumen

# Graficar las tasas atravez del tiempo
p <- plotDivRates(rates = rates) +
  xlab("Millions of years ago") +
  ylab("Rate per million years")

p

