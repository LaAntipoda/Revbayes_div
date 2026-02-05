# Cargar las bibliotecas necesarias
library(RevGadgets)
library(ggplot2)
library(ggtree)

# nombres de archivo
# archivo de salida para la figura final en PDF
plot_fn = "/home/daniela/GSM/Revbayes/DEC/out/98_2601.range.pdf" 
# árbol con estados ancestrales (generado en RevBayes)
tree_fn = "/home/daniela/GSM/Revbayes/DEC/out/98_2601.ase.tre" 
# archivo con las etiquetas de los estados (rango binario → etiquetas)
label_fn = "/home/daniela/GSM/Revbayes/DEC/out/98_2601.state_labels.txt"
# archivo con los colores para cada rango
color_fn = "/home/daniela/GSM/Revbayes/DEC/data/range_colors.txt"

# obtener etiquetas de estado y colores de estado
states = make_states(label_fn, color_fn)
state_labels = states$state_labels
state_colors = states$state_colors

# procesar los estados ancestrales
ase <- processAncStates(tree_fn,
                        # Estas etiquetas numéricas deben coincidir con
                        # las del archivo original de datos de entrada.
                        state_labels = state_labels)



# colores


# Función para añadir líneas verticales que marcan los tiempos de formación de las islas
add_island_times = function(p, x_offset) {
  t = "dashed"
  p = p + geom_vline(xintercept=x_offset-5.05, color="#7570b3", linetype=t)
  p = p + geom_vline(xintercept=x_offset-5.15, color="#7570b3", linetype=t)
  p = p + geom_vline(xintercept=x_offset-2.2, color="#e7298a", linetype=t)
  p = p + geom_vline(xintercept=x_offset-3.7, color="#e7298a", linetype=t)
  p = p + geom_vline(xintercept=x_offset-1.3, color="#66a61e", linetype=t)
  p = p + geom_vline(xintercept=x_offset-1.8, color="#66a61e", linetype=t)
  p = p + geom_vline(xintercept=x_offset-0.3, color="#e6ab02", linetype=t)
  p = p + geom_vline(xintercept=x_offset-0.7, color="#e6ab02", linetype=t)
  return(p)
}

# Función para construir etiquetas y colores de los estados ancestrales
make_states = function(label_fn, color_fn, fp="./") {
  
  # Leer la tabla de colores para cada rango geográfico
  range_color_list = read.csv(color_fn, header=T, sep=",", colClasses="character")
  
  # Extraer los nombres de las áreas individuales (por ejemplo: "K", "O", "M", "H")
  area_names = unlist(sapply(range_color_list$range, function(y) { if (nchar(y)==1) { return(y) } }))
  
  # Leer el archivo con las etiquetas de los estados codificados (binarios)
  state_descriptions = read.csv(label_fn, header=T, sep=",", colClasses="character")
  
  # Convertir los rangos codificados como "1010" en nombres como "KM"
  range_labels = sapply(state_descriptions$range[2:nrow(state_descriptions)],
                        function(x) {
                          present = as.vector(gregexpr(pattern="1", x)[[1]])
                          paste( area_names[present], collapse="")
                        })
  
  # Asignar un color a cada rango
  range_colors = range_color_list$color[ match(range_labels, range_color_list$range) ]
  
  # Generar listas con las etiquetas y los colores por estado
  idx = 1
  st_lbl = list()
  st_colors = c()
  for (j in 1:(nrow(state_descriptions)-1)) {
    st_lbl[[ as.character(j) ]] = range_labels[j]
    st_colors[j] = range_colors[j]
  }
  
  # Añadir color para estados misceláneos o desconocidos
  st_colors[ length(st_colors)+1 ] = "lightgray"
  st_lbl[["misc."]] = "misc."
  
  # Devolver etiquetas y colores
  return( list(state_labels=st_lbl, state_colors=st_colors) )
}



############## Plot


# graficar los estados ancestrales como gráficos de pastel (pie charts)
pp  <- plotAncStatesPie(t = ase,
                        # incluir los estados del nodo raíz
                        include_start_states=T, 
                        # mostrar etiquetas de los estados
                        state_labels=state_labels,
                        # asignar los colores definidos
                        state_colors=state_colors, 
                        # tamaño de las etiquetas en los extremos (tips)
                        tip_label_size=2.5,
                        # distancia de desplazamiento de las etiquetas
                        tip_label_offset=0.1,
                        # sin etiquetas numéricas en los nodos
                        node_label_size=0, 
                        # sin etiquetas en los hombros
                        shoulder_label_size=0, 
                        # mostrar leyenda de probabilidades posteriores
                        show_posterior_legend=T, 
                        # tamaño de los círculos en los extremos (tips)
                        tip_pie_diameter=0.5, 
                        # tamaño de los círculos en los nodos
                        node_pie_diameter=2.0,
                        # ajuste horizontal de los gráficos de pastel
                        pie_nudge_x=0.03, 
                        # ajuste vertical de los gráficos de pastel
                        pie_nudge_y=0.16, 
                        alpha=1) + # opacidad
  
  # mover la leyenda de estados
  theme(legend.position = c(0.1, 0.75))


# obtener las dimensiones del gráfico
# obtener el valor máximo del eje x (altura del árbol)
x_phy = max(pp$data$x)
# espacio extra para las etiquetas de las hojas
x_label = 3.5 
# edad inicial del eje (debe ser mayor a la raíz del árbol)
x_start = 7  
# inicio del eje x
x0 = -(x_start - x_phy) 
# fin del eje x
x1 = x_phy + x_label 

# añadir eje cronológico
pp = pp + theme_tree2()
pp = pp + labs(x="Edad (Ma)")

# modificar límites del eje x
pp = pp + coord_cartesian(xlim=c(x0,x1), expand=TRUE)

# agregar marcas al eje y etiquetas de islas como eje secundario
island_axis = sec_axis(~ ., breaks=x_phy-c(5.1, 2.95, 1.55, 0.5), labels=c("+K","+O","+M","+H") )
x_breaks = seq(0,x_start,1) + x0
x_labels = rev(seq(0,x_start,1))
pp = pp + scale_x_continuous(breaks=x_breaks, labels=x_labels, sec.axis=island_axis)

# dibujar líneas verticales que marcan los tiempos de formación de islas
pp = add_island_times(pp, x_phy)

# configurar la posición de la leyenda
pp = pp + theme(legend.position="left")

# mostrar gráfico final
pp

# guardar el gráfico como archivo PDF
ggsave(file=plot_fn, plot=pp, device="pdf", height=7, width=10, useDingbats=F)

